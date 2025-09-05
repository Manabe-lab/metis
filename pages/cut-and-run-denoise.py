import streamlit as st
import pandas as pd
import numpy as np
import gzip
from scipy.signal import find_peaks, savgol_filter
from scipy.ndimage import gaussian_filter1d
from numba import jit, prange
import io
import base64
import tempfile
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

st.set_page_config(page_title="CUT&RUN Denoising App", layout="wide")

st.title("CUT&RUN Bedgraph Denoising App (Parallel Processing)")

# File upload
uploaded_files = st.file_uploader("Upload your bedgraph files", accept_multiple_files=True, type=['bedgraph', 'bg', 'gz'])

# Algorithm selection
algorithm = st.selectbox("Select denoising algorithm", ["Enhanced Smoothing", "Savitzky-Golay Filter"])

# Common parameters
chunk_size = st.slider("Chunk size (rows)", 10000, 1000000, 100000, step=10000)
n_processes = st.slider("Number of processes", 1, multiprocessing.cpu_count(), min(4, multiprocessing.cpu_count()))

# Algorithm-specific parameters
if algorithm == "Enhanced Smoothing":
    peak_height_percentile = st.slider("Peak height percentile", 80, 99, 95)
    peak_width = st.slider("Peak width", 5, 20, 10)
    savgol_window = st.slider("Savitzky-Golay window", 3, 21, 11, step=2)
    savgol_order = st.slider("Savitzky-Golay order", 1, 5, 3)
    gaussian_sigma = st.slider("Gaussian sigma", 0.5, 3.0, 1.5, 0.1)
    smoothing_rounds = st.slider("Smoothing rounds", 1, 5, 2)
else:
    window_length = st.slider("Window length", 11, 101, 51, step=2)
    polyorder = st.slider("Polynomial order", 1, 5, 3)

# Helper functions
def read_bedgraph_chunk(file, chunk_size):
    open_func = gzip.open if file.name.endswith('.gz') else open
    with open_func(file, 'rt') as f:
        while True:
            chunk = pd.read_csv(f, sep='\t', header=None, names=['chrom', 'start', 'end', 'value'],
                                chunksize=chunk_size, iterator=True, nrows=chunk_size)
            try:
                yield next(chunk)
            except StopIteration:
                break

# Enhanced Smoothing algorithm
def smooth_and_preserve_peaks(values, peak_height_percentile, peak_width, 
                              savgol_window, savgol_order, 
                              gaussian_sigma, smoothing_rounds):
    peak_height = np.percentile(values, peak_height_percentile)
    peaks, _ = find_peaks(values, height=peak_height, distance=peak_width)
    
    peak_mask = np.zeros_like(values, dtype=bool)
    for peak in peaks:
        start = max(0, peak - peak_width // 2)
        end = min(len(values), peak + peak_width // 2)
        peak_mask[start:end] = True
    
    smoothed = np.copy(values)
    for _ in range(smoothing_rounds):
        smoothed_savgol = savgol_filter(smoothed, savgol_window, savgol_order)
        smoothed_gauss = gaussian_filter1d(smoothed_savgol, gaussian_sigma)
        smoothed[~peak_mask] = smoothed_gauss[~peak_mask]
    
    return smoothed

# Savitzky-Golay Filter algorithm
@jit(nopython=True)
def savgol_filter_numba(y, window_length, polyorder):
    half_window = (window_length - 1) // 2
    x = np.arange(-half_window, half_window+1, dtype=np.float64)
    order = np.arange(polyorder + 1).reshape(-1, 1)
    A = x ** order
    ATA = A.T @ A
    B = np.zeros(polyorder + 1)
    B[0] = 1
    coeffs = np.linalg.solve(ATA, B) @ A
    return np.convolve(y, coeffs, mode='same')

@jit(nopython=True)
def interpolate_and_denoise(starts, ends, values, all_positions, window_length, polyorder):
    all_values = np.zeros_like(all_positions, dtype=np.float64)
    for i in prange(len(starts)):
        mask = (all_positions >= starts[i]) & (all_positions < ends[i])
        all_values[mask] = values[i]
    return savgol_filter_numba(all_values, window_length, polyorder)

# Processing function
def process_chunk(chunk, algorithm, params):
    result_dfs = []

    for chrom in chunk['chrom'].unique():
        chrom_data = chunk[chunk['chrom'] == chrom]
        if chrom_data.empty:
            continue

        if algorithm == "Enhanced Smoothing":
            values = chrom_data['value'].values
            smoothed_values = smooth_and_preserve_peaks(values, **params)
            result = pd.DataFrame({
                'chrom': chrom_data['chrom'],
                'start': chrom_data['start'],
                'end': chrom_data['end'],
                'value': smoothed_values
            })
        else:  # Savitzky-Golay Filter
            all_positions = np.arange(chrom_data['start'].min(), chrom_data['end'].max())
            denoised_data = interpolate_and_denoise(
                chrom_data['start'].values, chrom_data['end'].values, chrom_data['value'].values,
                all_positions, params['window_length'], params['polyorder']
            )
            result = pd.DataFrame({
                'chrom': chrom,
                'start': all_positions[:-1],
                'end': all_positions[1:],
                'value': denoised_data[:-1]
            })

        result = result[(result['start'] < result['end']) & 
                        (~result['value'].isna()) & 
                        (~np.isinf(result['value']))]
        result_dfs.append(result)

    return pd.concat(result_dfs, ignore_index=True) if result_dfs else pd.DataFrame()

# Main app logic
if st.button("Process Files"):
    if not uploaded_files:
        st.error("Please upload at least one file.")
    else:
        for uploaded_file in uploaded_files:
            st.write(f"Processing {uploaded_file.name}...")
            
            if algorithm == "Enhanced Smoothing":
                params = {
                    'peak_height_percentile': peak_height_percentile,
                    'peak_width': peak_width,
                    'savgol_window': savgol_window,
                    'savgol_order': savgol_order,
                    'gaussian_sigma': gaussian_sigma,
                    'smoothing_rounds': smoothing_rounds
                }
            else:
                params = {
                    'window_length': window_length,
                    'polyorder': polyorder
                }
            
            # Create a temporary file to store the processed results
            with tempfile.NamedTemporaryFile(mode='w+t', delete=False, suffix='.bedgraph') as temp_file:
                # Parallel processing of chunks
                with ProcessPoolExecutor(max_workers=n_processes) as executor:
                    futures = []
                    for chunk in read_bedgraph_chunk(uploaded_file, chunk_size):
                        future = executor.submit(process_chunk, chunk, algorithm, params)
                        futures.append(future)
                    
                    for future in as_completed(futures):
                        result_df = future.result()
                        result_df.to_csv(temp_file, sep='\t', index=False, header=False)
                
                temp_file_path = temp_file.name

            # Create a download link for the processed file
            with open(temp_file_path, 'rb') as f:
                bytes_data = f.read()
            b64 = base64.b64encode(bytes_data).decode()
            href = f'<a href="data:file/txt;base64,{b64}" download="{uploaded_file.name}_processed.bedgraph">Download Processed File</a>'
            st.markdown(href, unsafe_allow_html=True)
            
            # Clean up the temporary file
            os.unlink(temp_file_path)
            
            st.success(f"Processed {uploaded_file.name}")

st.info("This app processes CUT&RUN bedgraph files to denoise the data using parallel processing. You can choose between two algorithms: Enhanced Smoothing and Savitzky-Golay Filter. Upload your files, select the desired algorithm and parameters, then click 'Process Files' to start.")