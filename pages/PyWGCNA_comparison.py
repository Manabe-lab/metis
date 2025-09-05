import streamlit as st
import PyWGCNA
import os
import tempfile
import zipfile
import shutil
from helper_func import clear_old_directories
from helper_func import clear_old_files
import time
import sys
import pymupdf
from PIL import Image
import anndata


@st.cache_data
def read_pywgcna(file):
    return PyWGCNA.readWGCNA(file)

@st.cache_data
def compare_networks(_pywgcnas):
    return PyWGCNA.compareNetworks(PyWGCNAs=_pywgcnas)

import matplotlib.colors as mcolors

@st.cache_data
def get_light_colors():
    # matplotlibの名前付き色から"light"で始まる色を抽出
    light_colors = [color for color in mcolors.CSS4_COLORS if color.startswith('light')]
    return light_colors


def assign_colors_to_objects(objects):
    light_colors = get_light_colors()
    color_dict = {}
    
    for i, obj in enumerate(objects):
        # 利用可能な色の数で循環させる
        color = light_colors[i % len(light_colors)]
        color_dict[obj.name] = color  
    return color_dict

@st.cache_data
def bubble(bubble_file):
    comparison.plotBubbleComparison(color=color, plot_format="pdf", file_name=bubble_file)



def convert_pdf_to_images(pdf_path):
    pdf = pymupdf.open(pdf_path)
    page = pdf[0]  # 最初のページを取得
    pix = page.get_pixmap()
    img = Image.frombytes("RGB", [pix.width, pix.height], pix.samples)
    return img

def main():
    st.title("PyWGCNA Objects Comparison")

    # temp内に保存する
    # --- Initialising SessionState ---
    if "temp_dir" not in st.session_state:
        st.session_state.temp_dir = None
        #古いdirecotryとファイルを削除する
        temp_dir = "temp/" + str(round(time.time()))
        if not os.path.exists('temp'):
            os.mkdir('temp')
        else:
            clear_old_directories("temp")
            clear_old_files("temp")
        os.mkdir(temp_dir)
        st.session_state.temp_dir = temp_dir
        res_dir = temp_dir + '/figures'
        st.session_state.res_dir = res_dir
        os.mkdir(res_dir)

    else:
        temp_dir = st.session_state.temp_dir
        res_dir = temp_dir + '/figures'
        st.session_state.res_dir = res_dir
        if not os.path.exists(temp_dir):
            os.mkdir(temp_dir)
            os.mkdir(res_dir)
        if not os.path.exists(res_dir):
            os.mkdir(res_dir)

    # File uploader
    uploaded_files = st.file_uploader("Upload PyWGCNA objects", type=["p"], accept_multiple_files=True)

    if uploaded_files is not None:
        if st.button('Run'):
            if len(uploaded_files) >= 2:
                pywgcnas = []
                file_list = []
                for file in uploaded_files:
                    temp_file_path = os.path.join(temp_dir, file.name)
                    with open(temp_file_path, "wb") as f:
                        f.write(file.getbuffer())
                    pywgcna_obj = read_pywgcna(temp_file_path)
                    pywgcnas.append(pywgcna_obj)
                    file_list.append(temp_file_path)

                comparison = compare_networks(pywgcnas)

                # Generate base filename
                base_filename = "_".join([obj.name for obj in pywgcnas])

                st.write('jaccard_similarity:')
                st.write(comparison.jaccard_similarity.head(5))


                # Jaccard Similarity Plot
                with st.sidebar:
                    cutoff = st.number_input("#### Jaccard Similarity Cutoff", value = 0.2)
                st.subheader("Jaccard Similarity Plot")

                # WGCNAオブジェクトのリストに対して色を割り当てる

                color = assign_colors_to_objects(pywgcnas)

                jaccard_file = os.path.join(temp_dir, f"jaccard_similarity_{base_filename}")
                comparison.plotJaccardSimilarity(color=color, cutoff=cutoff, plot_format="pdf", file_name=jaccard_file)

                try:
                    # PDFを画像に変換
                    img = convert_pdf_to_images(jaccard_file + '.pdf')
                    # 画像を表示
                    st.image(img, use_column_width=True)
                except Exception as e:
                    st.error(f"エラーが発生しました: {str(e)}")
                    st.write("Probably nothing pass the threshold")  

                # Heatmap Comparison
                st.subheader("Heatmap Comparison")
                heatmap_file = os.path.join(temp_dir, f"heatmap_comparison_{base_filename}")
                comparison.plotHeatmapComparison(plot_format="pdf", file_name=heatmap_file)

                try:
                    # PDFを画像に変換
                    img = convert_pdf_to_images(heatmap_file + '.pdf')
                    # 画像を表示
                    st.image(img, use_column_width=True)
                except Exception as e:
                    st.error(f"エラーが発生しました: {str(e)}")    

                # Bubble Comparison
                st.subheader("Bubble Comparison")
                bubble_file = os.path.join(temp_dir, f"bubble_comparison_{base_filename}")
                comparison.plotBubbleComparison(color=color, plot_format="pdf", file_name=bubble_file)

                try:
                    # PDFを画像に変換
                    img = convert_pdf_to_images(bubble_file + '.pdf')
                    # 画像を表示
                    st.image(img, use_column_width=True)
                except Exception as e:
                    st.error(f"エラーが発生しました: {str(e)}")    

                # Save Comparison
                comparison_file = os.path.join(temp_dir, f"comparison_{base_filename}")
                comparison.saveComparison(name=comparison_file)
                st.success("Comparison saved")

                for i in file_list:
                    os.remove(i)

                # Create ZIP file
                zip_path = os.path.join(temp_dir, "results.zip")
                with zipfile.ZipFile(zip_path, 'w') as zipf:
                    for root, dirs, files in os.walk(temp_dir):
                        for file in files:
                            if file != "results.zip":
                                zipf.write(os.path.join(root, file), file)

                # Download button
                with open(zip_path, "rb") as f:
                    st.download_button(
                        label="Download Results",
                        data=f,
                        file_name=f"pywgcna_comparison_results_{base_filename}.zip",
                        mime="application/zip"
                    )

            # Clean up temporary directory
            shutil.rmtree(temp_dir)

if __name__ == "__main__":
    main()