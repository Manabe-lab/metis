"""
RNA Velocity Data Filtering
Filter loom files based on h5ad cell barcodes and merge multiple loom files
"""

import streamlit as st
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import h5py
import loompy
import os
import io
import tempfile
import time
import difflib
from scipy.sparse import csr_matrix
from helper_func import clear_old_directories, clear_old_files

st.set_page_config(page_title="RNA Velocity Data Filtering", page_icon="🚀", layout="wide")

st.title("🚀 RNA Velocity Data Filtering")
st.markdown("""
koofapuriwithis、velocytowithgenbecomesaretaloomFiletheh5adFileofCelltobaseduiteFilteringshi、
multinumofloomFilethema-jiwithkimasu。

### wa-kufuro-
1. **h5adFileLoading**: CellofmetaDatatheincludemuanndatashapeformatofFile
2. **IdentitySelect**: loomFiletopairresponddoSample/GroupofinfoinfotheholdtsumetaDatacoltheSelect
3. **LoomFileUpload**: velocytowithgenbecomesaretaloomFile（multinumcan）
4. **Groupratioricurrentte**: eachloomFileisdoofGrouptopairresponddokaSelect
5. **Filtering＆ma-ji**: barcodetobaseduiteCelltheFilteringshi、multinumFilethema-ji
6. **Download**: FilteringafterofloomFiletheDownload
""")

# Initialize session state
if "velocity_temp_dir" not in st.session_state:
    velocity_temp_dir = os.path.join("temp", f"velocity_{round(time.time())}")
    os.makedirs("temp", exist_ok=True)
    clear_old_directories("temp")
    clear_old_files("temp")
    os.makedirs(velocity_temp_dir, exist_ok=True)
    st.session_state.velocity_temp_dir = velocity_temp_dir
else:
    velocity_temp_dir = st.session_state.velocity_temp_dir

if "filtered_loom_ready" not in st.session_state:
    st.session_state.filtered_loom_ready = False

@st.cache_data
def read_h5ad(file):
    """Read h5ad file"""
    adata = sc.read_h5ad(file)
    return adata

def extract_barcode(cell_id):
    """
    Extract barcode sequence from cell ID
    Extracts the actual DNA sequence part (ACGT characters) followed by optional suffix

    Examples:
    - "CTRL1:AAACCCAAGAAACCAT-1" -> "AAACCCAAGAAACCAT"
    - "MSC_control1_AAACCCAAGAAACCAT-1" -> "AAACCCAAGAAACCAT"
    - "possorted_genome_bam_L7XS7:AAACCCAAGAAACCAT-1x" -> "AAACCCAAGAAACCAT"
    """
    import re

    # Remove prefix before colon if present (velocyto format: sample:barcode)
    if ':' in cell_id:
        cell_id = cell_id.split(':')[1]

    # Remove any prefix before the last underscore if present
    if '_' in cell_id:
        parts = cell_id.split('_')
        cell_id = parts[-1]

    # Extract pure barcode sequence (DNA letters + optional dash + digit(s) + optional suffix)
    # Pattern: consecutive ACGT letters, optionally followed by -digit(s), ignore any suffix
    match = re.search(r'([ACGT]+)(?:-\d+)?', cell_id, re.IGNORECASE)

    if match:
        barcode = match.group(1).upper()
        return barcode

    # Fallback: just return the cleaned cell_id
    # Remove common suffixes like -1, -1x, x, etc.
    barcode = re.sub(r'[-x]\d*[x]?$', '', cell_id, flags=re.IGNORECASE)
    return barcode.upper()

def match_barcodes(loom_cells, h5ad_cells):
    """
    Match loom cell IDs to h5ad cell IDs based on barcodes
    Returns: matched indices in loom file, corresponding h5ad cell names
    """
    # Extract barcodes from both
    loom_barcodes = [extract_barcode(cell) for cell in loom_cells]
    h5ad_barcodes = [extract_barcode(cell) for cell in h5ad_cells]

    # Create mapping from barcode to h5ad cell name
    h5ad_barcode_to_name = dict(zip(h5ad_barcodes, h5ad_cells))

    # Find matches
    matched_loom_indices = []
    matched_h5ad_names = []

    for i, loom_bc in enumerate(loom_barcodes):
        if loom_bc in h5ad_barcode_to_name:
            matched_loom_indices.append(i)
            matched_h5ad_names.append(h5ad_barcode_to_name[loom_bc])

    return matched_loom_indices, matched_h5ad_names

def filter_loom_by_barcodes(loom_file, matched_indices, new_cell_names, output_path):
    """
    Filter loom file to keep only matched cells and rename them
    """
    with loompy.connect(loom_file, mode='r') as ds_in:
        # Get dimensions
        n_genes = ds_in.shape[0]
        n_cells_original = ds_in.shape[1]
        n_cells_filtered = len(matched_indices)

        st.info(f"Original cells: {n_cells_original}, Filtered cells: {n_cells_filtered}")

        if n_cells_filtered == 0:
            st.error("No matching cells found!")
            return None

        # Read layers for matched cells
        layers = {}
        for layer_name in ds_in.layers.keys():
            layers[layer_name] = ds_in.layers[layer_name][:, matched_indices]

        # Read row attributes (genes)
        row_attrs = {key: ds_in.ra[key] for key in ds_in.ra.keys()}

        # Create column attributes (cells) with new names
        col_attrs = {
            "CellID": np.array(new_cell_names, dtype=object)
        }

        # Copy other column attributes if they exist
        for key in ds_in.ca.keys():
            if key != "CellID":
                col_attrs[key] = ds_in.ca[key][matched_indices]

        # Create new loom file
        loompy.create(
            output_path,
            layers=layers,
            row_attrs=row_attrs,
            col_attrs=col_attrs
        )

    return output_path

def merge_loom_files(loom_files, output_path):
    """
    Merge multiple loom files by taking the union of genes across all files
    Missing genes in each file are filled with zeros
    """
    import shutil
    from scipy.sparse import csr_matrix, hstack

    if len(loom_files) == 1:
        # No merging needed, just copy
        shutil.copy(loom_files[0], output_path)
        return output_path

    st.write(f"Merging {len(loom_files)} loom files...")

    # Step 1: Collect all unique genes (union)
    all_genes = set()
    gene_list_per_file = []

    for loom_file in loom_files:
        with loompy.connect(loom_file, mode='r') as ds:
            genes = ds.ra.Gene
            gene_list_per_file.append(genes)
            all_genes.update(genes)

    # Create sorted list of all unique genes
    all_genes_sorted = sorted(list(all_genes))
    n_genes = len(all_genes_sorted)

    st.info(f"Total unique genes across all files: {n_genes}")

    # Create gene index mapping
    gene_to_idx = {gene: idx for idx, gene in enumerate(all_genes_sorted)}

    # Step 2: Read each file and align to the union gene set
    all_layers = {}  # layer_name -> list of matrices
    all_cell_ids = []
    all_col_attrs = {}
    first_file = True

    for file_idx, loom_file in enumerate(loom_files):
        with loompy.connect(loom_file, mode='r') as ds:
            file_genes = ds.ra.Gene
            n_cells_this_file = ds.shape[1]

            st.write(f"File {file_idx+1}: {n_cells_this_file} cells, {len(file_genes)} genes")

            # Create mapping from this file's genes to union gene indices
            gene_indices = np.array([gene_to_idx[g] for g in file_genes])

            # Process each layer
            for layer_name in ds.layers.keys():
                if first_file:
                    all_layers[layer_name] = []

                # Get layer data
                layer_data = ds.layers[layer_name][:, :]  # genes x cells

                # Create aligned matrix (fill with zeros for missing genes)
                aligned_matrix = np.zeros((n_genes, n_cells_this_file), dtype=layer_data.dtype)
                aligned_matrix[gene_indices, :] = layer_data

                all_layers[layer_name].append(aligned_matrix)

            # Collect cell IDs
            all_cell_ids.extend(ds.ca.CellID)

            # Collect other column attributes
            for key in ds.ca.keys():
                if key != "CellID":
                    if first_file:
                        all_col_attrs[key] = []
                    all_col_attrs[key].extend(ds.ca[key])

            first_file = False

    # Step 3: Concatenate all matrices horizontally
    merged_layers = {}
    for layer_name, matrices in all_layers.items():
        merged_layers[layer_name] = np.concatenate(matrices, axis=1)
        st.write(f"Layer '{layer_name}': {merged_layers[layer_name].shape}")

    # Step 4: Create row attributes (genes)
    # We need to get gene attributes from the first occurrence of each gene
    row_attrs = {'Gene': np.array(all_genes_sorted, dtype=object)}

    # Try to get other gene attributes from first file
    with loompy.connect(loom_files[0], mode='r') as ds:
        for key in ds.ra.keys():
            if key != 'Gene':
                # Create default values for all genes
                if ds.ra[key].dtype.kind in ['U', 'S', 'O']:  # String types
                    default_val = ''
                else:
                    default_val = 0

                attr_values = np.full(n_genes, default_val, dtype=ds.ra[key].dtype)

                # Fill in known values from each file
                for loom_file in loom_files:
                    with loompy.connect(loom_file, mode='r') as ds2:
                        if key in ds2.ra.keys():
                            file_genes = ds2.ra.Gene
                            for i, gene in enumerate(file_genes):
                                gene_idx = gene_to_idx[gene]
                                attr_values[gene_idx] = ds2.ra[key][i]

                row_attrs[key] = attr_values

    # Step 5: Create column attributes (cells)
    col_attrs = {
        'CellID': np.array(all_cell_ids, dtype=object)
    }
    for key, values in all_col_attrs.items():
        col_attrs[key] = np.array(values)

    # Step 6: Create merged loom file
    # Remove output file if it already exists
    if os.path.exists(output_path):
        os.remove(output_path)

    loompy.create(
        output_path,
        layers=merged_layers,
        row_attrs=row_attrs,
        col_attrs=col_attrs
    )

    st.success(f"✓ Merged {len(loom_files)} loom files successfully!")
    st.write(f"Final dimensions: {n_genes} genes × {len(all_cell_ids)} cells")

    return output_path

def find_best_matching_group(loom_filename, group_list):
    """
    Find the best matching group for a loom file based on string similarity

    Parameters:
    -----------
    loom_filename : str
        Name of the loom file (e.g., 'sample1_velocyto.loom')
    group_list : list
        List of available group names

    Returns:
    --------
    str : Best matching group name
    """
    # Remove file extension and common suffixes
    base_name = os.path.splitext(loom_filename)[0]

    # Remove common velocyto-related suffixes
    for suffix in ['_velocyto', '_loom', 'velocyto', '_original']:
        if base_name.endswith(suffix):
            base_name = base_name[:-len(suffix)]

    # Convert to lowercase for comparison
    base_name_lower = base_name.lower()

    # Strategy 1: Exact substring match (case-insensitive)
    for group in group_list:
        if group.lower() in base_name_lower or base_name_lower in group.lower():
            return group

    # Strategy 2: Use difflib to find close matches
    # Calculate similarity ratio for each group
    similarity_scores = []
    for group in group_list:
        ratio = difflib.SequenceMatcher(None, base_name_lower, group.lower()).ratio()
        similarity_scores.append((group, ratio))

    # Sort by similarity score (highest first)
    similarity_scores.sort(key=lambda x: x[1], reverse=True)

    # If best match has reasonable similarity (>0.3), use it
    best_match, best_score = similarity_scores[0]

    if best_score > 0.3:
        return best_match

    # Strategy 3: Check if group name is part of any word in filename
    base_name_words = base_name_lower.replace('_', ' ').replace('-', ' ').split()
    for group in group_list:
        group_lower = group.lower()
        for word in base_name_words:
            if group_lower in word or word in group_lower:
                return group

    # Default: return first group
    return group_list[0]

# ========================================
# Step 1: Upload h5ad file
# ========================================
st.header("Step 1: Upload h5ad file")
uploaded_h5ad = st.file_uploader("Upload h5ad file", type=['h5ad'], key="h5ad_upload")

if uploaded_h5ad is not None:
    adata = read_h5ad(uploaded_h5ad)

    # Store h5ad filename for later use
    h5ad_basename = os.path.splitext(uploaded_h5ad.name)[0]
    st.session_state.h5ad_basename = h5ad_basename

    st.success(f"✓ Loaded h5ad: {adata.n_obs} cells, {adata.n_vars} genes")

    # Display metadata
    st.subheader("Metadata preview")
    st.dataframe(adata.obs.head())

    # Get metadata columns
    metadata_cols = adata.obs.columns.tolist()

    # ========================================
    # Step 2: Select identity column or use all cells
    # ========================================
    st.header("Step 2: Select matching mode")

    # Option to use all cells without grouping
    use_all_cells = st.checkbox(
        "**Use all cells (no grouping)**",
        value=False,
        help="allofCelltheone包withmachingushimasu。multinumofloomFileisexistplacematch、allthesameh5adCellsetotopairshitemachingushimasu。"
    )

    if use_all_cells:
        st.info("ℹ️ allofh5adCellismachingupairobjtonarimasu。Groupdivratioisrowimasen。")

        # Set up for all cells mode
        st.session_state.use_all_cells = True
        st.session_state.identity_col = None
        st.session_state.group_cells = {"all": adata.obs_names.tolist()}
        st.session_state.adata_cells = adata.obs_names.tolist()

        unique_groups = ["all"]
        st.success(f"✓ Using all {adata.n_obs} cells for matching")

    else:
        st.session_state.use_all_cells = False

        with st.form("identity_form"):
            st.markdown("""
            loomFiletopairresponddoSampleyaGroupofinfoinfotheholdtsumetaDatacoltheSelectshitekudasai。
            Example: `sample`, `patient`, `condition` etc
            """)

            # Set default to 'orig.ident' if it exists
            default_index = 0
            if 'orig.ident' in metadata_cols:
                default_index = metadata_cols.index('orig.ident')

            identity_col = st.selectbox(
                "Identity column:",
                metadata_cols,
                index=default_index,
                help="eachloomFiletopairresponddoGroupinfoinfotheholdtsucol"
            )

            submit_identity = st.form_submit_button("Confirm identity column")

        if submit_identity or ("identity_col" in st.session_state and st.session_state.identity_col is not None):
            if submit_identity:
                st.session_state.identity_col = identity_col
            else:
                identity_col = st.session_state.identity_col

            st.success(f"✓ Selected identity: **{identity_col}**")

            # Show unique values
            unique_groups = adata.obs[identity_col].unique()
            st.write(f"Unique groups in {identity_col}:")
            st.write(unique_groups)

            # Get cell IDs for each group
            group_cells = {}
            for group in unique_groups:
                group_cells[group] = adata.obs_names[adata.obs[identity_col] == group].tolist()

            st.session_state.group_cells = group_cells
            st.session_state.adata_cells = adata.obs_names.tolist()
        else:
            unique_groups = None

    # Continue only if we have groups defined
    if use_all_cells or (not use_all_cells and "identity_col" in st.session_state and st.session_state.identity_col is not None):
        if not use_all_cells:
            unique_groups = adata.obs[st.session_state.identity_col].unique()

        # ========================================
        # Step 3: Upload loom files
        # ========================================
        st.header("Step 3: Upload loom files")

        uploaded_looms = st.file_uploader(
            "Upload loom file(s)",
            type=['loom'],
            accept_multiple_files=True,
            key="loom_upload",
            help="velocytowithgenbecomesaretaloomFilethe1tsuorupUpload"
        )

        if uploaded_looms:
            st.success(f"✓ Uploaded {len(uploaded_looms)} loom file(s)")

            # ========================================
            # Step 4: Assign groups to loom files (skip if using all cells)
            # ========================================
            if use_all_cells:
                # In "all cells" mode, assign all loom files to "all" group
                st.header("Step 4: Group assignment (skipped)")
                st.info("ℹ️ allCellmo-dooffor、Groupratioricurrenttethesukipushimasu。allofloomFileissameh5adCellsetotopairshitemachingusaremasu。")

                loom_to_group = {loom_file.name: "all" for loom_file in uploaded_looms}
                st.session_state.loom_to_group = loom_to_group

                # Show files to be processed
                st.write("**procprocpairobjFile:**")
                for loom_name in loom_to_group.keys():
                    st.write(f"- {loom_name}")

            else:
                st.header("Step 4: Assign groups to loom files")

                with st.form("loom_assignment_form"):
                    st.markdown("""
                    eachloomFileisdoofGrouptopairresponddokaSelectshitekudasai。
                    """)

                    loom_to_group = {}
                    for loom_file in sorted(uploaded_looms, key=lambda f: f.name):
                        # Find best matching group
                        default_group = find_best_matching_group(loom_file.name, unique_groups)
                        default_index = list(unique_groups).index(default_group)

                        selected_group = st.selectbox(
                            f"{loom_file.name} corresponds to:",
                            unique_groups,
                            index=default_index,
                            key=f"group_{loom_file.name}"
                        )
                        loom_to_group[loom_file.name] = selected_group

                    submit_assignment = st.form_submit_button("Confirm assignments")

                if submit_assignment:
                    st.session_state.loom_to_group = loom_to_group
                    st.success("✓ Group assignments confirmed")

                    # Show assignments
                    st.write("**Assignments:**")
                    for loom_name, group in loom_to_group.items():
                        st.write(f"- {loom_name} → {group}")

            # ========================================
            # Step 5: Process files
            # ========================================
            if "loom_to_group" in st.session_state:
                st.header("Step 5: Filter and merge")

                if st.button("🚀 Process files", type="primary"):
                    with st.spinner("Processing loom files..."):
                        filtered_loom_files = []

                        progress_bar = st.progress(0)

                        for idx, loom_file in enumerate(uploaded_looms):
                            st.subheader(f"Processing: {loom_file.name}")

                            # Get assigned group
                            assigned_group = st.session_state.loom_to_group[loom_file.name]

                            # Get h5ad cells for this group
                            group_h5ad_cells = st.session_state.group_cells[assigned_group]

                            st.info(f"Group: {assigned_group}, h5ad cells: {len(group_h5ad_cells)}")

                            # Save uploaded loom file temporarily
                            temp_loom_path = os.path.join(velocity_temp_dir, f"temp_{loom_file.name}")
                            with open(temp_loom_path, "wb") as f:
                                f.write(loom_file.read())

                            # Read loom cell IDs
                            with loompy.connect(temp_loom_path, mode='r') as ds:
                                loom_cells = ds.ca.CellID.tolist()
                                st.write(f"Loom cells: {len(loom_cells)}")

                            # Show sample cell IDs for debugging
                            with st.expander("🔍 Debug: View sample cell IDs"):
                                col1, col2 = st.columns(2)
                                with col1:
                                    st.write("**Loom cells (first 5):**")
                                    for i, cell in enumerate(loom_cells[:5]):
                                        barcode = extract_barcode(cell)
                                        st.code(f"{cell}\n→ {barcode}")
                                with col2:
                                    st.write("**h5ad cells (first 5):**")
                                    for i, cell in enumerate(group_h5ad_cells[:5]):
                                        barcode = extract_barcode(cell)
                                        st.code(f"{cell}\n→ {barcode}")

                            # Match barcodes
                            matched_indices, matched_names = match_barcodes(loom_cells, group_h5ad_cells)

                            st.write(f"Matched cells: {len(matched_indices)}")

                            if len(matched_indices) == 0:
                                st.warning(f"⚠ No matching cells found in {loom_file.name}")
                                continue

                            # Filter loom file
                            output_loom = os.path.join(
                                velocity_temp_dir,
                                f"filtered_{assigned_group}_{loom_file.name}"
                            )

                            filter_loom_by_barcodes(
                                temp_loom_path,
                                matched_indices,
                                matched_names,
                                output_loom
                            )

                            filtered_loom_files.append(output_loom)

                            # Clean up temp file
                            os.remove(temp_loom_path)

                            # Update progress
                            progress_bar.progress((idx + 1) / len(uploaded_looms))

                        progress_bar.empty()

                        if len(filtered_loom_files) == 0:
                            st.error("❌ No cells matched. Please check your data.")
                        else:
                            # Merge if multiple files
                            st.subheader("Merging files...")

                            final_loom_path = os.path.join(
                                velocity_temp_dir,
                                "filtered_merged.loom"
                            )

                            merge_loom_files(filtered_loom_files, final_loom_path)

                            # Get final statistics
                            with loompy.connect(final_loom_path, mode='r') as ds:
                                n_genes = ds.shape[0]
                                n_cells = ds.shape[1]

                            st.success(f"""
                            ✓ Processing complete!

                            **Final loom file:**
                            - Genes: {n_genes}
                            - Cells: {n_cells}
                            """)

                            st.session_state.final_loom_path = final_loom_path
                            st.session_state.filtered_loom_ready = True

                # ========================================
                # Step 6: Download
                # ========================================
                if st.session_state.filtered_loom_ready:
                    st.header("Step 6: Download filtered loom file")

                    with open(st.session_state.final_loom_path, "rb") as f:
                        loom_bytes = f.read()

                    # Create filename based on h5ad filename
                    if "h5ad_basename" in st.session_state:
                        output_filename = f"{st.session_state.h5ad_basename}.loom"
                    else:
                        output_filename = "filtered_velocity.loom"

                    st.download_button(
                        label="⬇️ Download filtered loom file",
                        data=loom_bytes,
                        file_name=output_filename,
                        mime="application/octet-stream",
                        type="primary"
                    )

                    st.info("""
                    ### nextofsutepu

                    koofFilteringsaretaloomFileandh5adFiletheuseuseshite、
                    RNA velocityAnalysistheRunwithkimasu（Example: scVelo）。

                    PythonwithofuseuseExample:
                    ```python
                    import scvelo as scv
                    import scanpy as sc

                    # Load files
                    adata = sc.read_h5ad('your_file.h5ad')
                    ldata = scv.read('filtered_velocity.loom')

                    # Merge
                    adata = scv.utils.merge(adata, ldata)

                    # Run velocity analysis
                    scv.pp.filter_and_normalize(adata)
                    scv.pp.moments(adata)
                    scv.tl.velocity(adata)
                    scv.tl.velocity_graph(adata)
                    scv.pl.velocity_embedding_stream(adata, basis='umap')
                    ```
                    """)

else:
    st.info("👆 h5adFiletheUploadshitestartshitekudasai")
