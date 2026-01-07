"""
Add Metadata to h5ad File
Merge metadata from another h5ad file or TSV file into a target h5ad file
"""

import streamlit as st
import scanpy as sc
import pandas as pd
import numpy as np
import os
import tempfile
from io import BytesIO

@st.cache_data
def load_h5ad_file(file_bytes, file_name):
    """
    Load h5ad file from uploaded bytes with caching

    Parameters
    ----------
    file_bytes : bytes
        File content in bytes
    file_name : str
        Original filename (used as cache key)

    Returns
    -------
    adata : AnnData
        Loaded AnnData object
    """
    with tempfile.NamedTemporaryFile(delete=False, suffix='.h5ad') as tmp:
        tmp.write(file_bytes)
        tmp_path = tmp.name

    try:
        adata = sc.read_h5ad(tmp_path)
        return adata
    finally:
        os.unlink(tmp_path)

st.set_page_config(page_title="Add Metadata to h5ad", page_icon="📝", layout="wide")

st.title("📝 Add Metadata to h5ad File")
st.markdown("""
koofapuriis、h5adFiletosepofh5adFilealsoisTSVFilefrommetaDatatheaddaddshimasu。

### useiway
1. **meinh5adFile**theUpload（metaDatatheaddadddopairobj）
2. **metaDataso-su**theUpload（h5adFilealsoisTSVFile）
3. addaddshitai**metaDatacol**theSelect
4. Cellname（cell barcode）withmachingushitemetaDatatheaddadd
5. modcorrectshitah5adFiletheDownload（sourceofFilename.mod.h5ad）
""")

st.markdown("---")

# Step 1: Upload main h5ad file
st.header("Step 1: Upload main h5ad file")
main_h5ad = st.file_uploader(
    "Upload main h5ad file (metaDatatheaddadddopairobj)",
    type=['h5ad'],
    key='main_h5ad',
    help="metaDatatheaddaddshitaih5adFiletheUpload"
)

if main_h5ad is not None:
    # Load main h5ad with caching
    with st.spinner("Loading main h5ad file..."):
        adata_main = load_h5ad_file(main_h5ad.getvalue(), main_h5ad.name)
        st.success(f"✓ Loaded main h5ad: {adata_main.n_obs} cells, {adata_main.n_vars} genes")

    # Show current metadata columns
    with st.expander("📊 Current metadata columns", expanded=False):
        st.write(f"**Total columns**: {len(adata_main.obs.columns)}")
        st.dataframe(adata_main.obs.head())

    # Show cell barcodes (first few)
    with st.expander("🔍 Cell barcodes (first 10)", expanded=False):
        st.write(adata_main.obs.index[:10].tolist())

    st.markdown("---")

    # Step 2: Upload metadata source
    st.header("Step 2: Upload metadata source")

    metadata_source_type = st.radio(
        "Metadata source type",
        ["h5ad file", "TSV file"],
        help="metaDataofso-sutaiputheSelect"
    )

    if metadata_source_type == "h5ad file":
        metadata_file = st.file_uploader(
            "Upload h5ad file (metadata source)",
            type=['h5ad'],
            key='metadata_h5ad',
            help="metaDatathegetgetdoh5adFile"
        )

        if metadata_file is not None:
            with st.spinner("Loading metadata h5ad file..."):
                adata_metadata = load_h5ad_file(metadata_file.getvalue(), metadata_file.name)
                metadata_df = adata_metadata.obs.copy()
                st.success(f"✓ Loaded metadata h5ad: {adata_metadata.n_obs} cells, {len(metadata_df.columns)} metadata columns")

            # Show metadata
            with st.expander("📊 Available metadata columns", expanded=True):
                st.write(f"**Total columns**: {len(metadata_df.columns)}")
                st.dataframe(metadata_df.head())

            # Show cell barcodes
            with st.expander("🔍 Cell barcodes (first 10)", expanded=False):
                st.write(metadata_df.index[:10].tolist())

    else:  # TSV file
        metadata_file = st.file_uploader(
            "Upload TSV file (metadata source)",
            type=['tsv', 'txt', 'csv'],
            key='metadata_tsv',
            help="metaDatatheincludemuTSVFile（1colidxiscell barcode）"
        )

        if metadata_file is not None:
            with st.spinner("Loading TSV file..."):
                # Try to detect separator
                file_content = metadata_file.read().decode('utf-8')
                metadata_file.seek(0)  # Reset file pointer

                # Detect separator
                first_line = file_content.split('\n')[0]
                if '\t' in first_line:
                    sep = '\t'
                elif ',' in first_line:
                    sep = ','
                else:
                    sep = '\t'

                metadata_df = pd.read_csv(metadata_file, sep=sep, index_col=0)
                st.success(f"✓ Loaded TSV file: {len(metadata_df)} cells, {len(metadata_df.columns)} metadata columns")

                # Show metadata
                with st.expander("📊 Available metadata columns", expanded=True):
                    st.write(f"**Total columns**: {len(metadata_df.columns)}")
                    st.dataframe(metadata_df.head())

                # Show cell barcodes
                with st.expander("🔍 Cell barcodes (first 10)", expanded=False):
                    st.write(metadata_df.index[:10].tolist())

    # Step 3: Select columns and merge
    if metadata_file is not None and 'metadata_df' in locals():
        st.markdown("---")
        st.header("Step 3: Cell barcode matching options")

        # Cell barcode preprocessing options
        st.subheader("🔧 Cell barcode preprocessing")

        col1, col2 = st.columns(2)

        with col1:
            st.write("**Main h5ad preprocessing**")
            remove_main_prefix = st.checkbox(
                "Remove prefix from main h5ad cell barcodes",
                value=False,
                help="Example: 'Data1_AAACCCACAAGACCTT-1' → 'AAACCCACAAGACCTT-1'"
            )

            if remove_main_prefix:
                main_prefix_delimiter = st.text_input(
                    "Delimiter for prefix removal (main)",
                    value="_",
                    help="purefuikusutheareacutrutextchar（Example: '_'）"
                )
                main_remove_parts = st.number_input(
                    "Number of prefix parts to remove (main)",
                    min_value=1,
                    max_value=5,
                    value=1,
                    help="areacutritextcharwithdivratioshitaocctodeleteremovedobeforepartdivofnum"
                )

            remove_main_suffix = st.checkbox(
                "Remove suffix from main h5ad cell barcodes",
                value=False,
                help="Example: 'AAACCCACAAGACCTT-1' → 'AAACCCACAAGACCTT'"
            )

            if remove_main_suffix:
                main_suffix_delimiter = st.text_input(
                    "Delimiter for suffix removal (main)",
                    value="-",
                    help="safuikusutheareacutrutextchar（Example: '-'）"
                )

        with col2:
            st.write("**Metadata preprocessing**")
            remove_meta_prefix = st.checkbox(
                "Remove prefix from metadata cell barcodes",
                value=False,
                help="Example: 'E17_EC_1_AGATGAATCGAGTTGT-1' → 'AGATGAATCGAGTTGT-1'"
            )

            if remove_meta_prefix:
                meta_prefix_delimiter = st.text_input(
                    "Delimiter for prefix removal (metadata)",
                    value="_",
                    help="purefuikusutheareacutrutextchar（Example: '_'）"
                )
                meta_remove_parts = st.number_input(
                    "Number of prefix parts to remove (metadata)",
                    min_value=1,
                    max_value=5,
                    value=3,
                    help="areacutritextcharwithdivratioshitaocctodeleteremovedobeforepartdivofnum"
                )

            remove_meta_suffix = st.checkbox(
                "Remove suffix from metadata cell barcodes",
                value=False,
                help="Example: 'AGATGAATCGAGTTGT-1' → 'AGATGAATCGAGTTGT'"
            )

            if remove_meta_suffix:
                meta_suffix_delimiter = st.text_input(
                    "Delimiter for suffix removal (metadata)",
                    value="-",
                    help="safuikusutheareacutrutextchar（Example: '-'）"
                )

        # Apply preprocessing
        main_index_processed = adata_main.obs.index.copy()
        meta_index_processed = metadata_df.index.copy()

        if remove_main_prefix:
            main_index_processed = pd.Index([
                main_prefix_delimiter.join(cell.split(main_prefix_delimiter)[main_remove_parts:])
                for cell in main_index_processed
            ])

        if remove_main_suffix:
            main_index_processed = pd.Index([
                cell.split(main_suffix_delimiter)[0] if main_suffix_delimiter in cell else cell
                for cell in main_index_processed
            ])

        if remove_meta_prefix:
            meta_index_processed = pd.Index([
                meta_prefix_delimiter.join(cell.split(meta_prefix_delimiter)[meta_remove_parts:])
                for cell in meta_index_processed
            ])

        if remove_meta_suffix:
            meta_index_processed = pd.Index([
                cell.split(meta_suffix_delimiter)[0] if meta_suffix_delimiter in cell else cell
                for cell in meta_index_processed
            ])

        # Preview processed barcodes
        with st.expander("👀 Preview processed cell barcodes", expanded=True):
            preview_col1, preview_col2 = st.columns(2)

            with preview_col1:
                st.write("**Main h5ad (first 5)**")
                st.write("Original:", adata_main.obs.index[:5].tolist())
                if remove_main_prefix or remove_main_suffix:
                    st.write("Processed:", main_index_processed[:5].tolist())

            with preview_col2:
                st.write("**Metadata (first 5)**")
                st.write("Original:", metadata_df.index[:5].tolist())
                if remove_meta_prefix or remove_meta_suffix:
                    st.write("Processed:", meta_index_processed[:5].tolist())

        st.markdown("---")
        st.header("Step 4: Select metadata columns to add")

        # Check for overlapping cell barcodes (using processed indices)
        common_cells = list(set(main_index_processed) & set(meta_index_processed))
        st.info(f"📊 **Matching cells**: {len(common_cells)} / {adata_main.n_obs} cells in main h5ad")

        if len(common_cells) == 0:
            st.error("❌ No matching cell barcodes found! Please check cell barcode preprocessing options.")
            st.write("**Main h5ad cell barcodes (processed, first 5):**", main_index_processed[:5].tolist())
            st.write("**Metadata cell barcodes (processed, first 5):**", meta_index_processed[:5].tolist())
        else:
            # Select columns to add
            available_columns = metadata_df.columns.tolist()

            selected_columns = st.multiselect(
                "Select metadata columns to add",
                available_columns,
                default=available_columns,
                help="addaddshitaimetaDatacoltheSelect（multinumSelectcan）"
            )

            if selected_columns:
                st.write(f"**Selected {len(selected_columns)} columns:**")
                st.write(selected_columns)

                # Check for column name conflicts
                existing_columns = adata_main.obs.columns.tolist()
                conflicts = [col for col in selected_columns if col in existing_columns]

                if conflicts:
                    st.warning(f"⚠️ **Column name conflicts detected**: {conflicts}")
                    overwrite_mode = st.radio(
                        "How to handle conflicts?",
                        ["Overwrite existing columns", "Skip conflicting columns", "Add suffix to new columns"],
                        help="alreadyexistofcolnameandweightmultidoplacematchofprocprocwaymethod"
                    )
                else:
                    overwrite_mode = None

                # Create mapping between processed and original indices
                # Create a temporary dataframe with processed indices for metadata
                metadata_df_processed = metadata_df.copy()
                metadata_df_processed.index = meta_index_processed

                # Preview merge
                with st.expander("👀 Preview merge result", expanded=False):
                    preview_df = adata_main.obs.copy()
                    preview_df.index = main_index_processed

                    for col in selected_columns:
                        if overwrite_mode == "Skip conflicting columns" and col in conflicts:
                            continue
                        elif overwrite_mode == "Add suffix to new columns" and col in conflicts:
                            new_col_name = f"{col}_new"
                            preview_df[new_col_name] = metadata_df_processed[col]
                        else:
                            preview_df[col] = metadata_df_processed[col]

                    # Restore original index for display
                    preview_df.index = adata_main.obs.index
                    st.dataframe(preview_df.head(10))

                # Merge button
                if st.button("🔀 Merge metadata", type="primary"):
                    with st.spinner("Merging metadata..."):
                        # Create a copy of main adata
                        adata_merged = adata_main.copy()

                        # Create temporary dataframes with processed indices for proper alignment
                        temp_main = pd.DataFrame(index=main_index_processed)
                        temp_main['_original_index'] = adata_main.obs.index

                        temp_meta = metadata_df.copy()
                        temp_meta.index = meta_index_processed

                        # Add selected columns
                        n_added = 0
                        n_skipped = 0
                        added_columns = []

                        for col in selected_columns:
                            if overwrite_mode == "Skip conflicting columns" and col in conflicts:
                                n_skipped += 1
                                continue
                            elif overwrite_mode == "Add suffix to new columns" and col in conflicts:
                                new_col_name = f"{col}_new"
                                # Align by processed index, then reindex to original
                                temp_main[new_col_name] = temp_meta[col]
                                temp_main_reindexed = temp_main.set_index('_original_index')
                                adata_merged.obs[new_col_name] = temp_main_reindexed[new_col_name]
                                added_columns.append(new_col_name)
                                n_added += 1
                            else:
                                # Align by processed index, then reindex to original
                                temp_main[col] = temp_meta[col]
                                temp_main_reindexed = temp_main.set_index('_original_index')
                                adata_merged.obs[col] = temp_main_reindexed[col]
                                added_columns.append(col)
                                n_added += 1

                        st.success(f"✅ Merged metadata successfully!")
                        st.write(f"- **Added columns**: {n_added}")
                        if n_skipped > 0:
                            st.write(f"- **Skipped columns**: {n_skipped}")
                        st.write(f"- **Total metadata columns**: {len(adata_merged.obs.columns)}")

                        # Show merged data
                        with st.expander("📊 Merged metadata", expanded=False):
                            st.dataframe(adata_merged.obs.head(10))

                        # Save to session state
                        st.session_state['adata_merged'] = adata_merged
                        st.session_state['original_filename'] = main_h5ad.name

        # Step 4: Download
        if 'adata_merged' in st.session_state:
            st.markdown("---")
            st.header("Step 4: Download merged h5ad file")

            adata_to_save = st.session_state['adata_merged']
            original_filename = st.session_state['original_filename']

            # Generate output filename
            if original_filename.endswith('.h5ad'):
                output_filename = original_filename.replace('.h5ad', '.mod.h5ad')
            else:
                output_filename = f"{original_filename}.mod.h5ad"

            st.write(f"**Output filename**: `{output_filename}`")
            st.write(f"**Cells**: {adata_to_save.n_obs}")
            st.write(f"**Genes**: {adata_to_save.n_vars}")
            st.write(f"**Metadata columns**: {len(adata_to_save.obs.columns)}")

            # Save to temporary file and prepare download
            with tempfile.NamedTemporaryFile(delete=False, suffix='.h5ad') as tmp:
                adata_to_save.write_h5ad(tmp.name)
                tmp_path = tmp.name

            try:
                with open(tmp_path, 'rb') as f:
                    h5ad_data = f.read()

                st.download_button(
                    label="📥 Download merged h5ad file",
                    data=h5ad_data,
                    file_name=output_filename,
                    mime="application/octet-stream",
                    help=f"Download {output_filename}"
                )
            finally:
                os.unlink(tmp_path)

else:
    st.info("👆 Please upload a main h5ad file to begin")
