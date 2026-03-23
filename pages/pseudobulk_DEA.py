# pseudobulk_DEA.py
# Pseudobulk DEA with RUV correction
# h5ad -> pseudobulk -> RUV -> edgeR/DESeq2

import streamlit as st
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import os
import re
import io
import time
import tempfile
import shutil
import zipfile
import scipy
import matplotlib.pyplot as plt
import decoupler as dc
from helper_func import clear_old_directories, clear_old_files, mk_temp_dir

from pathlib import Path
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import StrVector


st.set_page_config(page_title="Pseudobulk DEA (RUV)", page_icon="🧬")

# ============================================================
# Utility functions
# ============================================================

def clean_counts_layer(adata):
    if 'counts' in adata.layers:
        if scipy.sparse.issparse(adata.layers['counts']):
            data = adata.layers['counts'].data
            nan_count = np.sum(np.isnan(data))
            inf_count = np.sum(np.isinf(data))
            data[np.isnan(data)] = 0
            data[np.isinf(data)] = 0
            adata.layers['counts'].data = data
            return nan_count, inf_count
        else:
            nan_count = np.sum(np.isnan(adata.layers['counts']))
            inf_count = np.sum(np.isinf(adata.layers['counts']))
            adata.layers['counts'] = np.nan_to_num(
                adata.layers['counts'], nan=0.0, posinf=0.0, neginf=0.0)
            return nan_count, inf_count
    return 0, 0


def clean_column_names(adata):
    def clean_name(name):
        cleaned = re.sub(r'[^\w\s.]', '', name)
        cleaned = re.sub(r'\s+', '_', cleaned)
        if not cleaned[0].isalpha() and cleaned[0] != '_':
            cleaned = '_' + cleaned
        return cleaned

    name_mapping = {col: clean_name(col) for col in adata.obs.columns}
    seen = {}
    for orig, cleaned in name_mapping.items():
        if cleaned in seen:
            counter = 1
            while f"{cleaned}_{counter}" in seen:
                counter += 1
            name_mapping[orig] = f"{cleaned}_{counter}"
        seen[name_mapping[orig]] = orig

    changed = {o: c for o, c in name_mapping.items() if o != c}
    if changed:
        adata.obs.rename(columns=name_mapping, inplace=True)
    return adata, changed, bool(changed)


def reshape_pseudobulk_per_celltype(pdata, sample_col, groups_col):
    """Reshape pseudobulk AnnData to per-cell-type count DataFrames."""
    data = pdata.X.toarray() if scipy.sparse.issparse(pdata.X) else pdata.X
    gene_names = pdata.var_names
    unique_groups = pdata.obs[groups_col].unique()

    result = {}
    for group in unique_groups:
        group_mask = pdata.obs[groups_col] == group
        group_data = data[group_mask]
        group_samples = pdata.obs.loc[group_mask, sample_col]

        df = pd.DataFrame(group_data, columns=gene_names, index=group_samples)
        df = df.T  # genes x samples
        df.columns = [str(c) for c in df.columns]
        result[str(group)] = df

    return result


def remove_common_suffix(strings):
    if not strings or len(strings) == 0:
        return []
    min_length = min(len(s) for s in strings)
    suffix_length = 0
    for i in range(1, min_length + 1):
        suffix = strings[0][-i:]
        if all(s.endswith(suffix) for s in strings):
            suffix_length = i
        else:
            break
    if suffix_length == 0:
        return strings
    return [s[:-suffix_length] for s in strings]


def remove_sample_num(s):
    return re.sub(r'_?\d+$', '', s)


def get_group_from_sample_names(sample_names):
    """Auto-detect group names from sample column names."""
    names = [str(n) for n in sample_names]
    groups = remove_common_suffix(names)
    groups = [remove_sample_num(g) for g in groups]
    groups = [g.replace('_', '.') for g in groups]
    return groups


@st.cache_data
def read_h5ad(file):
    adata = sc.read_h5ad(file)
    return adata


@st.cache_data
def calc_pseudobulk(_adata, groups_col, sample_col, min_cells, min_counts,
                    cache_key=None, skip_checks=False):
    return dc.get_pseudobulk(
        _adata,
        groups_col=groups_col,
        sample_col=sample_col,
        layer='counts',
        mode='sum',
        min_cells=min_cells,
        min_counts=min_counts,
        skip_checks=skip_checks
    )


def generate_pbps(adata, pdata, groups_col, sample_col, batch_col,
                   bio_vars=None, n_pbps=10, seed=42):
    """
    Generate Pseudobulk Pseudosamples (PBPS) from single-cell data.

    Following ruvPBPS algorithm:
    1. Define biological subgroup (bsg) by cell type (groups_col) x bio_vars
    2. Define group by bsg x NVar (batch_col)
    3. Compute average cell count per bsg
    4. Sample cells with replacement from each group
    5. Aggregate to pseudobulk counts
    6. PBPS from same bsg share subject_id -> replicates in RUVIII

    Parameters:
        bio_vars: list of adata.obs column names for biological variables
                  (e.g., ['treatment', 'timepoint']). Cell type is always included.
        batch_col: adata.obs column name for NVar (technical/batch variable)

    Returns dict of {celltype: (pbps_counts_df, pbps_subject_ids)}
    """
    cell_meta = adata.obs.copy()

    if bio_vars is None:
        bio_vars = []

    # Get the counts matrix (same layer used for pseudobulk)
    if 'counts' in adata.layers:
        counts_mat = adata.layers['counts']
    else:
        counts_mat = adata.X

    if scipy.sparse.issparse(counts_mat):
        counts_mat = counts_mat.tocsc()

    gene_names = adata.var_names
    cell_types = cell_meta[groups_col].unique()

    rng = np.random.RandomState(seed)

    result = {}

    for ct in cell_types:
        ct_mask = cell_meta[groups_col] == ct
        ct_meta = cell_meta[ct_mask].copy()
        ct_indices = np.where(ct_mask.values)[0]

        # Define biological subgroups (bsg) within this cell type
        # bsg = bio_vars combinations (e.g., treatment x timepoint)
        if bio_vars:
            ct_meta['_bsg'] = ct_meta[bio_vars].astype(str).agg('_'.join, axis=1)
        else:
            ct_meta['_bsg'] = 'all'  # single bsg if no bio_vars

        # Define groups: bsg x batch
        ct_meta['_group'] = ct_meta['_bsg'] + '.' + ct_meta[batch_col].astype(str)

        bsg_levels = ct_meta['_bsg'].unique()

        pbps_counts_list = []
        pbps_names = []
        pbps_subject_ids = []

        for bsg in bsg_levels:
            bsg_mask = ct_meta['_bsg'] == bsg
            bsg_meta = ct_meta[bsg_mask]

            # NVar levels within this bsg
            batch_levels = bsg_meta[batch_col].unique()
            if len(batch_levels) < 2:
                # Need at least 2 NVar levels to create PBPS for this bsg
                continue

            # Average cell count per sample within this bsg
            sample_counts = bsg_meta.groupby(sample_col).size()
            avg_cells = max(1, int(np.round(sample_counts.mean())))

            for rep_i in range(n_pbps):
                for batch_val in batch_levels:
                    batch_mask = bsg_meta[batch_col] == batch_val
                    batch_cell_indices = ct_indices[bsg_mask.values][batch_mask.values]

                    if len(batch_cell_indices) == 0:
                        continue

                    # Sample with replacement
                    sampled_indices = rng.choice(
                        batch_cell_indices, size=avg_cells, replace=True)

                    # Aggregate counts
                    if scipy.sparse.issparse(counts_mat):
                        sampled_counts = np.asarray(
                            counts_mat[sampled_indices, :].sum(axis=0)).flatten()
                    else:
                        sampled_counts = counts_mat[sampled_indices, :].sum(axis=0).flatten()

                    pbps_name = f"pbps_{rep_i+1}.{ct}_{bsg}_{batch_val}"
                    pbps_counts_list.append(sampled_counts)
                    pbps_names.append(pbps_name)
                    # PBPS from same bsg share subject_id -> replicates
                    pbps_subject_ids.append(f"pbps_{ct}_{bsg}")

        if pbps_counts_list:
            pbps_df = pd.DataFrame(
                np.column_stack(pbps_counts_list),
                index=gene_names,
                columns=pbps_names
            )
            result[str(ct)] = (pbps_df, pbps_subject_ids)

    return result


def pass_df_to_r(df, r_varname):
    """Pass a pandas DataFrame to R as a matrix with row/col names."""
    # Avoid py2rpy conversion issues by passing data as plain Python floats
    # and constructing the matrix entirely in R via string evaluation.
    arr = df.values.astype(float)
    nr, nc = arr.shape
    # Plain Python float list (not numpy) -> rpy2 FloatVector (Sexp, no conversion needed)
    flat_values = [float(x) for x in arr.flatten()]
    r_vec = ro.FloatVector(flat_values)
    # Assign the FloatVector (Sexp type passes through py2rpy as-is)
    ro.globalenv['._tmp_pbdea_'] = r_vec
    # Reshape to matrix in R (avoids returning an R object to Python)
    ro.r(f'{r_varname} <- matrix(._tmp_pbdea_, nrow={nr}L, ncol={nc}L, byrow=TRUE)')
    ro.r('rm(._tmp_pbdea_)')
    # Set row and column names
    row_str = ','.join(repr(str(x)) for x in df.index)
    col_str = ','.join(repr(str(x)) for x in df.columns)
    ro.r(f'rownames({r_varname}) <- c({row_str})')
    ro.r(f'colnames({r_varname}) <- c({col_str})')


def run_dea_for_celltype(count_df, groups, celltype_name, dea_method,
                         ruv_method, ruv_k, ruv_alpha, res_dir,
                         use_glmtreat=False, treat_lfc=0.585,
                         deseq2_shrinkage="ashr", ref_group=None,
                         pbps_counts=None, pbps_subject_ids=None,
                         sva_n=None):
    """Run DEA for a single cell type."""

    ct_dir = os.path.join(res_dir, celltype_name.replace(' ', '_').replace('/', '_'))
    os.makedirs(ct_dir, exist_ok=True)

    # Pass count matrix to R (integer)
    count_int = count_df.round(0).astype(int)
    pass_df_to_r(count_int, 'pb_counts')

    # Pass group info
    ro.r.assign('pb_group', StrVector(groups))

    # Set res_dir in R
    ro.r.assign('pb_res_dir', ct_dir)

    # Source R functions (ensure conda library path)
    r_func_path = os.path.join(os.path.dirname(__file__), 'pseudobulk_dea_func.R')
    ro.r(f'.libPaths(c("{Path.home()}/anaconda3/envs/shiny/lib/R/library", .libPaths()))')
    ro.r(f'source("{r_func_path}")')

    ruv_W = None
    ruv_info = {}

    # ---- Batch correction (SVA / RUV) ----
    if ruv_method == "SVA":
        sva_n_arg = "NULL" if sva_n is None else str(int(sva_n))
        st.write(f"  Running SVA (n_sv={sva_n_arg})...")
        ro.r(f'''
        sva_result <- run_SVA(pb_counts, pb_group, n_sv={sva_n_arg},
                              res_dir=pb_res_dir)
        ruv_W_mat <- sva_result$SV
        pb_counts_filtered <- sva_result$filtered_counts
        ''')
        ruv_W = "ruv_W_mat"
        ro.r('pb_counts <- pb_counts_filtered')

        n_sv_used = int(ro.r('sva_result$n_sv')[0])
        n_sv_be = ro.r('sva_result$n_sv_be')[0]
        n_sv_leek = ro.r('sva_result$n_sv_leek')[0]
        ruv_info['n_sv'] = n_sv_used
        ruv_info['n_sv_be'] = int(n_sv_be) if not np.isnan(n_sv_be) else None
        ruv_info['n_sv_leek'] = int(n_sv_leek) if not np.isnan(n_sv_leek) else None
        ruv_info['method'] = 'SVA (svaseq)'

    elif ruv_method == "RUVg (RUVSeq)":
        st.write(f"  Running RUVg (k={ruv_k})...")
        ro.r.assign('ruv_k', ruv_k)
        ro.r.assign('ruv_alpha', ruv_alpha)
        ro.r('''
        ruv_result <- run_RUVg(pb_counts, pb_group, k=ruv_k,
                               alpha=ruv_alpha, res_dir=pb_res_dir)
        ruv_W_mat <- ruv_result$W
        pb_counts_filtered <- ruv_result$filtered_counts
        ''')
        ruv_W = "ruv_W_mat"
        # Update counts to filtered
        ro.r('pb_counts <- pb_counts_filtered')

        with localconverter(ro.default_converter + pandas2ri.converter):
            w_mat = ro.conversion.rpy2py(ro.r('ruv_result$W'))
        ruv_info['W'] = w_mat
        ruv_info['n_control_genes'] = int(ro.r('length(ruv_result$control_genes)')[0])

    elif ruv_method == "RUV2 (ruv)":
        st.write(f"  Running RUV2 (k={ruv_k})...")
        ro.r.assign('ruv_k', ruv_k)
        ro.r.assign('ruv_alpha', ruv_alpha)
        try:
            ro.r('''
            ruv_result <- run_RUV2(pb_counts, pb_group, k=ruv_k,
                                   alpha=ruv_alpha, res_dir=pb_res_dir)
            ruv_W_mat <- ruv_result$W
            pb_counts_filtered <- ruv_result$filtered_counts
            ''')
            ruv_W = "ruv_W_mat"
            ro.r('pb_counts <- pb_counts_filtered')

            with localconverter(ro.default_converter + pandas2ri.converter):
                w_mat = ro.conversion.rpy2py(ro.r('ruv_result$W'))
            ruv_info['W'] = w_mat
            ruv_info['n_control_genes'] = int(ro.r('length(ruv_result$control_genes)')[0])
        except Exception as e:
            st.warning(
                f"  RUV2 failed: {str(e).split(chr(10))[0]}\n\n"
                "  When the number of samples is small and there are many groups, "
                "the internal matrix computation of RUV2 may become singular. "
                "**RUVg (RUVSeq)** is recommended.\n\n"
                "  Proceeding with DEA without RUV correction.")
            ruv_W = None

    elif ruv_method == "RUVIII (ruv)":
        st.write(f"  Running RUVIII (k={ruv_k})...")
        ro.r.assign('ruv_k', ruv_k)
        try:
            ro.r('''
            ruv_result <- run_RUVIII(pb_counts, pb_group, k=ruv_k,
                                     res_dir=pb_res_dir)
            ruv_W_mat <- ruv_result$W
            pb_counts_filtered <- ruv_result$filtered_counts
            ''')
            ruv_W = "ruv_W_mat"
            ro.r('pb_counts <- pb_counts_filtered')

            with localconverter(ro.default_converter + pandas2ri.converter):
                w_mat = ro.conversion.rpy2py(ro.r('ruv_result$W'))
            ruv_info['W'] = w_mat
            ruv_info['n_control_genes'] = int(ro.r('length(ruv_result$control_genes)')[0])
        except Exception as e:
            st.warning(
                f"  RUVIII failed: {str(e).split(chr(10))[0]}\n\n"
                "  This may be due to insufficient replicates. "
                "**RUVg (RUVSeq)** is recommended.\n\n"
                "  Proceeding with DEA without RUV correction.")
            ruv_W = None

    elif ruv_method == "RUVIII PBPS":
        if pbps_counts is None or pbps_subject_ids is None:
            st.error("RUVIII PBPS: PBPS data has not been generated. Please check the batch variable.")
        else:
            st.write(f"  Running RUVIII PBPS (k={ruv_k}, {pbps_counts.shape[1]} PBPS)...")

            # Combine original + PBPS counts
            # Ensure same genes (intersection)
            common_genes = count_df.index.intersection(pbps_counts.index)
            combined_counts = pd.concat([
                count_df.loc[common_genes],
                pbps_counts.loc[common_genes]
            ], axis=1)

            # Create subject_id vector: original samples get unique IDs,
            # PBPS get shared IDs per biological subgroup
            orig_sample_names = list(count_df.columns)
            orig_subject_ids = orig_sample_names  # each original sample is its own "subject"
            all_subject_ids = orig_subject_ids + pbps_subject_ids

            # pbps_indicator: 0 = original, 1 = PBPS
            pbps_indicator = [0] * len(orig_sample_names) + [1] * len(pbps_subject_ids)

            # Pass combined data to R
            combined_int = combined_counts.round(0).astype(int)
            pass_df_to_r(combined_int, 'pbps_combined_counts')
            ro.r.assign('pbps_subject_id', StrVector(all_subject_ids))
            ro.r.assign('pbps_indicator', ro.IntVector(pbps_indicator))
            ro.r.assign('ruv_k', ruv_k)
            # Pass orig_group explicitly (pb_group is already set at line 296,
            # but re-assign to pbps_orig_group to be safe)
            ro.r.assign('pbps_orig_group', StrVector(groups))

            ro.r('''
            ruv_result <- run_RUVIII_PBPS(pbps_combined_counts, pbps_subject_id,
                                           pbps_indicator, k=ruv_k,
                                           res_dir=pb_res_dir,
                                           orig_group=pbps_orig_group)
            ruv_W_mat <- ruv_result$W
            pb_counts_filtered <- ruv_result$filtered_counts
            ''')
            ruv_W = "ruv_W_mat"
            # Replace pb_counts with filtered original-only counts
            ro.r('pb_counts <- pb_counts_filtered')

            with localconverter(ro.default_converter + pandas2ri.converter):
                w_mat = ro.conversion.rpy2py(ro.r('ruv_result$W'))
            ruv_info['W'] = w_mat
            ruv_info['n_genes'] = int(ro.r('ruv_result$n_genes')[0])
            ruv_info['n_orig'] = int(ro.r('ruv_result$n_orig')[0])
            ruv_info['n_pbps'] = int(ro.r('ruv_result$n_pbps')[0])
            ruv_info['method'] = 'RUVIII PBPS (all genes as negative controls)'

    # ---- QC plots (before/after RUV) ----
    ro.r('''
    plot_RLE(pb_counts, pb_group, title="RLE Plot (before DEA)",
             res_dir=pb_res_dir, filename="RLE_before_DEA.png")
    plot_PCA_comparison(pb_counts, pb_group, title="PCA (before DEA)",
                        res_dir=pb_res_dir, filename="PCA_before_DEA.png")
    ''')

    # ---- DEA ----
    results_dict = {}

    if dea_method == "edgeR":
        ruv_w_arg = ruv_W if ruv_W else "NULL"
        ro.r.assign('use_glmtreat', use_glmtreat)
        ro.r.assign('treat_lfc', treat_lfc)
        ro.r(f'''
        edger_results <- run_edgeR(pb_counts, pb_group,
                                   ruv_W = {ruv_w_arg},
                                   use_glmtreat = use_glmtreat,
                                   treat_lfc = treat_lfc,
                                   res_dir = pb_res_dir)
        ''')

        # Extract results
        comp_names = list(ro.r('names(edger_results)'))
        for comp in comp_names:
            ro.r.assign('comp_name', comp)
            r_tab = ro.r('edger_results[[comp_name]]')
            with localconverter(ro.default_converter + pandas2ri.converter):
                df_result = ro.conversion.rpy2py(r_tab)
            if 'gene' in df_result.columns:
                df_result = df_result.set_index('gene')
            results_dict[comp] = df_result

    elif dea_method == "DESeq2":
        ruv_w_arg = ruv_W if ruv_W else "NULL"
        ro.r.assign('shrinkage_type', deseq2_shrinkage)
        if ref_group:
            ro.r.assign('ref_group', ref_group)
            ref_r = "ref_group"
        else:
            ref_r = "NULL"
        ro.r(f'''
        deseq2_results <- run_DESeq2(pb_counts, pb_group,
                                     ruv_W = {ruv_w_arg},
                                     ref_group = {ref_r},
                                     shrinkage_type = shrinkage_type,
                                     res_dir = pb_res_dir)
        ''')

        comp_names = list(ro.r('names(deseq2_results)'))
        for comp in comp_names:
            ro.r.assign('comp_name', comp)
            r_tab = ro.r('deseq2_results[[comp_name]]')
            with localconverter(ro.default_converter + pandas2ri.converter):
                df_result = ro.conversion.rpy2py(r_tab)
            if 'gene' in df_result.columns:
                df_result = df_result.set_index('gene')
            results_dict[comp] = df_result

    # Save individual comparison TSVs
    for comp, df_res in results_dict.items():
        tsv_path = os.path.join(ct_dir, f"{comp}.tsv")
        df_res.to_csv(tsv_path, sep='\t')

    # Save combined TSV (all comparisons in one file, DESeq2.py style)
    if results_dict:
        combined = pd.DataFrame()
        for comp, df_res in results_dict.items():
            if dea_method == "edgeR":
                col_map = {'logFC': 'log2FC', 'PValue': 'pvalue', 'FDR': 'adj.pvalue'}
            else:  # DESeq2
                col_map = {'log2FoldChange': 'log2FC', 'pvalue': 'pvalue', 'padj': 'adj.pvalue'}
            for orig_col, short_name in col_map.items():
                if orig_col in df_res.columns:
                    combined[f'{comp}.{short_name}'] = df_res[orig_col]
        combined_path = os.path.join(ct_dir, f"{celltype_name}_all_comparisons.tsv")
        combined.to_csv(combined_path, sep='\t')

    return results_dict, ruv_info, ct_dir


# ============================================================
# Session state initialization
# ============================================================
if "pb_dea_temp_dir" not in st.session_state or \
        type(st.session_state.pb_dea_temp_dir) == bool:
    temp_dir = os.path.join("temp", str(round(time.time())) + "_pbdea")
    os.makedirs("temp", exist_ok=True)
    clear_old_directories("temp")
    clear_old_files("temp")
    os.makedirs(temp_dir, exist_ok=True)
    st.session_state.pb_dea_temp_dir = temp_dir
else:
    temp_dir = st.session_state.pb_dea_temp_dir
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir, exist_ok=True)

res_dir = os.path.join(temp_dir, "res")
os.makedirs(res_dir, exist_ok=True)

if "pb_dea_pseudobulk" not in st.session_state:
    st.session_state.pb_dea_pseudobulk = None  # dict of {celltype: DataFrame}
if "pb_dea_pdata" not in st.session_state:
    st.session_state.pb_dea_pdata = None
if "pb_dea_sample_col" not in st.session_state:
    st.session_state.pb_dea_sample_col = None
if "pb_dea_groups_col" not in st.session_state:
    st.session_state.pb_dea_groups_col = None
if "pb_dea_adata" not in st.session_state:
    st.session_state.pb_dea_adata = None  # keep adata for PBPS generation


# ============================================================
# STEP 1: Data Input & Pseudobulk
# ============================================================
st.markdown("## Pseudobulk DEA with Batch Correction")
st.markdown("h5ad -> Pseudobulk -> Batch correction (RUV / SVA) -> DEA (edgeR / DESeq2)")
st.markdown("""
**Reference:** Prieto Leon S, De Troyer E, Geys H, Van den Berge K, Thas O.
Removal of unwanted variation in pseudobulk analysis of single-cell RNA sequencing data
and the leveraging of pseudoreplicates.
*NAR Genomics and Bioinformatics*. 2025;7:lqaf179.
[doi:10.1093/nargab/lqaf179](https://doi.org/10.1093/nargab/lqaf179)

**Recommendations from the paper:**
1. **Trail 2** (estimate W per cell type) is recommended. This app adopts this approach.
2. **RUVIII PBPS** (recommended): Misspecification of negative control (NC) genes can lead to FDR inflation and TPR reduction (paper Appendix F).
   RUVIII PBPS uses all genes as NC genes to avoid this problem. Technical replicates are also not required.
3. **RUV2 / RUVIII**: Effective when true NC genes are available. Controls FDR even under confounding.

Note on NC genes: The paper validated using stably expressed gene sets [Deeke 2020, Lin 2019] but
reported inconsistent results. The RUV2/RUVIII/RUVg implementations in this app use **empirical control genes**
(genes with p > threshold in an initial DE test), which differs from the paper's validation approach.
""")
st.markdown("---")

st.markdown("### Step 1: Data Input & Pseudobulk Aggregation")

uploaded_file = st.file_uploader("Upload h5ad file", type=['h5ad'])

if uploaded_file is not None:
    adata = read_h5ad(uploaded_file)

    st.write(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")

    # ---- Layer information display ----
    with st.expander("Data structure overview", expanded=True):
        # adata.X info
        if hasattr(adata.X, 'toarray'):
            x_sample = adata.X[:100, :100].toarray()
        else:
            x_sample = adata.X[:100, :100]
        x_max = float(np.max(x_sample))
        x_has_decimals = not np.allclose(x_sample, np.round(x_sample))
        x_sparsity = float((x_sample == 0).sum() / x_sample.size * 100)
        st.write(f"**adata.X** -- Max: {x_max:.2f}, Sparsity: {x_sparsity:.1f}%, "
                 f"Decimals: {'Yes' if x_has_decimals else 'No'}")

        # adata.raw.X info
        if hasattr(adata, 'raw') and adata.raw is not None:
            if hasattr(adata.raw.X, 'toarray'):
                raw_sample = adata.raw.X[:100, :100].toarray()
            else:
                raw_sample = adata.raw.X[:100, :100]
            raw_max = float(np.max(raw_sample))
            raw_has_dec = not np.allclose(raw_sample, np.round(raw_sample))
            raw_sparsity = float((raw_sample == 0).sum() / raw_sample.size * 100)
            st.write(f"**adata.raw.X** -- Max: {raw_max:.2f}, Sparsity: {raw_sparsity:.1f}%, "
                     f"Decimals: {'Yes' if raw_has_dec else 'No'}")

        # layers info
        if adata.layers:
            for lname in adata.layers.keys():
                ldata = adata.layers[lname]
                if hasattr(ldata, 'toarray'):
                    l_sample = ldata[:100, :100].toarray()
                else:
                    l_sample = ldata[:100, :100]
                l_max = float(np.max(l_sample))
                l_dec = not np.allclose(l_sample, np.round(l_sample))
                l_sparsity = float((l_sample == 0).sum() / l_sample.size * 100)
                st.write(f"**adata.layers['{lname}']** -- Max: {l_max:.2f}, "
                         f"Sparsity: {l_sparsity:.1f}%, Decimals: {'Yes' if l_dec else 'No'}")
        else:
            st.write("layers: None")

    # ---- Build layer selection options ----
    st.markdown("**Please select raw counts (integer values) for DEA.**")
    layer_options = ["adata.X"]
    if hasattr(adata, 'raw') and adata.raw is not None:
        layer_options.append("adata.raw.X")
    if adata.layers:
        for lname in adata.layers.keys():
            layer_options.append(f"adata.layers['{lname}']")

    # Default: prefer raw.X or 'counts' layer
    default_idx = 0
    if "adata.raw.X" in layer_options:
        default_idx = layer_options.index("adata.raw.X")
    for i, opt in enumerate(layer_options):
        if "'counts'" in opt:
            default_idx = i
            break

    selected_layer = st.selectbox(
        "Counts data source:",
        options=layer_options,
        index=default_idx,
        help="Select the layer containing raw counts for pseudobulk + DEA"
    )

    # ---- Apply selected layer ----
    if selected_layer == "adata.X":
        adata.layers['counts'] = adata.X
    elif selected_layer == "adata.raw.X":
        adata.layers['counts'] = adata.raw.X
        adata.X = adata.raw.X
    elif "adata.layers['" in selected_layer:
        lname = selected_layer.split("'")[1]
        adata.layers['counts'] = adata.layers[lname]
        adata.X = adata.layers[lname]

    nan_count, inf_count = clean_counts_layer(adata)
    if nan_count > 0 or inf_count > 0:
        st.info(f"Replaced NaN: {nan_count}, Inf: {inf_count} with 0.")

    # Preview selected data
    st.markdown(f"**Selected data preview ({selected_layer}):**")
    preview_X = adata.layers['counts']
    temp_preview = pd.DataFrame(
        preview_X[:5, :8].toarray() if scipy.sparse.issparse(preview_X) else preview_X[:5, :8],
        index=adata.obs_names[:5],
        columns=adata.var_names[:8]
    )
    st.dataframe(temp_preview)

    # Clean column names
    adata, name_mapping, changes_made = clean_column_names(adata)
    if changes_made:
        st.write("Column names have been cleaned:")
        for orig, cleaned in name_mapping.items():
            st.write(f"  {orig} -> {cleaned}")

    st.markdown("**adata.obs:**")
    st.write(adata.obs.head())

    # Column selection
    meta_cols = adata.obs.columns.tolist()
    meta_cols = [c for c in meta_cols if c not in
                 ['nCount_RNA', 'nFeature_RNA', 'percent.mt', 'Cell_id']]

    with st.form("pseudobulk_settings"):
        col1, col2 = st.columns(2)
        with col1:
            sample_col = st.selectbox("Sample column:", meta_cols)
        with col2:
            groups_col = st.selectbox("Cell type column:", meta_cols)

        st.markdown("##### Quality threshold")
        c1, c2 = st.columns(2)
        with c1:
            min_cells = st.number_input('Minimum number of cells:', min_value=0, value=10)
        with c2:
            min_counts = st.number_input('Minimum total counts:', min_value=0, value=100)

        st.markdown("Using mode=sum for DEA.")
        run_pseudobulk = st.form_submit_button("Run pseudobulk aggregation")

    if run_pseudobulk:
        with st.spinner("Running pseudobulk aggregation..."):
            skip_checks = False
            try:
                pdata = calc_pseudobulk(adata, groups_col, sample_col,
                                        min_cells, min_counts,
                                        cache_key=time.time())
            except ValueError:
                st.warning("Counts contain decimal values. Proceeding with computation as-is.")
                pdata = calc_pseudobulk(adata, groups_col, sample_col,
                                        min_cells, min_counts,
                                        cache_key=time.time(), skip_checks=True)

            # Reshape to per-cell-type DataFrames
            pb_data = reshape_pseudobulk_per_celltype(pdata, sample_col, groups_col)

            # Cache
            st.session_state.pb_dea_pseudobulk = pb_data
            st.session_state.pb_dea_pdata = pdata
            st.session_state.pb_dea_sample_col = sample_col
            st.session_state.pb_dea_groups_col = groups_col
            st.session_state.pb_dea_adata = adata  # keep for PBPS generation

            st.success(f"Pseudobulk aggregation complete: {len(pb_data)} cell types")

    # Show cached pseudobulk data
    if st.session_state.pb_dea_pseudobulk is not None:
        pb_data = st.session_state.pb_dea_pseudobulk

        st.markdown("#### Pseudobulk Results (cached)")

        # QC plot
        if st.session_state.pb_dea_pdata is not None:
            pdata = st.session_state.pb_dea_pdata
            s_col = st.session_state.pb_dea_sample_col
            g_col = st.session_state.pb_dea_groups_col
            try:
                dc.plot_psbulk_samples(pdata, groupby=[s_col, g_col], figsize=(12, 4))
                fig = plt.gcf()
                st.pyplot(fig)
                plt.close(fig)
            except Exception as e:
                st.warning(f"QC plot error: {e}")

        for ct_name, ct_df in pb_data.items():
            with st.expander(f"{ct_name} ({ct_df.shape[0]} genes x {ct_df.shape[1]} samples)"):
                st.dataframe(ct_df.head())

        # ============================================================
        # STEP 2: DEA Settings
        # ============================================================
        st.markdown("---")
        st.markdown("### Step 2: DEA Settings")

        # Cell type selection
        all_celltypes = list(pb_data.keys())
        run_all = st.checkbox("Run all cell types at once", value=True)
        if not run_all:
            selected_celltypes = st.multiselect(
                "Select cell types:", all_celltypes, default=all_celltypes[:1])
        else:
            selected_celltypes = all_celltypes

        # Group editing
        st.markdown("##### Group Settings")
        st.markdown("Groups are automatically inferred from sample names across all selected cell types. Edit as needed.")

        # Collect all unique sample names across all selected cell types
        all_samples_set = []
        for ct in selected_celltypes:
            for s in pb_data[ct].columns:
                if s not in all_samples_set:
                    all_samples_set.append(s)
        auto_groups = get_group_from_sample_names(all_samples_set)

        group_df = pd.DataFrame({
            'Sample': all_samples_set,
            'Group': auto_groups
        })

        edited_groups = st.data_editor(group_df, key="group_editor",
                                       disabled=["Sample"])

        # Show per-cell-type sample counts
        with st.expander("Sample counts per cell type"):
            for ct in selected_celltypes:
                ct_samples = list(pb_data[ct].columns)
                st.write(f"**{ct}**: {len(ct_samples)} samples")

        # RUV settings
        st.markdown("##### Batch Correction")
        with st.sidebar:
            st.markdown("### Batch Correction Settings")
            ruv_method = st.selectbox(
                "Correction method:",
                ["RUVIII PBPS", "RUV2 (ruv)", "RUVIII (ruv)", "RUVg (RUVSeq)", "SVA", "None"],
                index=0,
                help="""
                **RUVIII PBPS** (most recommended): Uses all genes as NC genes, avoiding NC gene misspecification risk. Technical replicates are not required. Generates pseudosamples (PBPS). Requires a batch variable (NVar).
                **RUV2** (recommended): ruv package. Trail 2 + RUV2 controls FDR even under confounding. Uses empirical control genes (differs from paper's validation).
                **RUVIII** (recommended): ruv package. Trail 2 + RUVIII also controls FDR under confounding. Requires technical replicates. Uses empirical control genes (differs from paper).
                **RUVg**: RUVSeq package. Uses empirical control genes. Not directly recommended in the paper but available as an existing implementation.
                **SVA**: Surrogate Variable Analysis. Data-driven estimation of batch effects. Independent from RUV methods.
                **None**: No correction.
                """
            )

            ruv_k = 1
            ruv_alpha = 0.1
            sva_n = 2
            n_pbps = 10
            pbps_batch_col = None
            pbps_bio_vars = []
            pbps_seed = 42
            if ruv_method == "SVA":
                sva_n = st.number_input("Number of SVs:", min_value=1, max_value=10,
                                        value=2,
                                        help="Number of surrogate variables. BE/Leek estimates will also be displayed for reference.")
            elif ruv_method != "None":
                ruv_k = st.number_input("k (number of RUV factors):", min_value=1, max_value=5,
                                        value=1,
                                        help="Recommended: start with k=1. Check QC plots and adjust.")
                if ruv_method in ["RUVg (RUVSeq)", "RUV2 (ruv)"]:
                    ruv_alpha = st.number_input(
                        "Control gene p-value threshold:",
                        min_value=0.01, max_value=1.0, value=0.1, step=0.05,
                        help="Genes with p-values above this threshold are used as control genes")
                if ruv_method == "RUVIII PBPS":
                    st.markdown("---")
                    st.markdown("#### PBPS Settings")
                    if st.session_state.pb_dea_adata is not None:
                        obs_cols = st.session_state.pb_dea_adata.obs.columns.tolist()
                        # BioVar: biological variables (treatment, timepoint, etc.)
                        # Cell type (groups_col) is automatically included
                        g_col = st.session_state.pb_dea_groups_col
                        bio_candidates = [c for c in obs_cols if c != g_col]
                        pbps_bio_vars = st.multiselect(
                            "Experimental conditions (BioVar):",
                            bio_candidates,
                            help="Select columns representing experimental conditions (e.g., condition, treatment, timepoint)."
                                 f" Cell type ({g_col}) is automatically included."
                                 " Samples with the same experimental condition will be treated as replicates."
                        )
                        # NVar: batch/technical variable (default to sample_col)
                        s_col = st.session_state.pb_dea_sample_col
                        nvar_default = obs_cols.index(s_col) if s_col in obs_cols else 0
                        pbps_batch_col = st.selectbox(
                            "Sample ID / Batch variable (NVar):",
                            obs_cols,
                            index=nvar_default,
                            help="Select a sample ID or batch variable. "
                                 "PBPS will be generated from different samples within the same experimental condition."
                        )
                        if pbps_bio_vars:
                            st.info(
                                f"biological subgroup (bsg) = {g_col} x "
                                f"{' x '.join(pbps_bio_vars)}\n\n"
                                f"PBPS can be generated if there are 2 or more different {pbps_batch_col} within the same bsg")
                        else:
                            st.warning(
                                "No experimental conditions selected. Grouping will be by cell type only, "
                                "which may mix samples from different experimental conditions (e.g., CTRL/HFD).")
                    else:
                        st.warning("Please run pseudobulk aggregation in Step 1 first.")
                    n_pbps = st.number_input(
                        "Number of PBPS (per group):", min_value=1, max_value=50,
                        value=10,
                        help="The paper used n=10. More PBPS leads to more stable W estimation."
                    )
                    pbps_seed = st.number_input(
                        "Random seed:", min_value=1, value=42,
                        help="Seed value for reproducibility"
                    )

            st.markdown("---")

            # DEA method
            st.markdown("### DEA Settings")
            dea_method = st.selectbox(
                "DEA method:",
                ["edgeR", "DESeq2"],
                index=0
            )

            if dea_method == "edgeR":
                test_method = st.radio(
                    "Test method:",
                    ["glmQLFTest (standard)", "glmTreat (minimum FC)"],
                    index=0
                )
                use_glmtreat = (test_method == "glmTreat (minimum FC)")
                treat_lfc = 0.585
                if use_glmtreat:
                    treat_lfc = st.number_input(
                        "Log2FC threshold:", value=0.585, min_value=0.0, step=0.1)

            elif dea_method == "DESeq2":
                deseq2_shrinkage = st.selectbox(
                    "LFC shrinkage:",
                    ["ashr", "apeglm", "normal"],
                    index=0,
                    help="ashr recommended (compatible with contrast-based usage)"
                )
                ref_group_option = st.checkbox("Specify reference group", value=False)
                ref_group = None
                if ref_group_option:
                    unique_groups = list(edited_groups['Group'].unique())
                    ref_group = st.selectbox("Reference group:", unique_groups)

        # ============================================================
        # STEP 3: Run DEA
        # ============================================================
        st.markdown("---")
        st.markdown("### Step 3: Run DEA")

        # Summary
        if ruv_method == "SVA":
            correction_str = f"- Correction: SVA" + (f" (n_sv={sva_n})" if sva_n else " (auto-estimated)")
        elif ruv_method != "None":
            correction_str = f"- Correction: {ruv_method} (k={ruv_k})"
        else:
            correction_str = "- Correction: None"
        summary_lines = [
            f"- Cell types: {len(selected_celltypes)}",
            correction_str,
            f"- DEA method: {dea_method}",
            f"- Groups: {', '.join(edited_groups['Group'].unique())}",
        ]
        if ruv_method == "RUVIII PBPS":
            summary_lines.append(f"- Batch variable: {pbps_batch_col}")
            summary_lines.append(f"- PBPS per group: {n_pbps}")
        st.markdown("**Settings:**\n" + "\n".join(summary_lines))

        if st.button("Run DEA", type="primary"):
            # Build sample -> group mapping (name-based)
            group_map = {}
            for _, row in edited_groups.iterrows():
                g = row['Group']
                g = g.replace('_', '.')
                g = re.sub(r'[^a-zA-Z0-9\.]', '', g)
                # R requires syntactically valid names (can't start with digit)
                if g and g[0].isdigit():
                    g = 'X' + g
                group_map[row['Sample']] = g

            # Generate PBPS if RUVIII PBPS selected
            pbps_data = {}
            if ruv_method == "RUVIII PBPS":
                if pbps_batch_col is None or st.session_state.pb_dea_adata is None:
                    st.error("RUVIII PBPS: adata or batch variable is not configured.")
                    st.stop()
                with st.spinner("Generating PBPS (Pseudobulk Pseudosamples)..."):
                    s_col = st.session_state.pb_dea_sample_col
                    g_col = st.session_state.pb_dea_groups_col
                    pbps_data = generate_pbps(
                        st.session_state.pb_dea_adata,
                        st.session_state.pb_dea_pdata,
                        groups_col=g_col,
                        sample_col=s_col,
                        batch_col=pbps_batch_col,
                        bio_vars=pbps_bio_vars if pbps_bio_vars else None,
                        n_pbps=n_pbps,
                        seed=pbps_seed
                    )
                    if pbps_data:
                        total_pbps = sum(df.shape[1] for df, _ in pbps_data.values())
                        st.info(f"PBPS generation complete: {len(pbps_data)} cell types, "
                                f"total {total_pbps} pseudosamples")
                    else:
                        st.error(
                            "PBPS generation failed: At least 2 batch levels are required within each biological subgroup.\n\n"
                            "Possible causes:\n"
                            "- The batch variable and experimental condition are confounded (batch ~ condition)\n"
                            "- Incorrect selection of experimental conditions (BioVar)\n\n"
                            "**Recommendation**: Use RUVg (RUVSeq) instead.")
                        st.stop()

            all_results = {}
            all_ruv_info = {}
            all_ct_dirs = {}

            progress = st.progress(0)
            total = len(selected_celltypes)

            for idx, ct_name in enumerate(selected_celltypes):
                st.markdown(f"#### {ct_name}")

                ct_df = pb_data[ct_name]
                ct_samples = list(ct_df.columns)

                # Map groups by sample name
                ct_groups = []
                unmapped = []
                for s in ct_samples:
                    if s in group_map:
                        ct_groups.append(group_map[s])
                    else:
                        unmapped.append(s)

                if unmapped:
                    st.warning(
                        f"  {ct_name}: The following samples have no group assignment: "
                        f"{', '.join(unmapped)}. Skipping.")
                    continue

                # Check minimum samples per group (warning only)
                group_counts = pd.Series(ct_groups).value_counts()
                if (group_counts < 2).any():
                    low_groups = group_counts[group_counts < 2].to_dict()
                    st.warning(
                        f"  {ct_name}: Some groups have fewer than 2 samples: {low_groups}. "
                        f"Variance estimation accuracy may be reduced.")
                if len(set(ct_groups)) < 2:
                    st.warning(
                        f"  {ct_name}: Only one group found. Cannot perform comparison. Skipping.")
                    continue

                # Get PBPS data for this cell type
                ct_pbps_counts = None
                ct_pbps_subject_ids = None
                if ruv_method == "RUVIII PBPS" and ct_name in pbps_data:
                    ct_pbps_counts, ct_pbps_subject_ids = pbps_data[ct_name]

                try:
                    with st.spinner(f"Analyzing {ct_name}..."):
                        if dea_method == "edgeR":
                            results, ruv_info, ct_dir = run_dea_for_celltype(
                                ct_df, ct_groups, ct_name, "edgeR",
                                ruv_method, ruv_k, ruv_alpha, res_dir,
                                use_glmtreat=use_glmtreat,
                                treat_lfc=treat_lfc,
                                pbps_counts=ct_pbps_counts,
                                pbps_subject_ids=ct_pbps_subject_ids,
                                sva_n=sva_n
                            )
                        else:  # DESeq2
                            results, ruv_info, ct_dir = run_dea_for_celltype(
                                ct_df, ct_groups, ct_name, "DESeq2",
                                ruv_method, ruv_k, ruv_alpha, res_dir,
                                deseq2_shrinkage=deseq2_shrinkage,
                                ref_group=ref_group,
                                pbps_counts=ct_pbps_counts,
                                pbps_subject_ids=ct_pbps_subject_ids,
                                sva_n=sva_n
                            )

                        all_results[ct_name] = results
                        all_ruv_info[ct_name] = ruv_info
                        all_ct_dirs[ct_name] = ct_dir

                        # Display results
                        for comp_name, df_res in results.items():
                            st.markdown(f"**{comp_name}**")
                            n_sig = (df_res.get('FDR', df_res.get('padj', pd.Series())) < 0.05).sum()
                            st.write(f"  Significant genes (FDR < 0.05): {n_sig}")
                            st.dataframe(df_res.head(10))

                        # Display QC plots
                        for plot_name in ['RLE_before_DEA.png', 'PCA_before_DEA.png',
                                          'BCV_plot.png', 'DispersionEstimates.png']:
                            plot_path = os.path.join(ct_dir, plot_name)
                            if os.path.exists(plot_path):
                                st.image(plot_path, caption=f"{ct_name}: {plot_name}",
                                         use_container_width=True)

                        # Batch correction QC plots
                        if ruv_method != "None":
                            if ruv_method == "SVA":
                                factor_label = "SV (Surrogate Variable)"
                            else:
                                factor_label = "W (Unwanted Factor)"
                            st.info(
                                f"**{factor_label} Stripchart interpretation:** "
                                f"If each {factor_label.split(' ')[0]} **overlaps** between groups, "
                                f"it is properly capturing technical variation (good). "
                                f"If it **separates** between groups, it is confounded with the biological condition, "
                                f"and including it in the design matrix may also remove biological signal, "
                                f"reducing detection power. "
                                f"If confounding is observed, consider reducing k or changing the method."
                            )
                            ruv_plots = []
                            if ruv_method == "SVA":
                                n_sv_used = ruv_info.get('n_sv', 2)
                                ruv_plots = [f'SVA_stripchart_n{n_sv_used}.png',
                                             f'PCA_SVA_n{n_sv_used}.png']
                            elif "RUVg" in ruv_method:
                                ruv_plots = [f'W_stripchart_RUVg_k{ruv_k}.png',
                                             f'RLE_after_RUVg_k{ruv_k}.png',
                                             f'PCA_after_RUVg_k{ruv_k}.png']
                            elif "RUV2" in ruv_method:
                                ruv_plots = [f'W_stripchart_RUV2_k{ruv_k}.png',
                                             f'RLE_after_RUV2_k{ruv_k}.png',
                                             f'PCA_after_RUV2_k{ruv_k}.png']
                            elif ruv_method == "RUVIII PBPS":
                                ruv_plots = [f'W_stripchart_RUVIII_PBPS_k{ruv_k}.png',
                                             f'RLE_RUVIII_PBPS_k{ruv_k}.png',
                                             f'PCA_RUVIII_PBPS_orig_k{ruv_k}.png',
                                             f'PCA_RUVIII_PBPS_combined_k{ruv_k}.png']
                            elif "RUVIII" in ruv_method:
                                ruv_plots = [f'W_stripchart_RUVIII_k{ruv_k}.png',
                                             f'RLE_after_RUVIII_k{ruv_k}.png',
                                             f'PCA_after_RUVIII_k{ruv_k}.png']
                            for plot_name in ruv_plots:
                                plot_path = os.path.join(ct_dir, plot_name)
                                if os.path.exists(plot_path):
                                    st.image(plot_path,
                                             caption=f"{ct_name}: {plot_name}",
                                             use_container_width=True)

                            if ruv_info:
                                if 'n_sv' in ruv_info:
                                    st.write(f"  Number of SVs: {ruv_info['n_sv']} "
                                             f"(BE estimate: {ruv_info.get('n_sv_be', 'N/A')}, "
                                             f"Leek estimate: {ruv_info.get('n_sv_leek', 'N/A')})")
                                if 'n_control_genes' in ruv_info:
                                    st.write(f"  Control genes: {ruv_info['n_control_genes']}")
                                if 'n_pbps' in ruv_info:
                                    st.write(f"  PBPS samples: {ruv_info['n_pbps']}, "
                                             f"Original samples: {ruv_info['n_orig']}, "
                                             f"Genes used: {ruv_info['n_genes']}")
                                if 'method' in ruv_info:
                                    st.write(f"  Method: {ruv_info['method']}")

                        # MA plots
                        for comp_name in results.keys():
                            ma_path = os.path.join(ct_dir, f"MA_{comp_name}.png")
                            if os.path.exists(ma_path):
                                st.image(ma_path,
                                         caption=f"{ct_name}: MA plot ({comp_name})",
                                         use_container_width=True)

                except Exception as e:
                    st.error(f"  {ct_name}: Error - {str(e)}")
                    import traceback
                    st.code(traceback.format_exc())

                progress.progress((idx + 1) / total)

            # ---- Download ----
            if all_results:
                st.markdown("---")
                st.markdown("### Download")

                # Create ZIP
                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zf:
                    for ct_name, ct_dir in all_ct_dirs.items():
                        ct_clean = ct_name.replace(' ', '_').replace('/', '_')
                        for fname in os.listdir(ct_dir):
                            fpath = os.path.join(ct_dir, fname)
                            if os.path.isfile(fpath):
                                zf.write(fpath, os.path.join(ct_clean, fname))

                zip_buffer.seek(0)

                file_prefix = os.path.splitext(uploaded_file.name)[0]
                zip_name = f"{file_prefix}_pseudobulk_{dea_method}"
                if ruv_method == "SVA":
                    zip_name += f"_SVA_n{sva_n if sva_n else 'auto'}"
                elif ruv_method != "None":
                    zip_name += f"_{ruv_method.split(' ')[0]}_k{ruv_k}"
                zip_name += ".zip"

                st.download_button(
                    label="Download all results as ZIP",
                    data=zip_buffer.getvalue(),
                    file_name=zip_name,
                    mime="application/zip"
                )

                st.success("Analysis complete")

else:
    st.info("Please upload an h5ad file.")
