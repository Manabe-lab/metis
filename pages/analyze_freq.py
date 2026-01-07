import streamlit as st
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
from statsmodels.othermod.betareg import BetaModel
from statsmodels.genmod.families import Binomial
from statsmodels.genmod.generalized_linear_model import GLM
import plotly.express as px
from helper_func import remove_sample_num
import os
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

# rpy2 configuration
pandas2ri.activate()

def remove_common_suffix(strings):
    if not strings or len(strings) == 0:
        return []
    # Get the length of the shortest string
    min_length = min(len(s) for s in strings)
    # Find the length of the common suffix
    suffix_length = 0
    for i in range(1, min_length + 1):
        suffix = strings[0][-i:]
        if all(s.endswith(suffix) for s in strings):
            suffix_length = i
        else:
            break
    # If no common suffix is found, return the original list
    if suffix_length == 0:
        return strings
    # Create a new list with the common suffix removed
    return [s[:-suffix_length] for s in strings]


def read_file_with_options(file, file_extension):
    """Function to read files"""
    try:
        file_extension = file_extension.lower()

        # For Excel files
        if file_extension in ['xlsx', 'xls']:
            df = pd.read_excel(file, index_col=0)
        # For other text files (csv, tsv, txt)
        else:
            try:
                # Try tab-delimited first
                df = pd.read_csv(file, sep='\t', index_col=0)
            except:
                # Auto-detect delimiter
                df = pd.read_csv(file, sep=None, engine='python', index_col=0)

        return df

    except Exception as e:
        st.error(f'File reading error: {str(e)}')
        return None


def analyze_logit_t(data, groups, df):
    """Logit transformation + t-test"""
    # Ensure data is treated as a numpy array
    data = np.asarray(data, dtype=float)

    # Calculate column sums and convert to proportions
    column_sums = data.sum(axis=0)
    proportions = data / column_sums

    # Logit transformation
    logit_transformed = np.log(proportions / (1 - proportions))

    # Get group indices
    unique_groups = np.unique(groups)
    if len(unique_groups) != 2:
        st.error("Currently only supports comparison of 2 groups")
        return None

    group1, group2 = unique_groups
    group1_idx = np.where(groups == group1)[0]
    group2_idx = np.where(groups == group2)[0]

    # Split data by group
    group1_data = logit_transformed[:, group1_idx]
    group2_data = logit_transformed[:, group2_idx]

    # Perform t-test for each cluster
    t_stats = []
    p_values = []
    means_group1 = []
    means_group2 = []
    cell_types = df.index.tolist()

    for i in range(len(data)):
        t_stat, p_val = stats.ttest_ind(group1_data[i, :], group2_data[i, :])
        t_stats.append(t_stat)
        p_values.append(p_val)
        means_group1.append(np.mean(proportions[i, group1_idx] * 100))
        means_group2.append(np.mean(proportions[i, group2_idx] * 100))

    # Multiple comparison correction
    rejected_holm, p_adjusted_holm, _, _ = multipletests(p_values, method='holm')
    rejected_bh, p_adjusted_bh, _, _ = multipletests(p_values, method='fdr_bh')
    rejected_by, p_adjusted_by, _, _ = multipletests(p_values, method='fdr_by')

    # Compile results into a dataframe
    results = pd.DataFrame({
        'Cell_Type': cell_types,
        f'Mean_{group1}': means_group1,
        f'Mean_{group2}': means_group2,
        't_statistic': t_stats,
        'p_value': p_values,
        'p_adjusted_Holm': p_adjusted_holm,
        'p_adjusted_BH': p_adjusted_bh,
        'p_adjusted_BY': p_adjusted_by
    })

    return results


def analyze_beta_regression(data, groups, df):
    """Beta regression"""
    # Ensure data is treated as a numpy array
    data = np.asarray(data, dtype=float)

    # Calculate column sums and convert to proportions
    column_sums = data.sum(axis=0)
    proportions = data / column_sums

    unique_groups = np.unique(groups)
    if len(unique_groups) != 2:
        st.error("Currently only supports comparison of 2 groups")
        return None

    group1, group2 = unique_groups

    # Prepare data for beta regression
    cell_types_beta = []
    group_values = []
    values = []

    for i in range(len(data)):
        cell_name = df.index[i]
        for j in range(data.shape[1]):
            cell_types_beta.append(cell_name)
            group_values.append(groups[j])
            values.append(proportions[i, j])

    beta_df = pd.DataFrame({
        'CellType': cell_types_beta,
        'Group': group_values,
        'Proportion': values
    })

    # Store beta regression results
    beta_results = []

    for cell_type in np.unique(cell_types_beta):
        subset = beta_df[beta_df['CellType'] == cell_type].copy()

        # Create design matrix
        X = pd.get_dummies(subset['Group'], drop_first=True)
        X = sm.add_constant(X)
        X = np.asarray(X, dtype=float)

        y = np.asarray(subset['Proportion'], dtype=float)

        try:
            model = BetaModel(y, X)
            results_model = model.fit()

            mean_group1 = subset[subset['Group'] == group1]['Proportion'].mean() * 100
            mean_group2 = subset[subset['Group'] == group2]['Proportion'].mean() * 100

            beta_results.append({
                'Cell_Type': cell_type,
                f'Mean_{group1}': mean_group1,
                f'Mean_{group2}': mean_group2,
                'Estimate': results_model.params[1],
                'Std_Error': results_model.bse[1],
                'z_value': results_model.tvalues[1],
                'p_value': results_model.pvalues[1]
            })
        except Exception as e:
            st.warning(f"Beta regression failed for {cell_type}: {str(e)}")
            beta_results.append({
                'Cell_Type': cell_type,
                f'Mean_{group1}': np.nan,
                f'Mean_{group2}': np.nan,
                'Estimate': np.nan,
                'Std_Error': np.nan,
                'z_value': np.nan,
                'p_value': np.nan
            })

    if beta_results:
        results = pd.DataFrame(beta_results)

        # FDR correction
        valid_pvalues = results['p_value'].dropna()
        if len(valid_pvalues) > 0:
            _, p_adjusted_bh_valid = multipletests(valid_pvalues, method='fdr_bh')[:2]
            _, p_adjusted_by_valid = multipletests(valid_pvalues, method='fdr_by')[:2]

            results['p_adjusted_BH'] = np.nan
            results['p_adjusted_BY'] = np.nan
            results.loc[results['p_value'].notna(), 'p_adjusted_BH'] = p_adjusted_bh_valid
            results.loc[results['p_value'].notna(), 'p_adjusted_BY'] = p_adjusted_by_valid

    return results


def analyze_binomial_glm(data, groups, df, ref_group=None, use_robust_se=True, se_type='HC3'):
    """Binomial GLM (logistic regression) with robust standard errors"""
    # Ensure data is treated as a numpy array
    data = np.asarray(data, dtype=float)

    # Calculate column sums
    column_sums = data.sum(axis=0)

    unique_groups = np.unique(groups)
    if len(unique_groups) != 2:
        st.error("Currently only supports comparison of 2 groups")
        return None

    group1, group2 = unique_groups

    # Set reference group
    if ref_group is None:
        ref_group = group1
    elif ref_group not in unique_groups:
        st.error(f"Reference group '{ref_group}' does not exist in the data")
        return None

    # Determine comparison group
    comparison_group = group2 if ref_group == group1 else group1

    # Store results
    binomial_results = []

    for i in range(len(data)):
        cell_name = df.index[i]

        # Prepare count data (ensure numeric type)
        success_counts = np.asarray(data[i, :], dtype=float)
        total_counts = np.asarray(column_sums, dtype=float)
        failure_counts = total_counts - success_counts

        # Create dataframe (code reference group as 0)
        group_coded = np.where(groups == ref_group, 0, 1)

        glm_df = pd.DataFrame({
            'success': success_counts,
            'failure': failure_counts,
            'group_coded': group_coded
        })

        # Design matrix (ensure float type)
        X = sm.add_constant(glm_df['group_coded'], has_constant='add')
        X = np.asarray(X, dtype=float)

        # Response variable (pairs of success and failure counts)
        y = np.column_stack([
            np.asarray(glm_df['success'], dtype=float),
            np.asarray(glm_df['failure'], dtype=float)
        ])

        try:
            # Binomial GLM
            model = GLM(y, X, family=Binomial())
            result = model.fit()

            # Calculate standard errors
            if use_robust_se and se_type in ['HC0', 'HC1', 'HC2', 'HC3']:
                try:
                    # Get HC robust standard errors
                    result_robust = result.get_robustcov_results(cov_type=se_type)
                    std_error = result_robust.bse[1]
                    z_value = result_robust.tvalues[1]
                    p_value = result_robust.pvalues[1]
                    se_label = f'Std_Error_{se_type}'
                except:
                    # Fallback: normal standard errors
                    std_error = result.bse[1]
                    z_value = result.tvalues[1]
                    p_value = result.pvalues[1]
                    se_label = 'Std_Error'
            elif use_robust_se and se_type == 'quasi':
                # Quasi-binomial (estimate overdispersion parameter φ)
                # Estimate overdispersion parameter from residuals
                pearson_chi2 = result.pearson_chi2
                df_resid = result.df_resid
                dispersion = pearson_chi2 / df_resid if df_resid > 0 else 1.0

                # Correct for overdispersion
                std_error = result.bse[1] * np.sqrt(dispersion)
                z_value = result.params[1] / std_error
                p_value = 2 * (1 - stats.norm.cdf(abs(z_value)))
                se_label = 'Std_Error_quasi'
            else:
                # Normal standard errors
                std_error = result.bse[1]
                z_value = result.tvalues[1]
                p_value = result.pvalues[1]
                se_label = 'Std_Error'

            # Mean proportion for each group
            ref_idx = groups == ref_group
            comp_idx = groups == comparison_group
            mean_ref = (success_counts[ref_idx] / total_counts[ref_idx]).mean() * 100
            mean_comp = (success_counts[comp_idx] / total_counts[comp_idx]).mean() * 100

            # Odds ratio and CI
            estimate = result.params[1]
            odds_ratio = np.exp(estimate)

            # 95% confidence interval (using robust SE)
            z_critical = 1.96
            ci_low = np.exp(estimate - z_critical * std_error)
            ci_high = np.exp(estimate + z_critical * std_error)

            # Fitted values
            fitted_ref = result.fittedvalues[ref_idx].mean()
            fitted_comp = result.fittedvalues[comp_idx].mean()

            # Store results in reference group order
            if ref_group == group1:
                mean_g1 = mean_ref
                mean_g2 = mean_comp
            else:
                mean_g1 = mean_comp
                mean_g2 = mean_ref

            result_dict = {
                'Cell_Type': cell_name,
                f'Mean_{group1}': mean_g1,
                f'Mean_{group2}': mean_g2,
                'Comparison': f'{comparison_group} vs {ref_group}',
                'Estimate_logOR': estimate,
                'Odds_Ratio': odds_ratio,
                'OR_CI95_low': ci_low,
                'OR_CI95_high': ci_high,
                'z_value': z_value,
                'p_value': p_value,
                'fitted_p_ref': fitted_ref,
                'fitted_p_comp': fitted_comp,
                'converged': result.converged
            }
            result_dict[se_label] = std_error

            # Add overdispersion parameter (for quasi-binomial)
            if se_type == 'quasi':
                result_dict['dispersion'] = dispersion

            binomial_results.append(result_dict)
        except Exception as e:
            st.warning(f"Binomial GLM failed for {cell_name}: {str(e)}")
            binomial_results.append({
                'Cell_Type': cell_name,
                f'Mean_{group1}': np.nan,
                f'Mean_{group2}': np.nan,
                'Comparison': 'Failed',
                'Estimate_logOR': np.nan,
                'Odds_Ratio': np.nan,
                'OR_CI95_low': np.nan,
                'OR_CI95_high': np.nan,
                'Std_Error_HC3' if use_robust_se else 'Std_Error': np.nan,
                'z_value': np.nan,
                'p_value': np.nan,
                'fitted_p_ref': np.nan,
                'fitted_p_comp': np.nan,
                'converged': False
            })

    if binomial_results:
        results = pd.DataFrame(binomial_results)

        # FDR correction
        valid_pvalues = results['p_value'].dropna()
        if len(valid_pvalues) > 0:
            _, p_adjusted_bh_valid = multipletests(valid_pvalues, method='fdr_bh')[:2]
            _, p_adjusted_by_valid = multipletests(valid_pvalues, method='fdr_by')[:2]

            results['p_adjusted_BH'] = np.nan
            results['p_adjusted_BY'] = np.nan
            results.loc[results['p_value'].notna(), 'p_adjusted_BH'] = p_adjusted_bh_valid
            results.loc[results['p_value'].notna(), 'p_adjusted_BY'] = p_adjusted_by_valid

    return results


def test_overdispersion(data, groups, df):
    """
    Test for overdispersion
    For each cluster, test whether observed variance exceeds the theoretical variance of binomial distribution

    Variance Ratio = observed variance / theoretical variance
    - Close to 1: binomial distribution is sufficient
    - > 1.5: possible overdispersion (Beta-binomial recommended)
    - > 2.0: clear overdispersion (Beta-binomial required)
    """
    # Ensure data is treated as a numpy array
    data = np.asarray(data, dtype=float)

    column_sums = data.sum(axis=0)
    proportions = data / column_sums

    unique_groups = np.unique(groups)

    results = []

    for i in range(len(data)):
        cell_name = df.index[i]

        # Overall test
        obs_counts = data[i, :]
        obs_totals = column_sums
        obs_props = obs_counts / obs_totals

        mean_prop = obs_props.mean()
        obs_var = obs_props.var(ddof=1)

        # Theoretical variance (binomial distribution): Var(p) = p(1-p)/n
        theoretical_vars = mean_prop * (1 - mean_prop) / obs_totals
        theoretical_var = theoretical_vars.mean()

        if theoretical_var > 0:
            variance_ratio = obs_var / theoretical_var
        else:
            variance_ratio = np.nan

        # Classification
        if np.isnan(variance_ratio):
            status = "N/A"
        elif variance_ratio > 2.0:
            status = "Warning: Clear overdispersion"
        elif variance_ratio > 1.5:
            status = "Warning: Possible overdispersion"
        else:
            status = "OK: No issues"

        results.append({
            'Cluster': cell_name,
            'Mean_proportion': mean_prop * 100,
            'Observed_variance': obs_var,
            'Theoretical_variance': theoretical_var,
            'Variance_ratio': variance_ratio,
            'Status': status
        })

    results_df = pd.DataFrame(results)

    # Count presence of overdispersion
    n_overdispersed = sum(results_df['Variance_ratio'] > 1.5)
    n_severe = sum(results_df['Variance_ratio'] > 2.0)

    return results_df, n_overdispersed, n_severe


def analyze_propeller(data, groups, df, ref_group=None, transform='logit'):
    """
    Propeller method via R (speckle package)

    Parameters:
    -----------
    data : array-like
        Cell type count matrix (cell types x samples)
    groups : array-like
        Group labels for each sample
    df : pandas.DataFrame
        Original dataframe with cell type names as index
    ref_group : str, optional
        Reference group for comparison
    transform : str
        Transformation to use: 'asin' or 'logit'

    Returns:
    --------
    results : pandas.DataFrame
        Results with statistics for each cell type
    """
    try:
        # Load R packages
        ro.r('suppressPackageStartupMessages(library(speckle))')
        ro.r('suppressPackageStartupMessages(library(limma))')

        # Ensure data is treated as a numpy array
        data = np.asarray(data, dtype=float)

        unique_groups = np.unique(groups)
        if len(unique_groups) != 2:
            st.error("Currently only supports comparison of 2 groups")
            return None

        group1, group2 = unique_groups

        # Set reference group
        if ref_group is None:
            ref_group = group1
        elif ref_group not in unique_groups:
            st.error(f"Reference group '{ref_group}' does not exist in the data")
            return None

        # Determine comparison group
        comparison_group = group2 if ref_group == group1 else group1

        # Prepare data for R
        # Propeller receives count data per sample as DataFrame
        # Rows: cell types, columns: samples

        cell_type_names = df.index.tolist()
        sample_names = [f'Sample_{i}' for i in range(data.shape[1])]

        # Pass count data to R manually (avoid pandas2ri compatibility issues)
        # First pass as matrix to R
        data_matrix = ro.r.matrix(ro.FloatVector(data.flatten(order='F')),
                                   nrow=len(cell_type_names),
                                   ncol=len(sample_names))

        # Convert to dataframe
        ro.globalenv['data_matrix'] = data_matrix
        ro.globalenv['cell_type_names'] = ro.StrVector(cell_type_names)
        ro.globalenv['sample_names'] = ro.StrVector(sample_names)
        ro.globalenv['group_vector'] = ro.StrVector(groups.tolist())
        ro.globalenv['ref_group'] = ref_group
        ro.globalenv['transform_type'] = transform

        # Create DataFrame in R
        ro.r('''
            # Create count data DataFrame
            counts <- as.data.frame(data_matrix)
            rownames(counts) <- cell_type_names
            colnames(counts) <- sample_names

            # Create sample information DataFrame
            sample_info <- data.frame(
                sample = sample_names,
                group = group_vector,
                stringsAsFactors = FALSE
            )

            # === Check data structure ===
            cat("\n=== Data structure on R side ===\n")
            cat("\n1. counts (count data):\n")
            cat("  class:", class(counts), "\n")
            cat("  dim:", paste(dim(counts), collapse=" x "), "\n")
            cat("  rownames:", paste(head(rownames(counts), 3), "..."), "\n")
            cat("  colnames:", paste(head(colnames(counts), 3), "..."), "\n")
            cat("  Data type:\n")
            print(str(counts))
            cat("\n  First rows:\n")
            print(head(counts[, 1:min(3, ncol(counts))], 3))

            cat("\n2. sample_names:\n")
            cat("  class:", class(sample_names), "\n")
            cat("  length:", length(sample_names), "\n")
            cat("  content:", paste(head(sample_names, 5), collapse=", "), "\n")

            cat("\n3. group_vector:\n")
            cat("  class:", class(group_vector), "\n")
            cat("  length:", length(group_vector), "\n")
            cat("  content:", paste(head(group_vector, 10), collapse=", "), "\n")
            cat("  unique:", paste(unique(group_vector), collapse=", "), "\n")

            cat("\n4. ref_group:\n")
            cat("  class:", class(ref_group), "\n")
            cat("  content:", ref_group, "\n")

            cat("\n5. transform_type:\n")
            cat("  class:", class(transform_type), "\n")
            cat("  content:", transform_type, "\n")
            cat("\n===================\n\n")
        ''')

        # Execute Propeller analysis
        try:
            ro.r('''
                library(speckle)

                cat("\n=== Propeller execution: expand count matrix to cell-level data ===\n")

                # Expand count matrix to cell-level vectors
                # propeller expects cell-level data (cluster, sample, group for each cell)
                clusters_vec <- c()
                sample_vec <- c()
                group_vec <- c()

                # Loop through each sample and each cell type
                for(j in 1:ncol(counts)) {
                    for(i in 1:nrow(counts)) {
                        n_cells <- counts[i, j]
                        if(n_cells > 0) {
                            # Add n_cells entries
                            clusters_vec <- c(clusters_vec, rep(rownames(counts)[i], n_cells))
                            sample_vec <- c(sample_vec, rep(colnames(counts)[j], n_cells))
                            group_vec <- c(group_vec, rep(group_vector[j], n_cells))
                        }
                    }
                }

                cat("Total cells expanded:", length(clusters_vec), "\n")
                cat("Unique clusters:", paste(unique(clusters_vec), collapse=", "), "\n")
                cat("Unique samples:", paste(unique(sample_vec), collapse=", "), "\n")
                cat("Unique groups:", paste(unique(group_vec), collapse=", "), "\n")

                # Convert to factors (reference group as first level)
                clusters_factor <- factor(clusters_vec)
                sample_factor <- factor(sample_vec)
                group_factor <- factor(group_vec, levels = c(ref_group, setdiff(unique(group_vec), ref_group)))

                # Execute Propeller
                props <- propeller(
                    clusters = clusters_factor,
                    sample = sample_factor,
                    group = group_factor,
                    transform = transform_type,
                    trend = FALSE,
                    robust = TRUE
                )

                cat("Propeller completed successfully!\n")
                cat("Props class:", class(props), "\n")
                cat("Props columns:", paste(colnames(props), collapse=", "), "\n")
            ''')
        except Exception as e:
            st.error(f"Error during Propeller execution: {str(e)}")
            # Display more detailed error information
            import traceback
            st.code(traceback.format_exc())
            raise

        # Retrieve results (manually convert to pandas)
        ro.r('''
            # Prepare results in list format
            # Propeller output columns: BaselineProp.clusters, PropMean.X, PropMean.Y, PropRatio, Tstatistic, P.Value, FDR
            result_list <- list(
                cell_types = as.character(props$BaselineProp.clusters),
                baseline_prop = props$BaselineProp.Freq,
                prop_mean_ctrl = props[[paste0("PropMean.", levels(group_factor)[1])]],
                prop_mean_tac = props[[paste0("PropMean.", levels(group_factor)[2])]],
                prop_ratio = props$PropRatio,
                t_statistics = props$Tstatistic,
                p_values = props$P.Value,
                fdr = props$FDR
            )
        ''')

        # Get results from R
        result_cell_types = list(ro.r('result_list$cell_types'))
        result_baseline_prop = list(ro.r('result_list$baseline_prop'))
        result_prop_mean_ctrl = list(ro.r('result_list$prop_mean_ctrl'))
        result_prop_mean_tac = list(ro.r('result_list$prop_mean_tac'))
        result_prop_ratio = list(ro.r('result_list$prop_ratio'))
        result_t_stats = list(ro.r('result_list$t_statistics'))
        result_p_values = list(ro.r('result_list$p_values'))
        result_fdr = list(ro.r('result_list$fdr'))

        # Organize column names
        results = pd.DataFrame()
        results['Cell_Type'] = result_cell_types

        # Add statistics from Propeller results
        results['BaselineProp'] = result_baseline_prop

        # PropMean is already in percentage (0-100)
        # Check group order and place appropriately
        if ref_group == group1:
            results[f'Mean_{group1}'] = result_prop_mean_ctrl
            results[f'Mean_{group2}'] = result_prop_mean_tac
        else:
            results[f'Mean_{group1}'] = result_prop_mean_tac
            results[f'Mean_{group2}'] = result_prop_mean_ctrl

        results['Comparison'] = f'{comparison_group} vs {ref_group}'
        results['PropRatio'] = result_prop_ratio
        results['t_statistic'] = result_t_stats
        results['p_value'] = result_p_values
        results['p_adjusted_BH'] = result_fdr

        # Sort by cell type name
        return results.sort_values('Cell_Type').reset_index(drop=True)

    except Exception as e:
        st.error(f"Propeller analysis error: {str(e)}")
        st.info("The speckle package may not be installed. Run install.packages('speckle') in R.")
        import traceback
        st.code(traceback.format_exc())
        return None


def analyze_beta_binomial(data, groups, df, ref_group=None):
    """Beta-binomial regression via R (VGAM package)"""
    try:
        # Load R package
        ro.r('suppressPackageStartupMessages(library(VGAM))')

        # Ensure data is treated as a numpy array
        data = np.asarray(data, dtype=float)

        unique_groups = np.unique(groups)
        if len(unique_groups) != 2:
            st.error("Currently only supports comparison of 2 groups")
            return None

        group1, group2 = unique_groups

        # Set reference group
        if ref_group is None:
            ref_group = group1
        elif ref_group not in unique_groups:
            st.error(f"Reference group '{ref_group}' does not exist in the data")
            return None

        # Calculate column sums
        column_sums = data.sum(axis=0)

        # Store results
        bb_results = []

        for i in range(len(data)):
            cell_name = df.index[i]

            # Prepare count data (ensure integer type)
            success_counts = np.asarray(data[i, :], dtype=int)
            total_counts = np.asarray(column_sums, dtype=int)

            # Pass data to R
            ro.globalenv['success'] = ro.IntVector(success_counts.tolist())
            ro.globalenv['total'] = ro.IntVector(total_counts.tolist())
            ro.globalenv['group'] = ro.StrVector(groups.tolist())
            ro.globalenv['ref_group'] = ref_group

            try:
                # Beta-binomial regression
                ro.r('''
                    # Create dataframe
                    df_r <- data.frame(
                        success = success,
                        total = total,
                        group = factor(group)
                    )

                    # Fix reference group
                    df_r$group <- relevel(df_r$group, ref = ref_group)

                    # Null model (intercept only) for LRT
                    model_null <- vglm(cbind(success, total - success) ~ 1,
                                      family = betabinomial,
                                      data = df_r)

                    # Full model (with group)
                    model_full <- vglm(cbind(success, total - success) ~ group,
                                      family = betabinomial,
                                      data = df_r)

                    # Wald test results from full model
                    cs <- coef(summary(model_full))

                    # Get group coefficient row name (rows starting with "group")
                    group_coef_name <- rownames(cs)[grep("^group", rownames(cs))]

                    if (length(group_coef_name) > 0) {
                        est <- cs[group_coef_name, "Estimate"]
                        se  <- cs[group_coef_name, "Std. Error"]
                        z_wald <- cs[group_coef_name, "z value"]
                        p_wald <- cs[group_coef_name, "Pr(>|z|)"]
                    } else {
                        est <- NA
                        se  <- NA
                        z_wald <- NA
                        p_wald <- NA
                    }

                    # Likelihood Ratio Test (LRT)
                    tryCatch({
                        lrt_stat <- 2 * (logLik(model_full) - logLik(model_null))
                        lrt_df <- 1  # 1 parameter difference (group coefficient)
                        p_lrt <- pchisq(lrt_stat, df = lrt_df, lower.tail = FALSE)
                    }, error = function(e) {
                        lrt_stat <- NA
                        p_lrt <- NA
                    })

                    # Get overdispersion parameter (rho)
                    tryCatch({
                        rho <- Coef(model_full)["loge(rho)"]
                        rho <- exp(rho)  # log(rho) so convert back with exp
                    }, error = function(e) {
                        rho <- NA
                    })

                    # Mean proportion for each group
                    mean_ref <- mean(df_r$success[df_r$group == ref_group] /
                                    df_r$total[df_r$group == ref_group]) * 100
                    mean_other <- mean(df_r$success[df_r$group != ref_group] /
                                      df_r$total[df_r$group != ref_group]) * 100

                    # Odds ratio (exp(estimate))
                    odds_ratio <- exp(est)

                    # Debug information
                    group_levels <- levels(df_r$group)
                    comparison_name <- paste(group_levels[2], "vs", group_levels[1])
                ''')

                # Get results from R
                estimate = float(ro.r('est')[0]) if not pd.isna(ro.r('est')[0]) else np.nan
                std_error = float(ro.r('se')[0]) if not pd.isna(ro.r('se')[0]) else np.nan
                z_wald = float(ro.r('z_wald')[0]) if not pd.isna(ro.r('z_wald')[0]) else np.nan
                p_wald = float(ro.r('p_wald')[0]) if not pd.isna(ro.r('p_wald')[0]) else np.nan
                lrt_stat = float(ro.r('lrt_stat')[0]) if not pd.isna(ro.r('lrt_stat')[0]) else np.nan
                p_lrt = float(ro.r('p_lrt')[0]) if not pd.isna(ro.r('p_lrt')[0]) else np.nan
                rho_est = float(ro.r('rho')[0]) if not pd.isna(ro.r('rho')[0]) else np.nan
                odds_ratio = float(ro.r('odds_ratio')[0]) if not pd.isna(ro.r('odds_ratio')[0]) else np.nan

                mean_ref_val = float(ro.r('mean_ref')[0])
                mean_other_val = float(ro.r('mean_other')[0])
                comparison_name = str(ro.r('comparison_name')[0])

                # Place reference and other appropriately
                if ref_group == group1:
                    mean_g1 = mean_ref_val
                    mean_g2 = mean_other_val
                else:
                    mean_g1 = mean_other_val
                    mean_g2 = mean_ref_val

                bb_results.append({
                    'Cell_Type': cell_name,
                    f'Mean_{group1}': mean_g1,
                    f'Mean_{group2}': mean_g2,
                    'Comparison': comparison_name,
                    'Estimate_logOR': estimate,
                    'Odds_Ratio': odds_ratio,
                    'Std_Error': std_error,
                    'z_wald': z_wald,
                    'p_value_Wald': p_wald,
                    'LRT_stat': lrt_stat,
                    'p_value_LRT': p_lrt,
                    'rho': rho_est
                })

            except Exception as e:
                st.warning(f"Beta-binomial regression failed for {cell_name}: {str(e)}")
                bb_results.append({
                    'Cell_Type': cell_name,
                    f'Mean_{group1}': np.nan,
                    f'Mean_{group2}': np.nan,
                    'Comparison': 'Failed',
                    'Estimate_logOR': np.nan,
                    'Odds_Ratio': np.nan,
                    'Std_Error': np.nan,
                    'z_wald': np.nan,
                    'p_value_Wald': np.nan,
                    'LRT_stat': np.nan,
                    'p_value_LRT': np.nan,
                    'rho': np.nan
                })

        if bb_results:
            results = pd.DataFrame(bb_results)

            # FDR correction - Wald test
            valid_pvalues_wald = results['p_value_Wald'].dropna()
            if len(valid_pvalues_wald) > 0:
                _, p_adjusted_bh_wald = multipletests(valid_pvalues_wald, method='fdr_bh')[:2]
                _, p_adjusted_by_wald = multipletests(valid_pvalues_wald, method='fdr_by')[:2]

                results['p_adjusted_BH_Wald'] = np.nan
                results['p_adjusted_BY_Wald'] = np.nan
                results.loc[results['p_value_Wald'].notna(), 'p_adjusted_BH_Wald'] = p_adjusted_bh_wald
                results.loc[results['p_value_Wald'].notna(), 'p_adjusted_BY_Wald'] = p_adjusted_by_wald

            # FDR correction - LRT
            valid_pvalues_lrt = results['p_value_LRT'].dropna()
            if len(valid_pvalues_lrt) > 0:
                _, p_adjusted_bh_lrt = multipletests(valid_pvalues_lrt, method='fdr_bh')[:2]
                _, p_adjusted_by_lrt = multipletests(valid_pvalues_lrt, method='fdr_by')[:2]

                results['p_adjusted_BH_LRT'] = np.nan
                results['p_adjusted_BY_LRT'] = np.nan
                results.loc[results['p_value_LRT'].notna(), 'p_adjusted_BH_LRT'] = p_adjusted_bh_lrt
                results.loc[results['p_value_LRT'].notna(), 'p_adjusted_BY_LRT'] = p_adjusted_by_lrt

        return results

    except Exception as e:
        st.error(f"Beta-binomial regression error: {str(e)}")
        st.info("The VGAM package may not be installed. Run install.packages('VGAM') in R.")
        return None


def main():
    st.title('Comparison of cluster/cell type proportions')

    # Analysis method selection
    st.sidebar.title("Analysis Method")
    analysis_method = st.sidebar.radio(
        "Select analysis method:",
        ["Logit-transformed t-test",
         "Beta regression",
         "Binomial GLM",
         "Beta-binomial regression",
         "Propeller (arcsin)",
         "Propeller (logit)"],
        index=3,  # Default is Beta-binomial regression
        help="""
        **Logit-transformed t-test:**
        - Data: Uses proportion data
        - Method: t-test after logit transformation
        - Feature: Simple but ignores compositional constraints
        - Recommended: As reference value

        **Beta regression:**
        - Data: Uses proportion data (0-1)
        - Method: Regression assuming Beta distribution
        - Feature: Caution needed for 0/1 boundary values
        - Recommended: As reference value

        **Binomial GLM (recommended):**
        - Data: Uses count data (cell numbers)
        - Method: Logistic regression with binomial distribution
        - Standard errors: Normal/Quasi-binomial/HC robust (selectable)
        - Feature: Select appropriate standard error based on overdispersion
        - Recommended: Appropriate for main analysis

        **Beta-binomial regression (recommended):**
        - Data: Count data required (cell numbers)
        - Method: Beta-binomial distribution accounts for overdispersion
        - Feature: Most rigorous model handling overdispersion
        - Recommended: Optimal for main analysis (requires R package VGAM)

        **Propeller (arcsin) (recommended):**
        - Data: Uses count data (cell numbers)
        - Method: arcsin square root transformation + empirical Bayes moderated test
        - Feature: Robust to outliers, best performance with large samples
        - Recommended: n≥10/group, or when outliers are a concern (requires R speckle)
        - Use case: **Optimal for scRNA-seq cell type proportion comparisons**

        **Propeller (logit) (recommended):**
        - Data: Uses count data (cell numbers)
        - Method: logit transformation + empirical Bayes moderated test
        - Feature: Superior Type I error control in small samples
        - Recommended: Small samples with n=3-5/group (requires R speckle)
        - Use case: **Optimal for scRNA-seq cell type proportion comparisons**

        💡 **For scRNA-seq data**: Also consider **scCODA** for direct analysis from h5ad files
        (in the "scRNA-seq" section of METIS)
        """
    )

    # Method selection guide
    st.info("""
    **📌 Method Selection Guide:**

    **Single cluster (one proportion) comparison between 2 groups:**
    - ✅ Logit-t / Beta regression: Practically fine (simple)
    - ✅ Binomial GLM: More accurate (recommended)
    - ✅ Beta-binomial: Optimal when overdispersion exists
    - ✅ Propeller: Stabilizes variance estimation with empirical Bayes (recommended)

    **Multiple clusters (compositional data) comparison between 2 groups:**
    - ⚠️ Logit-t / Beta regression: Inappropriate as they ignore compositional constraints (reference value only)
    - ✅ Binomial GLM: Appropriate (recommended)
    - ✅ Beta-binomial: Optimal solution accounting for overdispersion (most recommended)
    - ✅ **Propeller (most recommended)**: Shares information across cell types, stabilizing variance estimation

    **Recommendations by sample size:**
    - **Small samples (n=3-5/group)**: Propeller (logit) or Beta-binomial
    - **Medium to large samples (n≥10/group)**: Propeller (arcsin) or Beta-binomial
    - **When outliers are a concern**: Propeller (arcsin) is robust

    💡 **When total cell counts vary greatly between samples**, count-based methods (Binomial GLM/Beta-binomial/Propeller) are required.

    ---

    **🔬 About scRNA-seq data:**
    - **Propeller**: A method designed for group comparisons of cell type proportions in single-cell data
    - **scCODA**: For compositional analysis of scRNA-seq data, **scCODA (Compositional Data Analysis)** is also available
      - 📊 **scCODA app**: Available in the "scRNA-seq" section of METIS
      - ✨ Hierarchical Bayesian model that appropriately handles compositional constraints
      - ✨ Direct analysis from h5ad files
    """)

    # Standard error selection explanation
    with st.expander("🔧 About Standard Error Selection (Important)"):
        st.markdown("""
        ## Standard Error Selection: An Important Decision

        ### ⚠️ Important Conclusions

        **Standard error selection in binomial GLM:**

        1. **No overdispersion** → **Normal standard error** (most efficient)
        2. **Overdispersion present** → **Quasi-binomial** (theoretically optimal)
        3. **More rigorous** → **Beta-binomial regression** (separate method)

        **HC3 robust standard errors:**
        - ❌ **Theoretically inappropriate** for binomial GLM
        - ⚠️ May be overly conservative
        - ✅ HC0/HC1 are acceptable (but Quasi-binomial is better)

        ---

        ## Details of Each Standard Error

        ### 1. **Normal Standard Error** (recommended: no overdispersion)

        **Theory:**
        - Standard inference for binomial GLM
        - Variance structure: Var(Y) = np(1-p) built-in

        **Usage conditions:**
        - Variance ratio < 1.5 in overdispersion test
        - Binomial distribution assumption is valid

        **Advantages:**
        - Most efficient (high power)
        - Theoretically correct (when no overdispersion)

        ---

        ### 2. **Quasi-binomial** (recommended: overdispersion present)

        **Theory:**
        - Estimates overdispersion parameter φ
        - Var(Y) = φ × np(1-p)
        - Corrects standard error by √φ

        **Usage conditions:**
        - Variance ratio ≥ 1.5 in overdispersion test
        - Variation larger than binomial distribution variance

        **Advantages:**
        - Natural extension of binomial GLM
        - **Theoretically most appropriate** (when overdispersion exists)
        - Simpler calculation than Beta-binomial

        **Reporting example in literature:**
        ```
        "We used binomial generalized linear models with
        quasi-binomial variance to account for overdispersion
        (dispersion parameter φ = X.XX)."
        ```

        ---

        ### 3. **HC Robust Standard Errors** (⚠️ Not recommended for binomial GLM)

        **HC series:**
        - **HC0**: White standard errors
        - **HC1**: HC0 × √(n/(n-k)) - small sample correction
        - **HC2**: Accounts for leverage
        - **HC3**: Most conservative

        **Original use:**
        - Heteroscedasticity in linear regression
        - Overdispersion in Poisson regression

        **Problems in binomial GLM:**
        1. Binomial GLM already models variance structure
           - Var(p) = p(1-p)/n is built-in
        2. If overdispersion exists, **Quasi-binomial is theoretically correct**
        3. HC may be overly conservative
           - Especially HC3 is extreme in small samples

        **Cases to use:**
        - ⚠️ Only in special situations (strong influence of outliers, etc.)
        - ⚠️ HC0/HC1 are "still" acceptable
        - ❌ HC3 not recommended

        ---

        ## Recommended Flow (Important)

        ```
        1. Run overdispersion test
           ├─ Variance ratio < 1.5
           │  └→ Normal standard error
           │
           └─ Variance ratio ≥ 1.5 (overdispersion present)
              ├→ Quasi-binomial (recommended)
              │
              └→ More rigorous inference needed
                 └→ Beta-binomial regression
        ```

        **HC robust standard errors:**
        - Usually unnecessary for binomial GLM
        - Quasi-binomial is theoretically more appropriate

        ---

        ## Comparison Table

        | Method | Overdispersion | Theoretical Appropriateness | Recommendation for Binomial GLM |
        |------|-----------|--------------|-----------------|
        | Normal SE | ❌ | ✅✅✅ | ⭐⭐⭐ (no OD) |
        | Quasi-binomial | ✅ | ✅✅✅✅ | ⭐⭐⭐⭐ (with OD) |
        | Beta-binomial | ✅ | ✅✅✅✅✅ | ⭐⭐⭐⭐⭐ (most rigorous) |
        | HC0/HC1 | △ | ⚠️ | ⭐ (not recommended) |
        | HC2/HC3 | △ | ❌ | ❌ (inappropriate) |

        OD = overdispersion

        ---

        ## Why Experts Accept HC0/HC1

        **Reasons:**
        1. HC0/HC1 are relatively moderate (not as extreme as HC3)
        2. Widely used as general approach to heteroscedasticity
        3. Practical convenience (software availability)

        **However:**
        - Quasi-binomial is theoretically more correct
        - Should use Quasi-binomial if overdispersion exists

        ---

        ## References

        - McCullagh, P., & Nelder, J. A. (1989). *Generalized Linear Models* (2nd ed.).
          - Chapter 4: Quasi-likelihood and overdispersion
        - Ver Hoef, J. M., & Boveng, P. L. (2007). Quasi-Poisson vs. negative binomial regression.
          - Comparison of overdispersion handling
        - Zuur, A. F., et al. (2009). *Mixed Effects Models and Extensions in Ecology with R*.
          - Practical GLM selection

        ---

        ## Conclusion

        **In this app:**

        1. Default: **Quasi-binomial**
           - Handles overdispersion
           - Theoretically appropriate
           - Safe choice

        2. No overdispersion: **Normal standard error**
           - Most efficient

        3. Most rigorous: **Beta-binomial regression**
           - Provided as separate analysis method

        **HC robust standard errors:**
        - Provided as option but not recommended
        - Use HC0/HC1 if using at all
        - Avoid HC3
        """)

    # Data format explanation
    with st.expander("📊 About Data Format"):
        st.markdown("""
        #### Input Data Format

        Recommended data format for each method:

        | Method | Data Format | Example |
        |------|-----------|-----|
        | Logit-t | Proportion or count | Either (internally converted) |
        | Beta regression | Proportion or count | Either (internally converted) |
        | **Binomial GLM** | **Count** | Cell numbers for each cluster |
        | **Beta-binomial** | **Count required** | Cell numbers for each cluster |

        #### File Format

        Rows: Cluster/cell type names
        Columns: Sample names
        Values: Count numbers (cell numbers)

        ```
        |          | Sample1 | Sample2 | Sample3 |
        |----------|---------|---------|---------|
        | Cluster0 | 150     | 200     | 180     |
        | Cluster1 | 300     | 250     | 280     |
        | Cluster2 | 100     | 150     | 120     |
        ```

        ⚠️ **Important**: Binomial GLM and Beta-binomial **must use count data**.
        Using proportion data will produce different results.
        """)

    st.markdown(f'### Current method: {analysis_method}')

    uploaded_file = st.file_uploader(
        "File upload (txt, csv, tsv, excel)",
        type=['txt', 'csv', 'tsv', 'xlsx', 'xls']
    )

    transout = st.checkbox('Transpose the data?')

    if uploaded_file is not None:
        try:
            file_extension = uploaded_file.name.split('.')[-1]
            df = read_file_with_options(uploaded_file, file_extension)

            if df is not None:
                if transout:
                    df = df.T
                    st.markdown("#### Data are transposed.")

                st.subheader("Uploaded data")
                st.dataframe(df)

                condition = [str(i) for i in df.columns.tolist()]
                group_condition = remove_common_suffix(condition)
                group_condition = [remove_sample_num(x) for x in group_condition]
                group_condition = [x.replace('_', '.') for x in group_condition]

                # Create dataframe for group settings
                group_df = pd.DataFrame({
                    'Sample': df.columns,
                    'Group': group_condition
                })

                with st.form("Set_groups"):
                    edited_group_df = st.data_editor(
                        group_df,
                        column_config={
                            "Sample": st.column_config.TextColumn("Sample", disabled=True),
                            "Group": st.column_config.TextColumn("Group")
                        }
                    )
                    submitted = st.form_submit_button("Submit")

                # Reference group selection
                st.markdown("---")
                st.markdown("### Reference Group Selection")
                st.markdown("""
                Select which group to use as reference (baseline).
                Results will be expressed as log(OR) for "comparison group vs reference group".
                """)

                unique_groups_list = edited_group_df['Group'].unique().tolist()
                if len(unique_groups_list) == 2:
                    ref_group = st.selectbox(
                        "Select reference group (baseline):",
                        unique_groups_list,
                        index=0,
                        help="Select the group to use as baseline for comparison. Results will be expressed as 'the other group vs this reference group'."
                    )
                else:
                    ref_group = None

                # Standard error type selection (for Binomial GLM)
                if analysis_method == "Binomial GLM":
                    st.markdown("---")
                    st.markdown("### Standard Error Selection (Binomial GLM)")

                    se_option = st.radio(
                        "Standard error type:",
                        ["Normal standard error (recommended)",
                         "Quasi-binomial (when overdispersion exists)",
                         "HC0 robust standard error",
                         "HC1 robust standard error (small sample)",
                         "HC2 robust standard error",
                         "HC3 robust standard error (most conservative)"],
                        index=1,  # Default is Quasi-binomial
                        help="""
                        **Recommended:**
                        1. First run overdispersion test
                        2. No overdispersion → Normal standard error
                        3. Overdispersion present → Quasi-binomial

                        **HC robust standard errors:**
                        - Not always necessary for binomial GLM
                        - HC0/HC1: Relatively moderate
                        - HC3: May be overly conservative
                        """
                    )

                    # Map options to codes
                    if se_option == "Normal standard error (recommended)":
                        use_robust_se = False
                        se_type = 'normal'
                    elif se_option == "Quasi-binomial (when overdispersion exists)":
                        use_robust_se = True
                        se_type = 'quasi'
                    elif se_option == "HC0 robust standard error":
                        use_robust_se = True
                        se_type = 'HC0'
                    elif se_option == "HC1 robust standard error (small sample)":
                        use_robust_se = True
                        se_type = 'HC1'
                    elif se_option == "HC2 robust standard error":
                        use_robust_se = True
                        se_type = 'HC2'
                    else:  # HC3
                        use_robust_se = True
                        se_type = 'HC3'

                    with st.expander("📖 Explanation of Each Standard Error"):
                        st.markdown("""
                        ### 1. Normal Standard Error (recommended)

                        **Usage conditions:**
                        - No overdispersion
                        - Binomial distribution assumption is valid
                        - Sufficient sample size (n ≥ 5 per group)

                        **Features:**
                        - Simplest
                        - Standard inference for binomial GLM
                        - Most efficient when no overdispersion

                        ---

                        ### 2. Quasi-binomial (when overdispersion exists)

                        **Usage conditions:**
                        - Overdispersion detected in overdispersion test
                        - Variation larger than binomial distribution variance

                        **Features:**
                        - Estimates overdispersion parameter φ: `Var(Y) = φ × np(1-p)`
                        - Corrects standard error by √φ
                        - Natural extension of binomial GLM
                        - **Theoretically most appropriate** (when overdispersion exists)

                        **Recommended:**
                        - ✅ First choice when overdispersion exists
                        - ✅ Simpler calculation than Beta-binomial

                        ---

                        ### 3. HC Robust Standard Errors (HC0, HC1, HC2, HC3)

                        **Original use:**
                        - Heteroscedasticity in linear regression
                        - Overdispersion in Poisson regression

                        **Problems in binomial GLM:**
                        - Binomial GLM already models variance structure
                        - If overdispersion exists, quasi-binomial is theoretically more appropriate
                        - May be overly conservative (especially HC3)

                        **Differences between types:**
                        - **HC0**: Basic White standard errors
                        - **HC1**: HC0 × √(n/(n-k)) - small sample correction
                        - **HC2**: Accounts for leverage
                        - **HC3**: Most conservative (strongly corrects influential observations)

                        **Cases to use:**
                        - ⚠️ Special situations (outliers have strong influence)
                        - ⚠️ HC0/HC1 still acceptable
                        - ❌ HC3 overly conservative for binomial GLM

                        ---

                        ### Recommended Flow

                        ```
                        1. Run overdispersion test
                           ↓
                        2. No overdispersion
                           → Normal standard error
                           ↓
                        3. Overdispersion present (variance ratio > 1.5)
                           → Quasi-binomial
                           ↓
                        4. More rigorous inference needed
                           → Beta-binomial regression
                        ```

                        ---

                        ### Comparison: Quasi vs HC vs Beta-binomial

                        | Method | Overdispersion | Theoretical Appropriateness | Computation | Recommendation |
                        |------|--------|--------------|------|--------|
                        | Quasi-binomial | ✅ | ✅✅✅ | Simple | ⭐⭐⭐ |
                        | Beta-binomial | ✅ | ✅✅✅✅ | Complex | ⭐⭐⭐⭐ |
                        | HC0/HC1 | △ | ⚠️ | Simple | ⭐ |
                        | HC3 | △ | ❌ | Simple | ❌ |
                        | Normal SE | ❌ | ✅ (if no OD) | Simple | ⭐⭐ |

                        OD = overdispersion

                        ---

                        ### References

                        - McCullagh, P., & Nelder, J. A. (1989). *Generalized Linear Models* (2nd ed.).
                        - Ver Hoef, J. M., & Boveng, P. L. (2007). Quasi-Poisson vs. negative binomial regression.
                        - Zuur, A. F., et al. (2009). *Mixed Effects Models and Extensions in Ecology with R*.
                        """)
                else:
                    use_robust_se = False
                    se_type = 'normal'

                # Overdispersion test
                st.markdown("---")
                st.subheader("🔍 Overdispersion Test")
                st.markdown("""
                Before running the analysis, you can check if overdispersion exists in the data.
                If overdispersion is detected, we recommend using Beta-binomial regression.
                """)

                if st.button("Test for overdispersion"):
                    groups = edited_group_df['Group'].values
                    data = df.values

                    with st.spinner('Testing for overdispersion...'):
                        overdispersion_results, n_overdispersed, n_severe = test_overdispersion(data, groups, df)

                    st.subheader("Overdispersion Test Results")

                    # Summary display
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Total clusters", len(overdispersion_results))
                    with col2:
                        st.metric("Possible overdispersion", n_overdispersed,
                                 delta="Caution" if n_overdispersed > 0 else None,
                                 delta_color="inverse")
                    with col3:
                        st.metric("Clear overdispersion", n_severe,
                                 delta="Warning" if n_severe > 0 else None,
                                 delta_color="inverse")

                    # Detailed results (formatted)
                    display_df = overdispersion_results.copy()
                    display_df['Mean_proportion'] = display_df['Mean_proportion'].apply(lambda x: f"{x:.2f}%")
                    display_df['Observed_variance'] = display_df['Observed_variance'].apply(lambda x: f"{x:.6f}")
                    display_df['Theoretical_variance'] = display_df['Theoretical_variance'].apply(lambda x: f"{x:.6f}")
                    display_df['Variance_ratio'] = display_df['Variance_ratio'].apply(lambda x: f"{x:.2f}" if not pd.isna(x) else "N/A")
                    st.dataframe(display_df)

                    # Recommendation message (considering sample size)
                    n_samples_per_group = pd.Series(groups).value_counts()
                    min_n = n_samples_per_group.min()

                    if n_severe > 0:
                        st.error(f"""
                        ⚠️ **Clear overdispersion detected in {n_severe} clusters!**

                        **Recommendations (overdispersion handling required):**
                        """)

                        if min_n <= 5:
                            st.markdown("""
                            **For small samples (n≤5/group):**
                            - ✅ **Propeller (logit)** - Most recommended (stabilizes variance with empirical Bayes, optimal for small samples)
                            - ✅ **Beta-binomial regression** - Recommended (rigorously models overdispersion)
                            - ⚠️ Binomial GLM (quasi) - Usable but unstable in small samples
                            - ❌ Binomial GLM (normal SE) - Inappropriate (ignores overdispersion)
                            """)
                        else:
                            st.markdown("""
                            **For sufficient sample size (n>5/group):**
                            - ✅ **Propeller (arcsin)** - Most recommended (robust to outliers, stable variance estimation)
                            - ✅ **Beta-binomial regression** - Recommended (rigorously models overdispersion)
                            - ✅ **Propeller (logit)** - Recommended (high power)
                            - ⚠️ Binomial GLM (quasi) - Usable
                            - ❌ Binomial GLM (normal SE) - Inappropriate (ignores overdispersion)
                            """)
                    elif n_overdispersed > 0:
                        st.warning(f"""
                        ⚠️ **Possible overdispersion detected in {n_overdispersed} clusters.**

                        **Recommendations (overdispersion handling recommended):**
                        """)

                        if min_n <= 5:
                            st.markdown("""
                            **For small samples (n≤5/group):**
                            - ✅ **Propeller (logit)** - Most recommended (stable even in small samples with empirical Bayes)
                            - ✅ **Beta-binomial regression** - Recommended
                            - △ Binomial GLM (quasi) - Usable but variance estimation may be unstable
                            - ⚠️ Binomial GLM (normal SE) - Not recommended as it ignores overdispersion
                            """)
                        else:
                            st.markdown("""
                            **For sufficient sample size (n>5/group):**
                            - ✅ **Propeller (arcsin)** - Most recommended (robust to outliers)
                            - ✅ **Propeller (logit)** - Recommended (high power)
                            - ✅ **Beta-binomial regression** - Recommended
                            - ✅ Binomial GLM (quasi) - Usable (corrects overdispersion)
                            - ⚠️ Binomial GLM (normal SE) - Not recommended as it ignores overdispersion
                            """)
                    else:
                        st.success("""
                        ✅ **No overdispersion detected.**

                        **Recommendations:**
                        """)

                        if min_n <= 5:
                            st.markdown("""
                            **For small samples (n≤5/group):**
                            - ✅ **Propeller (logit)** - Most recommended (stabilizes variance estimation with empirical Bayes)
                            - ✅ Binomial GLM (normal SE) - Usable
                            - ✅ Beta-binomial regression - Usable (may be somewhat conservative)

                            💡 Propeller provides the most stable estimation in small samples with empirical Bayes.
                            """)
                        else:
                            st.markdown("""
                            **For sufficient sample size (n>5/group):**
                            - ✅ **Propeller (arcsin/logit)** - Recommended (stable variance estimation with empirical Bayes)
                            - ✅ Binomial GLM (normal SE) - Recommended (most efficient)
                            - ✅ Beta-binomial regression - Usable (similar results)
                            - ✅ Binomial GLM (quasi) - Usable (somewhat conservative)

                            💡 When no overdispersion, all methods return similar results, but Propeller has the most stable variance estimation.
                            """)


                    # Add explanation
                    with st.expander("📖 What is overdispersion?"):
                        st.markdown("""
                        **Overdispersion** refers to a state where data variance is larger than the theoretical distribution's (here, binomial distribution) variance.

                        **Interpretation of Variance Ratio:**
                        - **Around 1.0**: Binomial distribution assumption is valid → Binomial GLM or Propeller
                        - **1.5 or higher**: Possible overdispersion → Beta-binomial or Propeller recommended
                        - **2.0 or higher**: Clear overdispersion → Beta-binomial or Propeller required

                        **Causes of overdispersion:**
                        - Large variation between samples
                        - Heterogeneity within groups
                        - Measurement error
                        - Biological variation

                        **Solutions:**
                        When overdispersion is detected, the following methods are recommended:

                        1. **Beta-binomial regression**: Explicitly models overdispersion parameter (ρ)
                           - Most rigorous modeling of overdispersion
                           - Somewhat high computational cost

                        2. **Propeller (arcsin/logit)**: Stabilizes variance with empirical Bayes after transformation
                           - Shares information across cell types (especially effective when ≥3 cell types)
                           - Stable estimation even in small samples (advantage of empirical Bayes)
                           - arcsin is robust to outliers, logit has superior Type I error control in small samples
                           - Fast computation

                        3. **Binomial GLM (quasi)**: Corrects with overdispersion parameter (φ)
                           - Simplified version of Beta-binomial
                           - Usable if sample size is sufficient

                        **Priority of recommendations (with overdispersion):**
                        - Small sample (n≤5): **Propeller (logit)** > Beta-binomial > Binomial GLM (quasi)
                        - Large sample (n>5): **Propeller (arcsin)** ≈ Beta-binomial > Propeller (logit) > Binomial GLM (quasi)
                        """)

                st.markdown("---")

                if st.button("Run analysis"):
                    groups = edited_group_df['Group'].values
                    data = df.values

                    # Execute according to analysis method
                    with st.spinner(f'Running {analysis_method}...'):
                        if analysis_method == "Logit-transformed t-test":
                            results = analyze_logit_t(data, groups, df)
                        elif analysis_method == "Beta regression":
                            results = analyze_beta_regression(data, groups, df)
                        elif analysis_method == "Binomial GLM":
                            results = analyze_binomial_glm(data, groups, df, ref_group=ref_group, use_robust_se=use_robust_se, se_type=se_type)
                        elif analysis_method == "Beta-binomial regression":
                            results = analyze_beta_binomial(data, groups, df, ref_group=ref_group)
                        elif analysis_method == "Propeller (arcsin)":
                            results = analyze_propeller(data, groups, df, ref_group=ref_group, transform='asin')
                        elif analysis_method == "Propeller (logit)":
                            results = analyze_propeller(data, groups, df, ref_group=ref_group, transform='logit')

                    if results is not None:
                        group1, group2 = np.unique(groups)

                        # Display results
                        st.markdown("---")
                        st.subheader(f"📊 {analysis_method} Results")

                        # Data summary
                        col1, col2, col3 = st.columns(3)
                        with col1:
                            st.metric("Clusters analyzed", len(results))
                        with col2:
                            st.metric("Group 1", group1)
                        with col3:
                            st.metric("Group 2", group2)

                        st.markdown("---")

                        # Format for display
                        display_results = results.copy()

                        # Format numeric columns to strings
                        for col in [f'Mean_{group1}', f'Mean_{group2}']:
                            if col in display_results.columns:
                                display_results[col] = display_results[col].apply(lambda x: f"{x:.4f}" if not pd.isna(x) else "N/A")

                        if 't_statistic' in display_results.columns:
                            display_results['t_statistic'] = display_results['t_statistic'].apply(lambda x: f"{x:.4f}" if not pd.isna(x) else "N/A")

                        for col in ['Estimate', 'Std_Error', 'z_value', 'z_wald', 'Estimate_logOR', 'Odds_Ratio', 'rho',
                                   'Std_Error_HC0', 'Std_Error_HC1', 'Std_Error_HC2', 'Std_Error_HC3', 'Std_Error_quasi',
                                   'OR_CI95_low', 'OR_CI95_high', 'fitted_p_ref', 'fitted_p_comp', 'dispersion', 'LRT_stat']:
                            if col in display_results.columns:
                                display_results[col] = display_results[col].apply(lambda x: f"{x:.4f}" if not pd.isna(x) else "N/A")

                        for col in ['p_value', 'Original p-value', 'p_adjusted_Holm',
                                   'p_adjusted_BH', 'p_adjusted_BY',
                                   'Holm p-value', 'BH FDR p-value', 'BY FDR p-value',
                                   'p_value_Wald', 'p_value_LRT',
                                   'p_adjusted_BH_Wald', 'p_adjusted_BY_Wald',
                                   'p_adjusted_BH_LRT', 'p_adjusted_BY_LRT']:
                            if col in display_results.columns:
                                display_results[col] = display_results[col].apply(lambda x: f"{x:.4e}" if not pd.isna(x) else "N/A")

                        # Main results table
                        st.markdown("#### 📋 Analysis Results Table")
                        st.dataframe(display_results, use_container_width=True)

                        # Original numeric data also available for inspection
                        with st.expander("🔢 Show Original Numeric Data"):
                            st.dataframe(results, use_container_width=True)
                            st.caption("※ This is the original numeric data. This is what will be downloaded.")

                        # Display additional information for Binomial GLM
                        if analysis_method == "Binomial GLM":
                            st.markdown("---")
                            st.markdown("### 📊 Binomial GLM Additional Information")

                            # Display standard error type
                            if 'dispersion' in results.columns:
                                # Quasi-binomial
                                dispersion_values = results['dispersion'].dropna()
                                if len(dispersion_values) > 0:
                                    st.success(f"""
                                    **Standard error used:** Quasi-binomial

                                    **Overdispersion parameter (φ):**
                                    - Mean: {dispersion_values.mean():.4f}
                                    - Range: {dispersion_values.min():.4f} - {dispersion_values.max():.4f}

                                    φ ≈ 1: Binomial distribution assumption is valid
                                    φ > 1: Overdispersion present (standard error corrected by √φ)
                                    """)
                            elif any(col.startswith('Std_Error_HC') for col in results.columns):
                                # HC robust standard errors
                                hc_type = [col.replace('Std_Error_', '') for col in results.columns if col.startswith('Std_Error_HC')][0]
                                st.info(f"""
                                **Standard error used:** {hc_type} robust standard error

                                **Note:**
                                - {hc_type} is robust to heteroscedasticity
                                - For binomial GLM with overdispersion, Quasi-binomial is theoretically more appropriate
                                """)
                            else:
                                # Normal standard errors
                                st.success("""
                                **Standard error used:** Normal standard error

                                Most efficient when no overdispersion.
                                """)

                            # Comparison explanation
                            if 'Comparison' in results.columns:
                                comparison = results['Comparison'].iloc[0]
                                st.info(f"""
                                **Comparison:** {comparison}

                                **Interpretation of results:**
                                - **Estimate_logOR**: log(odds ratio) - positive means former is higher, negative means latter is higher
                                - **Odds_Ratio**: Odds ratio - greater than 1 means former is higher, less than 1 means latter is higher
                                - **OR_CI95**: 95% confidence interval of odds ratio (significant if does not include 1)
                                """)

                            # Check convergence status
                            if 'converged' in results.columns:
                                n_converged = results['converged'].sum()
                                n_total = len(results)
                                if n_converged == n_total:
                                    st.success(f"✅ All models converged successfully ({n_converged}/{n_total})")
                                else:
                                    st.warning(f"⚠️ {n_total - n_converged} models did not converge")

                        # Display additional information for Beta-binomial
                        elif analysis_method == "Beta-binomial regression" and 'rho' in results.columns:
                            st.markdown("---")
                            st.markdown("### 📊 Beta-binomial Additional Information")

                            # Summary statistics for rho
                            rho_values = results['rho'].dropna()
                            if len(rho_values) > 0:
                                st.markdown(f"""
                                **Overdispersion parameter (ρ):**
                                - Mean: {rho_values.mean():.4f}
                                - Range: {rho_values.min():.4f} - {rho_values.max():.4f}

                                ρ close to 0 → Binomial distribution is sufficient
                                ρ large → Large overdispersion (Beta-binomial needed)
                                """)

                            # Comparison explanation
                            if 'Comparison' in results.columns:
                                comparison = results['Comparison'].iloc[0]
                                st.info(f"""
                                **Comparison:** {comparison}

                                - **Estimate_logOR**: log(odds ratio) - positive means former is higher, negative means latter is higher
                                - **Odds_Ratio**: Odds ratio - greater than 1 means former is higher, less than 1 means latter is higher
                                """)

                        # Display additional information for Propeller
                        elif "Propeller" in analysis_method:
                            st.markdown("---")
                            st.markdown(f"### 📊 {analysis_method} Additional Information")

                            # Transformation method explanation
                            transform_type = "arcsin square root transformation" if "arcsin" in analysis_method else "logit transformation"
                            st.success(f"""
                            **Transformation used:** {transform_type}

                            **Propeller features:**
                            - Uses Empirical Bayes framework
                            - Shares information across cell types to stabilize variance estimation
                            - Especially effective: few replicates and ≥3 cell types
                            - Multiple testing correction with FDR (False Discovery Rate)
                            - **Optimized for scRNA-seq cell type proportion comparisons**
                            """)

                            st.info("""
                            💡 **If using scRNA-seq data:**

                            Also consider **scCODA** for direct analysis from h5ad files.
                            - METIS "scRNA-seq" section → "scCODA compositional analysis"
                            - Hierarchical Bayesian model that appropriately handles compositional constraints
                            - Automatically aggregates samples from cell-level data
                            """)

                            # Sample size information
                            n_samples_per_group = pd.Series(groups).value_counts()
                            st.info(f"""
                            **Sample size:**
                            - {n_samples_per_group.index[0]}: {n_samples_per_group.values[0]} samples
                            - {n_samples_per_group.index[1]}: {n_samples_per_group.values[1]} samples

                            **Recommended transformation:**
                            - Small sample (n=3-5): logit transformation - Good Type I error control
                            - Large sample (n≥10): arcsin transformation - Best balance of power and Type I error control
                            - Outliers present: arcsin transformation - More robust
                            """)

                            # Comparison explanation
                            if 'Comparison' in results.columns:
                                comparison = results['Comparison'].iloc[0]
                                st.info(f"""
                                **Comparison:** {comparison}

                                **Interpretation of results:**
                                - **t_statistic**: t-statistic (moderated test)
                                - **p_value**: p-value (variance stabilized with empirical Bayes)
                                - **p_adjusted_BH**: FDR-corrected p-value (Benjamini-Hochberg method)
                                """)

                        # Display number of significant results
                        st.subheader("Number of significant results (α = 0.05)")

                        # For Beta-binomial, display both p-values
                        if 'p_value_Wald' in results.columns and 'p_value_LRT' in results.columns:
                            st.write(f"**Wald test p-value:** {sum(results['p_value_Wald'] < 0.05)}")
                            st.write(f"**LRT p-value (recommended):** {sum(results['p_value_LRT'] < 0.05)}")
                        else:
                            p_col = 'p_value' if 'p_value' in results.columns else 'Original p-value'
                            if p_col in results.columns:
                                st.write(f"Original p-value: {sum(results[p_col] < 0.05)}")

                        if 'p_adjusted_Holm' in results.columns:
                            st.write(f"Holm-adjusted: {sum(results['p_adjusted_Holm'] < 0.05)}")
                        if 'p_adjusted_BH' in results.columns:
                            st.write(f"BH FDR: {sum(results['p_adjusted_BH'] < 0.05)}")
                        if 'p_adjusted_BY' in results.columns:
                            st.write(f"BY FDR: {sum(results['p_adjusted_BY'] < 0.05)}")

                        # Beta-binomial FDR
                        if 'p_adjusted_BH_Wald' in results.columns:
                            st.write(f"BH FDR (Wald): {sum(results['p_adjusted_BH_Wald'] < 0.05)}")
                        if 'p_adjusted_BH_LRT' in results.columns:
                            st.write(f"**BH FDR (LRT, recommended):** {sum(results['p_adjusted_BH_LRT'] < 0.05)}")
                        if 'p_adjusted_BY_Wald' in results.columns:
                            st.write(f"BY FDR (Wald): {sum(results['p_adjusted_BY_Wald'] < 0.05)}")
                        if 'p_adjusted_BY_LRT' in results.columns:
                            st.write(f"**BY FDR (LRT, recommended):** {sum(results['p_adjusted_BY_LRT'] < 0.05)}")

                        # Visualize by cluster
                        st.subheader("Distribution of cluster proportions")

                        # Calculate column sums and convert to proportions
                        column_sums = data.sum(axis=0)
                        proportions = data / column_sums

                        # Prepare data for plotting
                        plot_data = []
                        for i in range(len(data)):
                            cell_name = df.index[i]
                            for j in range(data.shape[1]):
                                plot_data.append({
                                    'Cluster': cell_name,
                                    'Group': groups[j],
                                    'Proportion': proportions[i,j] * 100
                                })
                        plot_df = pd.DataFrame(plot_data)

                        # Create box plot for each cluster
                        for cluster in plot_df['Cluster'].unique():
                            cluster_data = plot_df[plot_df['Cluster'] == cluster]
                            fig = px.box(
                                cluster_data,
                                x='Group',
                                y='Proportion',
                                title=f'{cluster} proportion distribution (%)',
                                points='all'
                            )

                            fig.update_layout(
                                yaxis_title='Proportion (%)',
                                showlegend=True,
                                boxmode='group'
                            )

                            st.plotly_chart(fig, use_container_width=True)

                        # Download results
                        st.subheader("Download results")

                        file_name_head = os.path.splitext(uploaded_file.name)[0]
                        method_suffix = analysis_method.replace(' ', '_').replace('-', '_')

                        # Save Excel file to in-memory buffer
                        from io import BytesIO
                        buffer = BytesIO()
                        with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                            results.to_excel(writer, sheet_name='results', index=False)

                        # Download button
                        st.download_button(
                            label="Download analysis results",
                            data=buffer.getvalue(),
                            file_name=f'{file_name_head}_{method_suffix}_results.xlsx',
                            mime='application/vnd.ms-excel'
                        )

        except Exception as e:
            st.error(f'An error occurred: {str(e)}')
            import traceback
            st.code(traceback.format_exc())

if __name__ == '__main__':
    main()
