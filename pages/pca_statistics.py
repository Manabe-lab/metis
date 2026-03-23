"""
PCA Score Statistical Analysis with General Linear Models
=========================================================
Statistical analysis tool for PCA scores - supports complex designs with multiple covariates

This tool performs statistical analysis using PC scores as the response variable
for pseudobulk or sample-level data.

Main features:
- OLS (Ordinary Least Squares) with robust standard errors
- LMM (Linear Mixed Models) with random effects
- Freedman-Lane permutation testing
- Estimated Marginal Means (EMM)

Applications:
- Association between PC scores and sex/cell type in scRNA-seq pseudobulk data
- Analysis adjusted for batch effects, age, and other covariates
- Hierarchical models accounting for donor-level variation
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
from statsmodels.regression.mixed_linear_model import MixedLM
from patsy import dmatrix
import warnings
import zipfile
import io
from pathlib import Path
from datetime import datetime
warnings.filterwarnings('ignore')

# Import helper functions for sample name processing
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from helper_func import remove_after_space, remove_sample_num

# Page configuration
st.set_page_config(page_title="PCA Statistical Analysis", page_icon="📊", layout="wide")

# Initialize session state for temp directory
if "pca_temp_dir" not in st.session_state:
    temp_dir = f"temp/pca_{int(time.time())}"
    os.makedirs(temp_dir, exist_ok=True)
    st.session_state.pca_temp_dir = temp_dir

st.title("📊 PCA Score Statistical Analysis")
st.markdown("""
### General Linear Model Analysis with PC Scores as Response Variable

This app performs statistical tests on PC scores obtained from PCA (Principal Component Analysis)
to examine effects of **sex, cell subtype, batch, age**, and other factors.

---
""")

with st.expander("📚 Detailed Feature Description (Click to expand)", expanded=False):
    st.markdown("""
#### 📌 Main Use Cases:
- **scRNA-seq pseudobulk analysis**: Association between PC scores and phenotypes in sample-aggregated data
- **Unbalanced design support**: Handles unequal group sizes and missing data
- **Hierarchical structure**: Accounts for donor and technical batch hierarchies
- **Covariate adjustment**: Adjusts for age, sequencing depth, quality metrics, etc.

#### 🔬 Implemented Statistical Methods:

**1. OLS (Ordinary Least Squares)**
- Fixed effects model (sex, subtype, batch, age, etc.)
- **HC3 robust SE**: Heteroscedasticity-robust standard errors (recommended)
- **Cluster-robust SE**: Accounts for within-cluster correlation (e.g., by donor)
- Type II ANOVA: Significance test for each variable

**2. LMM (Linear Mixed Model)**
- **Random effects**: Supports hierarchical structures like `(1|donor)`, `(1|batch)`
- **REML estimation**: Unbiased variance estimation (recommended)
- **ML estimation**: For model comparison (AIC/BIC)
- **ICC calculation**: Intraclass correlation coefficient (random effect contribution)
- Warning displayed when number of donors is small (<5)

**3. Permutation Test**
- **Freedman-Lane method**: Rigorous test controlling for covariates (recommended)
- **Simple permutation**: Simple label shuffling
- Effective for small samples or when distributional assumptions are questionable
- Iterations: 1,000-50,000 (customizable)

**4. EMM (Estimated Marginal Means)**
- Estimation of group means adjusted for covariates
- Visualizes sex differences by subtype, etc.
- Easy interpretation with confidence intervals

#### ⚙️ Data Quality Checks:
- **Perfect collinearity detection**: Pre-warns about non-estimable effects
- **VIF calculation**: Multicollinearity diagnosis for numeric covariates
- **Sample size verification**: Sample counts per group and cross-tabulation
- **Missing value handling**: Automatic removal with reporting

#### 📊 Visualizations:
- **Forest plot**: Coefficient estimates with 95% confidence intervals
- **Diagnostic plots**: Residuals, Q-Q plot, Scale-Location
- **EMM plot**: Estimated marginal means by group
- **Permutation histogram**: Null distribution and observed statistic

#### 💾 Download Results:
- Coefficient table (CSV)
- ANOVA table (CSV)
- Variance components (LMM)
- Permutation distribution (CSV)
- All figures (saveable from Streamlit)
""")

# =============================================================================
# Practical Guide
# =============================================================================
with st.expander("📖 Practical Guide: Finding Sex-Associated Axis in Cell Type Pseudobulk Data", expanded=False):
    st.markdown("""
### 🎯 Use Case: Finding the "Sex Axis" with 4 Cell Types × Sex (8 Samples)

#### **Input Data Requirements**

**1 row = 1 sample (pseudobulk)**

Example columns:
```
sample_id    sex       cell.type    PC1      PC2      PC3      PC4      ...
sample_001   Female    art         -2.34     1.52    -0.87     0.45
sample_002   Male      art          1.87    -1.23     0.92    -0.35
sample_003   Female    cap1         0.23     2.11     1.34    -1.12
sample_004   Male      cap1        -1.45    -0.98     0.67     1.23
sample_005   Female    cap2         1.12     0.34    -1.56     0.89
sample_006   Male      cap2        -0.89     1.67     1.12    -0.67
sample_007   Female    vein         0.56    -1.34     0.23     1.78
sample_008   Male      vein        -1.23     0.45    -0.89     0.34
```

---

#### **App Settings (Column Mapping)**

1. **Main Variable 1 (main_var)** = `sex`
2. **Main Variable 2 / Blocking Variable (blocking_var)** = `cell.type`
3. **Analysis PC** = Select PC1 to PCk sequentially (or create a "sex composite axis" directly)

---

#### **Procedure A: Screening "Sex Main Effect" for Each PC (One Dimension at a Time)**

**Objective:** Determine which PC best explains sex while adjusting for cell.type as a covariate.

**Steps:**

1. **Select one PC** (e.g., PC3)

2. **Select OLS in the Analysis tab** and construct the model formula:
   ```
   PC3 ~ C(sex) + C(cell.type)
   ```
   - No interaction (main effects only)
   - **HC3 Robust** standard error is recommended (robust for small samples)

3. **Check the output**:
   - Verify the **sex coefficient** (Female vs Male) and its HC3 robust SE **p-value**
   - Record the coefficient **sign** and **effect size** (absolute value of coefficient)

4. **Turn ON Permutation Test** for more accurate p-values:
   - Select **Freedman-Lane method** (recommended)
   - Shuffles sex labels while maintaining strata (= cell.type)
   - Robust testing possible even with small samples (n=8)

5. **Check consistency with EMM (Estimated Marginal Means) plot**:
   - Verify if the **sign of Female vs Male difference is consistent** across cell.types
   - If all cell.types show the same direction of difference, it's evidence of consistent sex main effect

6. **Repeat for PC1 to PCk** and evaluate comprehensively:
   - ✅ **Effect size** (magnitude of sex coefficient)
   - ✅ **Consistency** (sign agreement across subtypes: check EMM plot)
   - ✅ **P-values** (both HC3 and Permutation)

   → Use these criteria to select the "**PC most associated with sex**"

---

#### **📝 Statistical Interpretation**

This OLS model is essentially equivalent to a **paired t-test (Δ = Female − Male mean)**:

- **Degrees of freedom**: df = number of cell.types − 1 = 3
- **Test target**: Sex main effect while controlling for cell.type as strata
- **Robust even for small samples**: Using HC3 Robust SE and Permutation together enables reliable inference even with n=8

**Notes:**
- Statistical power is highest when sample sizes are equal across cell.types (balanced design)
- When including interaction `C(sex) * C(cell.type)`, at least 3 samples per cell (sex × cell.type) are needed
- In this use case with only 1 sample per cell, interaction cannot be estimated (complete separation)

---

#### **🔍 Example Result Interpretation**

| PC  | sex coef | HC3 p-val | Perm p-val | Consistency | Overall |
|-----|----------|-----------|------------|-------------|---------|
| PC1 | 0.45     | 0.234     | 0.251      | ✅ All +    | Weak    |
| PC2 | -1.87    | 0.032     | 0.041      | ✅ All -    | **Strong** |
| PC3 | 0.23     | 0.678     | 0.712      | ± Mixed     | None    |
| PC4 | 1.12     | 0.098     | 0.112      | ✅ All +    | Moderate |

→ In this example, **PC2 is the sex axis** (separates in negative direction)

---

#### **💡 Next Steps (Applications)**

1. **Creating a "sex score" by combining multiple PCs**
   - Linear combination of PCs with large effect sizes (e.g., PC2, PC4)
   - Weighted average like `sex_score = -1.87 × PC2 + 1.12 × PC4`
   - Re-analyze using this composite axis as the new response variable

2. **Drilling down to individual gene level**
   - Extract genes with large loadings on the sex axis (e.g., PC2)
   - Interpret as gene groups contributing to sex differences

3. **Multiple testing correction**
   - When exploring multiple PCs (e.g., PC1-PC10), apply Benjamini-Hochberg FDR correction
   - When testing 10 PCs, adjust p < 0.05 → FDR < 0.05

---

#### **📚 Features Used in This Use Case**

- ✅ OLS with HC3 Robust SE
- ✅ Freedman-Lane Permutation Test
- ✅ EMM (Estimated Marginal Means) plot
- ✅ Diagnostic plots (residual normality and homoscedasticity)
- ✅ Coefficient table download (combine results from multiple PCs)

---
""")

# =============================================================================
# Practical Guide 2: Time Series × Genotype
# =============================================================================
with st.expander("📖 Practical Guide 2: WT vs KO × 3 Time Points Interaction Analysis", expanded=False):
    st.markdown("""
### 🎯 Use Case: Testing "Is Genotype Effect Time-Dependent" with WT/KO × 3 Time Points (n=3 each)

#### **Input Data Requirements**

**1 row = 1 sample (assuming independent samples)**

- **Design**: 2 genotypes (WT/KO) × 3 time points (t1/t2/t3) × n=3 per cell = **18 samples total**
- **Assumption**: Different individuals at each time point (independent samples)
- **For repeated measures on same individuals** → See "Repeated Measures" section below

Example columns:
```
sample_id    genotype    time    donor_id    PC1      PC2      PC3      PC4      ...
WT_t1_1      WT          t1      D01        -1.23     0.45    -0.67     1.12
WT_t1_2      WT          t1      D02         0.87    -1.34     0.92    -0.45
WT_t1_3      WT          t1      D03        -0.45     1.23    -1.01     0.78
WT_t2_1      WT          t2      D04         1.56    -0.23     1.34    -0.89
...(KO_t1, KO_t2, KO_t3 follow similarly)
```

---

#### **Statistical Model Strategy**

**Important**: First **test for interaction** → then examine simple effects if needed

1. **Interaction is significant** → genotype effect varies by time point (time-dependent)
   - Test WT vs KO at each time point separately (simple effects)
   - Test time-series changes within each genotype

2. **Interaction is not significant** → genotype effect is constant across time points
   - Drop interaction and interpret main effects only model
   - Evaluate genotype main effect and time main effect separately

---

#### **Approach A: Treat time as categorical (general solution, recommended)**

**Model formula:**
```
PC ~ C(genotype) * C(time)
```

**Characteristics:**
- Treats each time point as an independent category
- Handles non-linear time changes
- Most flexible and safe when only 3 time points

---

##### **Settings in the App (Categorical version)**

**1. Data Preparation**
- `genotype` column: "WT", "KO"
- `time` column: "t1", "t2", "t3" (as strings)
- `PC1`, `PC2`, ... columns

**2. Column Mapping**
- **Primary Variable 1**: `genotype`
- **Primary Variable 2 / Blocking Variable**: `time`
- **Analysis PC**: Select from PC1 to PCk
- **Donor/Subject ID**: (not required for independent samples)

**3. Model Settings**
- **Model Type**: OLS
- **Standard Error**: **HC3 Robust** (recommended, robust for small samples)
- **Interaction Term**: Add `genotype:time` <- **Important!**

**4. Test Options**
- **ANOVA Type**: Type II (recommended)
- **Permutation Test**: ON (Freedman-Lane method, recommended)
  - n=18 is a small sample, so supplement with permutation test

**5. Visualization**
- **EMM (Estimated Marginal Means)**: ON
  - Select both genotype and time
  - Confirm group differences in time-series trends with plots
- **Diagnostic Plots**: ON (check residual normality and homoscedasticity)

---

##### **Interpreting Results (Categorical version)**

**Step 1: Check interaction in ANOVA table**

| Term | Sum Sq | df | F value | p value |
|----|--------|----|---------| --------|
| C(genotype) | 15.3 | 1 | 8.23 | 0.012 |
| C(time) | 45.7 | 2 | 12.34 | 0.001 |
| **C(genotype):C(time)** | **23.4** | **2** | **6.32** | **0.011** <- Important |
| Residual | 22.3 | 12 | - | - |

**Interpretation:**
- **Interaction p=0.011 < 0.05** -> **Significant**
- -> **Genotype effect varies by time point** (time-dependent KO effect)

**Step 2: Check group differences at each time point in EMM table**

| genotype | time | mean | SE | 95% CI lower | 95% CI upper |
|----------|------|------|----|--------------| -------------|
| WT | t1 | -0.45 | 0.23 | -0.95 | 0.05 |
| KO | t1 | -0.52 | 0.23 | -1.02 | -0.02 |
| WT | t2 | 0.87 | 0.23 | 0.37 | 1.37 |
| KO | t2 | 1.89 | 0.23 | 1.39 | 2.39 |
| WT | t3 | 1.34 | 0.23 | 0.84 | 1.84 |
| KO | t3 | 2.78 | 0.23 | 2.28 | 3.28 |

**Interpretation:**
- **t1**: KO - WT = -0.07 (almost no difference)
- **t2**: KO - WT = 1.02 (KO increased)
- **t3**: KO - WT = 1.44 (difference further expanded)

-> **KO effect increases over time** (no difference at t1 -> large difference at t3)

**Step 3: Visualize with EMM plot**

With time on x-axis and PC score on y-axis:
- WT (blue line) rises gradually
- KO (red line) rises sharply
- Two lines not parallel = interaction present

---

##### **Simple Effects Test (when interaction is significant)**

If interaction is significant, **test WT vs KO at each time point separately**:

**Method 1: Subset Analysis (manual)**
1. Split data by time point (t1 only, t2 only, t3 only)
2. Test `PC ~ C(genotype)` for each subset
3. **Bonferroni correction**: Set p-value threshold to 0.05/3 = 0.0167
   or apply **BH-FDR correction**

**Method 2: EMM contrast (recommended, planned for future implementation)**
- Automatically calculates pairwise comparisons at each time point
- Automatically applies multiple testing correction

**Current implementation in the app:**
1. Split original file by time point in Excel (t1.tsv, t2.tsv, t3.tsv)
2. Upload each file separately and test with `PC ~ C(genotype)`
3. Manually apply BH correction to the 3 p-values

---

#### **Approach B: Treat time as numeric trend (linear)**

**Model formula:**
```
PC ~ C(genotype) * time_numeric
```

**Characteristics:**
- Treats time as a continuous variable (0, 1, 2)
- **Interaction coefficient = difference in slopes** (whether KO and WT have different time gradients)
- Can summarize time dependency with a single coefficient
- **Note**: Assumes linearity (risk of overfitting with quadratic for only 3 points)

---

##### **Settings in the App (Numeric version)**

**1. Data Preparation**
- Add `time_numeric` column: t1->0, t2->1, t3->2
- Or change `time` column directly to 0, 1, 2

**2. Column Mapping**
- **Primary Variable 1**: `genotype` (categorical)
- **Additional Covariate**: `time_numeric` (recognized as numeric)
- **Interaction Term**: Add `genotype:time_numeric`

**3. Model Settings**
- Same as Approach A for other settings

---

##### **Interpreting Results (Numeric version)**

**Coefficient Table:**

| Term | Coef | SE | t | p value |
|----|------|----|---|---------|
| Intercept | -0.45 | 0.18 | -2.50 | 0.027 |
| C(genotype)[T.KO] | -0.07 | 0.25 | -0.28 | 0.783 |
| time_numeric | 0.90 | 0.15 | 6.00 | <0.001 |
| **C(genotype)[T.KO]:time_numeric** | **0.75** | **0.21** | **3.57** | **0.003** <- Important |

**Interpretation:**
- **time_numeric coefficient = 0.90**: In WT group, PC score increases by 0.90 per unit time
- **Interaction coefficient = 0.75**: In KO group, PC score increases by **an additional 0.75** per unit time
  - KO group slope = 0.90 + 0.75 = **1.65**
  - WT group slope = 0.90
- **Interaction p=0.003** -> **Time trends differ significantly between KO and WT**

**Advantages:**
- Evaluates "difference in slopes" with a single test
- Effect size is intuitive (difference in slopes)

**Disadvantages:**
- Assumes linearity (difficult to verify with 3 points)
- May miss non-linear changes

---

#### **Approach C: Repeated measures (tracking same individuals)**

**Assumption:** Same donors (D01, D02, D03) measured at 3 time points

**Model formula (LMM):**
```
PC ~ C(genotype) * C(time) + (1|donor)
```

Or with random time slopes:
```
PC ~ C(genotype) * time_numeric + (1 + time_numeric|donor)
```

**Note:**
- With few donors (<5), random slope `(1+time|donor)` is unstable
- Start with `(1|donor)` first

---

##### **Settings in the App (LMM version)**

**1. Data Preparation**
- `donor` column is required
- Same donor ID appears in multiple rows (3 rows)

**2. Column Mapping**
- **Primary Variable 1**: `genotype`
- **Primary Variable 2**: `time`
- **Donor/Subject ID**: `donor` <- **Important!**
- **Interaction Term**: `genotype:time`

**3. Model Settings**
- **Model Type**: **LMM**
- **Estimation Method**: REML (recommended)
- Other settings same as OLS

---

##### **Interpreting Results (LMM version)**

**Fixed Effects Table:**

| Term | Coef | SE | z | p value |
|----|------|----|---|---------|
| Intercept | -0.45 | 0.28 | -1.61 | 0.108 |
| C(genotype)[T.KO] | -0.07 | 0.40 | -0.18 | 0.860 |
| C(time)[T.t2] | 1.32 | 0.20 | 6.60 | <0.001 |
| C(time)[T.t3] | 1.79 | 0.20 | 8.95 | <0.001 |
| C(genotype)[T.KO]:C(time)[T.t2] | 1.09 | 0.28 | 3.89 | <0.001 |
| C(genotype)[T.KO]:C(time)[T.t3] | 1.51 | 0.28 | 5.39 | <0.001 |

**Random Effects:**
- `Var(donor)` = 0.45 (between-donor variability)
- `Var(Residual)` = 0.23 (within-measurement variability)
- **ICC = 0.45 / (0.45 + 0.23) = 0.66**
  - 66% of variation is explained by individual differences between donors
  - Strong correlation in repeated measures -> LMM is appropriate

**Interpretation:**
- Interaction is significant -> genotype effect is time-dependent
- Conclusion after appropriately accounting for between-donor variation

---

#### **Reporting Guidelines (recommended)**

**1. Prioritize reporting interaction conclusions**
```
"The genotype x time interaction was significant (F(2,12)=6.32, p=0.011, HC3 robust SE;
Freedman-Lane permutation p=0.015). This indicates the KO effect is time-dependent."
```

**2. Create EMM figure**
- Time on x-axis, PC score on y-axis
- WT and KO means +/- 95%CI as line plots
- If possible, add difference plot **Delta(t) = KO - WT**

**3. Simple effects (when interaction is significant)**
```
"Post-hoc comparisons of WT vs KO at each time point:
- t1: Delta=-0.07 (95%CI: -0.68, 0.54), p=0.78 (not significant)
- t2: Delta=1.02 (95%CI: 0.41, 1.63), p=0.004 (significant)
- t3: Delta=1.44 (95%CI: 0.83, 2.05), p<0.001 (significant)
t2 and t3 remain significant after Bonferroni correction."
```

**4. Report effect sizes**
- Estimated differences and 95%CI at each time point
- Overall genotype effect (eta-squared for categorical, slope difference for numeric)

**5. Include permutation p-value**
- Consistency with HC3 provides evidence of robustness

---

#### **Cautions**

**1. Sample Size**
- n=3 per cell is minimal (low statistical power)
- n>=5 recommended if possible
- **Strongly recommend supplementing with permutation test**

**2. Multiple Testing Correction**
- Apply **Bonferroni** or **BH-FDR** for simple effects (3 comparisons)
- Same applies when testing multiple PCs

**3. Linear vs Categorical**
- **With only 3 points, linear assumption is risky** -> categorical version recommended
- Use linear version for exploratory analysis, confirm with categorical version

**4. Perfect Separation**
- With interaction term: 2x3=6 cells x n=3 = 18 samples
- df = 18 - 6 - 1 = 11 (limited degrees of freedom)
- Missing data or outliers can cause unstable estimates -> check with diagnostic plots

**5. For Repeated Measures**
- With few donors (<5), `(1+time|donor)` is unstable
- Start with `(1|donor)`, try adding slope if converges

---

#### **Features Used in This Use Case**

- OLS with HC3 Robust SE (independent samples)
- LMM with (1|donor) (repeated measures)
- Type II ANOVA (testing interaction)
- Freedman-Lane Permutation Test (supplementing small samples)
- EMM (Estimated Marginal Means) plot (visualizing time-series trends)
- Diagnostic plots (checking assumptions)
- Simple effects automatic calculation (planned for future, currently manual subset analysis)

---

#### **Practical Examples Summary**

| Situation | Recommended Model | Interaction Term | Test Method | EMM |
|------|-----------|----------|--------|-----|
| Independent samples, categorical | `PC ~ C(genotype) * C(time)` | Yes | OLS + HC3 + Perm | genotype x time |
| Independent samples, linear | `PC ~ C(genotype) * time_numeric` | Yes | OLS + HC3 + Perm | slope difference |
| Repeated measures | `PC ~ C(genotype) * C(time) + (1\|donor)` | Yes | LMM (REML) | genotype x time |

---

**In summary:**
1. **First test interaction `genotype x time`**
2. **If significant, proceed to simple effects at each time point** (currently manual subset analysis)
3. **If not significant, interpret with main effects model**
4. **Small sample (n=3/cell) requires permutation test for robustness**

With this workflow, you can statistically evaluate **time-dependent genotype effects** even with realistic sample sizes like WT/KO x 3 time points x n=3.

---
""")

st.markdown("---")

# =============================================================================
# Helper Functions
# =============================================================================

def remove_common_suffix(strings):
    """Remove common suffix from strings

    Same algorithm as DESeq2-LRT.py
    """
    if not strings or len(strings) == 0:
        return []
    # Get the length of the shortest string
    min_length = min(len(s) for s in strings)
    # Find the length of common suffix
    suffix_length = 0
    for i in range(1, min_length + 1):
        suffix = strings[0][-i:]
        if all(s.endswith(suffix) for s in strings):
            suffix_length = i
        else:
            break
    # If no common suffix found, return original list
    if suffix_length == 0:
        return strings
    # Create new list with common suffix removed
    return [s[:-suffix_length] for s in strings]

def detect_collinearity(df, formula_terms):
    """Detect perfect collinearity

    Detects cases where categorical variables are perfectly confounded
    (one determines the other). In such cases, the model cannot be estimated.

    Example: If a subtype is entirely Female, sex effect and subtype effect cannot be separated.
    """
    issues = []

    # Check only categorical variables
    cat_vars = [col for col in formula_terms if col in df.columns and df[col].dtype == 'object']

    for i, var1 in enumerate(cat_vars):
        for var2 in cat_vars[i+1:]:
            cross_tab = pd.crosstab(df[var1], df[var2])
            # Perfect confounding if each row or column has only one non-zero element
            if (cross_tab > 0).sum(axis=0).min() == 1 or (cross_tab > 0).sum(axis=1).min() == 1:
                issues.append(f"Warning: **Perfect confounding**: `{var1}` and `{var2}` are perfectly confounded")

    return issues

def calculate_vif(df, numeric_cols):
    """Calculate VIF (Variance Inflation Factor)

    VIF > 10: Strong multicollinearity present
    VIF > 5: Moderate multicollinearity present
    VIF < 5: No problem
    """
    from statsmodels.stats.outliers_influence import variance_inflation_factor

    vif_data = pd.DataFrame()
    vif_data["Variable"] = numeric_cols
    vif_data["VIF"] = [variance_inflation_factor(df[numeric_cols].values, i)
                       for i in range(len(numeric_cols))]
    return vif_data

def freedman_lane_permutation(y, X_full_df, interest_var_name, n_perm=10000, random_state=42,
                               stratify_by=None, one_sided=None, design_info=None):
    """Freedman-Lane permutation test (robust version)

    Tests the effect of a specific variable (e.g., sex) while controlling for covariates.

    Procedure:
    1. Fit full model (all variables) and null model (excluding test variable)
    2. Permute residuals from null model (with stratification option)
    3. Construct pseudo-response variable from permuted residuals
    4. Refit full model and calculate test statistic
    5. Iterate to create null distribution

    Parameters:
    -----------
    y : array-like
        Response variable (PC scores)
    X_full_df : DataFrame
        Full design matrix (with column names)
    interest_var_name : str
        Name of the variable to test (e.g., 'sex'). All columns containing this variable
        are excluded to create the null model
    n_perm : int
        Number of permutations
    random_state : int
        Random seed
    stratify_by : array-like, optional
        Stratification variable (e.g., subtype). If specified, residuals are permuted
        only within each stratum
    one_sided : str, optional
        Direction for one-sided test. 'greater' (obs > null) or 'less' (obs < null)
        None for two-sided test
    design_info : patsy.DesignInfo, optional
        Design matrix information (for robust column identification)

    Returns:
    --------
    dict: observed (observed statistic), null_distribution, p_value,
          method_info (information for reproducibility)
    """
    import streamlit as st
    np.random.seed(random_state)

    # Fit full model
    model_full = sm.OLS(y, X_full_df).fit()

    # Identify columns containing the test variable (robust with design_info)
    if design_info is not None:
        # Robust column identification using patsy design_info
        interest_cols = []
        for term in design_info.terms:
            term_name = term.name()
            # Identify terms containing the test variable (main effect + interactions)
            if term_name != 'Intercept' and interest_var_name in term_name:
                # Get columns corresponding to this term
                slice_obj = design_info.term_name_slices[term_name]
                cols = X_full_df.columns[slice_obj]
                interest_cols.extend(cols)
    else:
        # Fallback: string matching (for backward compatibility)
        interest_cols = [col for col in X_full_df.columns
                         if f'C({interest_var_name})' in col and col != 'Intercept']

    if not interest_cols:
        raise ValueError(f"No columns found for test variable '{interest_var_name}'")

    # Observed statistic (using absolute t-value of main effect column - robust)
    main_effect_col = [col for col in interest_cols if ':' not in col]
    if not main_effect_col:
        raise ValueError(f"No main effect column found: {interest_cols}")

    # Use absolute t-statistic (to avoid sign dependence)
    obs_t = float(model_full.tvalues[main_effect_col[0]])
    obs_stat = np.abs(obs_t)

    # Null model (excluding all columns of test variable)
    X_reduced = X_full_df.drop(columns=interest_cols)

    if X_reduced.shape[1] == 1:  # Intercept only
        residuals = y - np.mean(y)
        fitted_reduced = np.full_like(y, np.mean(y))
    else:
        model_reduced = sm.OLS(y, X_reduced).fit()
        residuals = np.array(model_reduced.resid)  # Convert to numpy array
        fitted_reduced = np.array(model_reduced.predict())  # Convert to numpy array


    # Save null and alternative model formulas (for reporting)
    reduced_formula = "y ~ " + " + ".join([col for col in X_reduced.columns if col != 'Intercept'])
    if reduced_formula == "y ~ ":
        reduced_formula = "y ~ 1"  # Intercept only
    full_formula = "y ~ " + " + ".join([col for col in X_full_df.columns if col != 'Intercept'])

    # Permutation
    null_dist = []
    null_dist_abs = []  # Also store absolute values
    for _ in range(n_perm):
        # Permute residuals (with stratification option)
        if stratify_by is not None:
            # Stratified permutation: shuffle residuals only within each stratum
            perm_residuals = np.zeros_like(residuals)
            for stratum in np.unique(stratify_by):
                stratum_mask = (stratify_by == stratum)
                stratum_indices = np.where(stratum_mask)[0]
                perm_indices = np.random.permutation(stratum_indices)
                perm_residuals[stratum_indices] = residuals[perm_indices]
        else:
            # Non-stratified permutation: shuffle all
            perm_idx = np.random.permutation(len(residuals))
            perm_residuals = residuals[perm_idx]

        y_perm = fitted_reduced + perm_residuals

        # Refit full model
        model_perm = sm.OLS(y_perm, X_full_df).fit()
        perm_t = float(model_perm.tvalues[main_effect_col[0]])
        null_dist.append(perm_t)
        null_dist_abs.append(np.abs(perm_t))

    null_dist = np.array(null_dist)
    null_dist_abs = np.array(null_dist_abs)

    # Calculate p-value (one-sided/two-sided) + continuity correction
    if one_sided == 'greater':
        # One-sided (greater direction): use obs_t
        p_value = (np.sum(null_dist >= obs_t) + 1) / (n_perm + 1)
    elif one_sided == 'less':
        # One-sided (less direction): use obs_t
        p_value = (np.sum(null_dist <= obs_t) + 1) / (n_perm + 1)
    else:
        # Two-sided test: compare absolute values
        p_value = (np.sum(null_dist_abs >= obs_stat) + 1) / (n_perm + 1)

    return {
        'observed': obs_t,  # Return original t-value (with sign)
        'observed_abs': obs_stat,  # Also return absolute value
        'null_distribution': null_dist,  # Distribution of signed t-values
        'null_distribution_abs': null_dist_abs,  # Distribution of absolute values
        'p_value': p_value,
        'stratified': stratify_by is not None,
        'one_sided': one_sided,
        'interest_cols': interest_cols,
        'main_effect_col': main_effect_col[0],
        'method_info': {  # Information for reproducibility
            'n_permutations': n_perm,
            'random_seed': random_state,
            'test_statistic': '|t|' if one_sided is None else 't',
            'continuity_correction': True,
            'reduced_model': reduced_formula,
            'full_model': full_formula,
            'stratification_var': 'None' if stratify_by is None else 'Provided'
        }
    }

def calculate_emm(model, df, factors, numeric_vars_at_mean=True):
    """Calculate Estimated Marginal Means (EMM)

    Calculates predicted mean values for each level of main factors
    while fixing covariates at specific values (usually the mean).

    Example: Calculate predicted PC score means for each sex x subtype combination
             while fixing age at its mean value

    Parameters:
    -----------
    model : fitted statsmodels model
    df : DataFrame
        Original data
    factors : list
        List of factors to calculate EMM for
    numeric_vars_at_mean : bool
        Whether to fix numeric variables at their mean

    Returns:
    --------
    DataFrame: Predicted mean, SE, 95% CI for each combination
    """
    from itertools import product

    # Create reference grid
    ref_grid = {}

    for col in df.columns:
        if col in factors:
            ref_grid[col] = df[col].unique()
        elif pd.api.types.is_numeric_dtype(df[col]):
            if numeric_vars_at_mean:
                ref_grid[col] = [df[col].mean()]
            else:
                ref_grid[col] = df[col].unique()

    # Generate all combinations
    grid_combos = list(product(*[ref_grid[col] for col in df.columns if col in ref_grid]))

    # Create prediction dataframe
    pred_df = pd.DataFrame(grid_combos, columns=[col for col in df.columns if col in ref_grid])

    # Predict
    predictions = model.get_prediction(pred_df)
    pred_summary = predictions.summary_frame(alpha=0.05)

    result_df = pred_df.copy()
    result_df['mean'] = pred_summary['mean']
    result_df['se'] = pred_summary['mean_se']
    result_df['ci_lower'] = pred_summary['mean_ci_lower']
    result_df['ci_upper'] = pred_summary['mean_ci_upper']

    return result_df

def plot_forest(coef_df, title="Forest Plot"):
    """Create Forest plot for coefficients"""
    fig, ax = plt.subplots(figsize=(8, max(6, len(coef_df) * 0.5)))

    y_pos = np.arange(len(coef_df))
    ax.errorbar(coef_df['Coef'], y_pos,
                xerr=[coef_df['Coef'] - coef_df['CI_lower'],
                      coef_df['CI_upper'] - coef_df['Coef']],
                fmt='o', markersize=8, capsize=5, capthick=2, color='steelblue')

    ax.axvline(0, color='red', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(coef_df.index)
    ax.set_xlabel('Coefficient Estimate (95% CI)', fontsize=11)
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    plt.tight_layout()

    return fig

def plot_diagnostic(model, title="Diagnostic Plots"):
    """Model diagnostic plots

    Creates 4 diagnostic plots:
    1. Residuals vs Fitted: Check for non-linearity, homoscedasticity
    2. Q-Q plot: Check normality of residuals
    3. Scale-Location: Check homoscedasticity (standardized residuals)
    4. Residual histogram: Distribution of residuals
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    fitted = model.fittedvalues
    resid = model.resid

    # Residuals vs Fitted
    axes[0, 0].scatter(fitted, resid, alpha=0.6, s=50, edgecolors='k', linewidths=0.5)
    axes[0, 0].axhline(0, color='red', linestyle='--', linewidth=1.5)
    axes[0, 0].set_xlabel('Fitted values', fontsize=11)
    axes[0, 0].set_ylabel('Residuals', fontsize=11)
    axes[0, 0].set_title('Residuals vs Fitted', fontsize=12, fontweight='bold')
    axes[0, 0].grid(alpha=0.3)

    # Q-Q plot
    stats.probplot(resid, dist="norm", plot=axes[0, 1])
    axes[0, 1].set_title('Normal Q-Q Plot', fontsize=12, fontweight='bold')
    axes[0, 1].grid(alpha=0.3)

    # Scale-Location
    standardized_resid = resid / np.std(resid)
    axes[1, 0].scatter(fitted, np.sqrt(np.abs(standardized_resid)),
                      alpha=0.6, s=50, edgecolors='k', linewidths=0.5)
    axes[1, 0].set_xlabel('Fitted values', fontsize=11)
    axes[1, 0].set_ylabel('√|Standardized residuals|', fontsize=11)
    axes[1, 0].set_title('Scale-Location', fontsize=12, fontweight='bold')
    axes[1, 0].grid(alpha=0.3)

    # Residuals histogram
    axes[1, 1].hist(resid, bins=min(30, len(resid)//3),
                   edgecolor='black', alpha=0.7, color='steelblue')
    axes[1, 1].set_xlabel('Residuals', fontsize=11)
    axes[1, 1].set_ylabel('Frequency', fontsize=11)
    axes[1, 1].set_title('Residual Distribution', fontsize=12, fontweight='bold')
    axes[1, 1].grid(alpha=0.3, axis='y')

    plt.suptitle(title, fontsize=14, fontweight='bold', y=1.00)
    plt.tight_layout()
    return fig

def create_results_zip(pc_col):
    """
    Bundle analysis results and plots into a zip file

    Parameters:
    -----------
    pc_col : str
        PC column name (used for filename)

    Returns:
    --------
    bytes : Binary data of zip file
    """
    zip_buffer = io.BytesIO()

    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        base_name = f"pca_stats_{pc_col}_{timestamp}"

        # 0. Welch t-test results
        if 'welch_result' in st.session_state:
            welch_res = st.session_state.welch_result
            welch_df = pd.DataFrame({
                'test': ['Welch\'s t-test'],
                'group1': [welch_res['group1']],
                'group2': [welch_res['group2']],
                'group1_mean': [welch_res['group1_mean']],
                'group2_mean': [welch_res['group2_mean']],
                't_value': [welch_res['t']],
                'p_value': [welch_res['p']],
                'cohens_d': [welch_res['cohens_d']],
                'hedges_g': [welch_res['hedges_g']]
            })
            zip_file.writestr(f"{base_name}/tables/welch_ttest.csv", welch_df.to_csv(index=False))

        # 1. OLS coefficient table
        if 'ols_coef' in st.session_state:
            csv_data = st.session_state.ols_coef.to_csv()
            zip_file.writestr(f"{base_name}/tables/ols_coefficients.csv", csv_data)

        # 2. LMM coefficient table
        if 'lmm_coef' in st.session_state:
            csv_data = st.session_state.lmm_coef.to_csv()
            zip_file.writestr(f"{base_name}/tables/lmm_coefficients.csv", csv_data)

        # 3. Permutation test results
        if 'perm_result' in st.session_state:
            perm_result = st.session_state.perm_result

            # Statistics summary
            perm_summary = pd.DataFrame({
                'observed_stat': [perm_result['observed']],
                'p_value': [perm_result['p_value']],
                'n_permutations': [len(perm_result['null_distribution'])],
                'method': [perm_result.get('method', 'unknown')],
                'stratified': [perm_result.get('stratified', False)]
            })
            zip_file.writestr(f"{base_name}/tables/permutation_summary.csv",
                            perm_summary.to_csv(index=False))

            # Null distribution
            perm_df = pd.DataFrame({
                'permutation_stat': perm_result['null_distribution']
            })
            zip_file.writestr(f"{base_name}/tables/permutation_null_distribution.csv",
                            perm_df.to_csv(index=False))

        # 4. ANOVA table (OLS)
        if 'ols_anova' in st.session_state:
            csv_data = st.session_state.ols_anova.to_csv()
            zip_file.writestr(f"{base_name}/tables/ols_anova.csv", csv_data)

        # 5. EMM results
        if 'ols_emm' in st.session_state:
            csv_data = st.session_state.ols_emm.to_csv(index=False)
            zip_file.writestr(f"{base_name}/tables/estimated_marginal_means.csv", csv_data)

        # 6. Save plots
        # Forest plot (OLS)
        if 'fig_forest_ols' in st.session_state:
            img_buffer = io.BytesIO()
            st.session_state.fig_forest_ols.savefig(img_buffer, format='png', dpi=300, bbox_inches='tight')
            zip_file.writestr(f"{base_name}/figures/ols_forest_plot.png", img_buffer.getvalue())

            # Also save PDF version
            pdf_buffer = io.BytesIO()
            st.session_state.fig_forest_ols.savefig(pdf_buffer, format='pdf', bbox_inches='tight')
            zip_file.writestr(f"{base_name}/figures/ols_forest_plot.pdf", pdf_buffer.getvalue())

        # Diagnostic plot
        if 'fig_diagnostic' in st.session_state:
            img_buffer = io.BytesIO()
            st.session_state.fig_diagnostic.savefig(img_buffer, format='png', dpi=300, bbox_inches='tight')
            zip_file.writestr(f"{base_name}/figures/ols_diagnostic_plots.png", img_buffer.getvalue())

            pdf_buffer = io.BytesIO()
            st.session_state.fig_diagnostic.savefig(pdf_buffer, format='pdf', bbox_inches='tight')
            zip_file.writestr(f"{base_name}/figures/ols_diagnostic_plots.pdf", pdf_buffer.getvalue())

        # EMM plot
        if 'fig_emm' in st.session_state:
            img_buffer = io.BytesIO()
            st.session_state.fig_emm.savefig(img_buffer, format='png', dpi=300, bbox_inches='tight')
            zip_file.writestr(f"{base_name}/figures/estimated_marginal_means.png", img_buffer.getvalue())

            pdf_buffer = io.BytesIO()
            st.session_state.fig_emm.savefig(pdf_buffer, format='pdf', bbox_inches='tight')
            zip_file.writestr(f"{base_name}/figures/estimated_marginal_means.pdf", pdf_buffer.getvalue())

        # Forest plot (LMM)
        if 'fig_forest_lmm' in st.session_state:
            img_buffer = io.BytesIO()
            st.session_state.fig_forest_lmm.savefig(img_buffer, format='png', dpi=300, bbox_inches='tight')
            zip_file.writestr(f"{base_name}/figures/lmm_forest_plot.png", img_buffer.getvalue())

            pdf_buffer = io.BytesIO()
            st.session_state.fig_forest_lmm.savefig(pdf_buffer, format='pdf', bbox_inches='tight')
            zip_file.writestr(f"{base_name}/figures/lmm_forest_plot.pdf", pdf_buffer.getvalue())

        # Null distribution plot (Permutation test)
        if 'fig_null_dist' in st.session_state:
            img_buffer = io.BytesIO()
            st.session_state.fig_null_dist.savefig(img_buffer, format='png', dpi=300, bbox_inches='tight')
            zip_file.writestr(f"{base_name}/figures/permutation_null_distribution.png", img_buffer.getvalue())

            pdf_buffer = io.BytesIO()
            st.session_state.fig_null_dist.savefig(pdf_buffer, format='pdf', bbox_inches='tight')
            zip_file.writestr(f"{base_name}/figures/permutation_null_distribution.pdf", pdf_buffer.getvalue())

        # Create README
        readme_content = f"""PCA Statistics Analysis Results
================================

Analysis Date: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
PC Column: {pc_col}

Directory Structure:
-------------------
tables/          - CSV files with statistical results
  - ols_coefficients.csv           : OLS regression coefficients
  - lmm_coefficients.csv           : LMM fixed effects coefficients
  - ols_anova.csv                  : ANOVA table
  - estimated_marginal_means.csv   : EMM estimates
  - permutation_summary.csv        : Permutation test summary
  - permutation_null_distribution.csv : Null distribution values

figures/         - Publication-quality figures (PNG 300dpi + PDF)
  - ols_forest_plot.*              : OLS coefficient forest plot
  - lmm_forest_plot.*              : LMM coefficient forest plot
  - ols_diagnostic_plots.*         : Model diagnostic plots
  - estimated_marginal_means.*     : EMM visualization
  - permutation_null_distribution.*: Permutation test null distribution

Generated by: PCA Statistics Streamlit App
"""
        zip_file.writestr(f"{base_name}/README.txt", readme_content)

    zip_buffer.seek(0)
    return zip_buffer.getvalue()

# =============================================================================
# Sidebar Options
# =============================================================================

with st.sidebar:
    st.title("⚙️ Analysis Options")

    st.markdown("### Select Statistical Method")
    analysis_method = st.radio(
        "Analysis method:",
        ["OLS (Ordinary Least Squares)",
         "LMM (Linear Mixed Model)",
         "Both (OLS + LMM)"],
        index=0,
        help="OLS: Fixed effects only. LMM: Hierarchical model with random effects (donor, batch, etc.)"
    )

    if "OLS" in analysis_method:
        st.markdown("#### OLS Options")
        se_type = st.radio(
            "Standard error type:",
            ["Classical", "HC3 (Robust, recommended)", "Cluster-robust"],
            index=1,
            help="HC3: Robust to heteroscedasticity. Cluster-robust: Accounts for within-cluster correlation"
        )

        if se_type == "Cluster-robust":
            st.info("📌 Cluster variable will be selected after data upload")

        # ANOVA Type selection
        anova_type = st.radio(
            "ANOVA Type:",
            ["Type II (recommended)", "Type III"],
            index=0,
            help="Type II: Tests main effects adjusted for other variables (recommended when no interactions)\nType III: Tests including all interactions (recommended when interactions are present)"
        )

        # WLS option
        use_wls = st.checkbox(
            "Use WLS (Weighted Least Squares)",
            value=False,
            help="When variance is heterogeneous, estimates with weights for each observation. Weights are automatically calculated from inverse residual variance"
        )

        if use_wls:
            with st.expander("⚠️ Important Notes on WLS Usage (Required Reading)"):
                st.markdown("""
                ### Appropriate Use of WLS (Weighted Least Squares)

                **⚠️ Important**: WLS can produce unstable estimates depending on how weights are determined, and is **not recommended for small samples**

                ---

                ### 🚫 When to Avoid WLS (Important)

                #### **Generally NOT recommended for small sample designs**

                Example: subtype=4 groups, each group n=1-2 (very small design)
                - ❌ **Difficult to create justified weights** (variance estimation is unstable)
                - ❌ **Estimates are prone to fluctuation** (high uncertainty in weights)
                - ❌ **Risk of overfitting** (low degrees of freedom)

                **Methods to use instead**:
                - ✅ **Welch-type + Type II ANOVA** (primary analysis)
                - ✅ **OLS + HC3 robust SE** (coefficients & 95% CI)
                - ✅ **Freedman-Lane permutation test** (conditional null distribution)

                These methods are robust to heteroscedasticity **without assuming weights**, and are less prone to fluctuation even with small samples.

                ---

                ### ✅ When WLS is Effective (Limited Cases)

                **Only when there are sufficient replicates and between-group variance differences can be clearly estimated**

                #### Required conditions (ALL must be met):

                1. **Sufficient sample size per group**
                   - Minimum n >= 10-15 per group
                   - Variance can be stably estimated

                2. **Clear variance differences between groups**
                   - Detectable by Levene's or Bartlett's test
                   - Visually apparent in diagnostic plots

                3. **Clear basis for weights**
                   - Within-group inverse variance known from prior information
                   - Based on external information such as number of technical replicates

                4. **Model assumptions are satisfied**
                   - Linearity, independence
                   - Weighted residuals are normally distributed

                #### Recommendations in this case:

                - **Primary analysis candidate**: WLS (inverse variance weights, weights explicitly stated)
                - **HC3 usually not needed** (if WLS achieves homoscedasticity)
                - Include as **"reference"** if needed (sensitivity analysis)

                ---

                ### 📊 Guidelines for Method Selection

                #### If your design has small samples (e.g., 4 groups x n=1-2 each):

                ```
                [Primary Analysis]
                ✅ Welch-type + Type II ANOVA
                   PC_k ~ C(sex) + C(subtype)
                   -> Robust to heteroscedasticity, stable even with small samples

                [Also Report]
                ✅ OLS + HC3 robust SE
                   -> Coefficient estimates & 95% CI
                   -> Robust p-values

                ✅ Freedman-Lane permutation test
                   -> Conditional null distribution
                   -> No distributional assumptions

                [Sensitivity Analysis Only]
                ⚠️ WLS
                   -> Weights = within-group inverse variance, etc. (explicitly stated)
                   -> Do not use as primary analysis
                ```

                #### When sufficient replicates are available (n >= 15 per group):

                ```
                [Primary Analysis Candidates]
                ✅ WLS (inverse variance weights)
                   -> Efficient estimation
                   -> Weights explicitly stated

                [Not Needed/Reference]
                - HC3 usually not needed
                  (homoscedasticity achieved with WLS)
                - Can include as sensitivity analysis if needed
                ```

                ---

                ### ⚠️ Risks of WLS

                #### When using WLS with small samples:

                1. **Large estimation error in weights**
                   - Cannot accurately estimate true variance
                   - Estimates are distorted by incorrect weights

                2. **Overfitting**
                   - Fitting to observations is too strong
                   - Generalization performance decreases

                3. **Extreme weights**
                   - Extreme weights assigned to outliers
                   - Numerical instability

                4. **Difficulty in interpretation**
                   - Meaning of weighted average effect is unclear
                   - Low reproducibility

                ---

                ### 🔬 Implementation in This App

                **Method**: Residual-based inverse variance weights

                ```python
                1. Calculate residuals with standard OLS
                2. Weights = 1 / (residuals^2)
                3. Normalize weights (mean=1)
                4. Re-estimate with WLS
                ```

                **⚠️ Problems with this method**:
                - Based on residuals, so differs from true variance
                - Residuals are unstable with small samples
                - **Recommended**: Determine weights from prior information (technical replicates, known variance)

                ---

                ### 💡 Recommendations (Summary)

                #### Small samples (n < 15/group):
                1. ❌ **Do not use WLS**
                2. ✅ **Use Welch-type + HC3 as primary analysis**
                3. ✅ **Supplement with permutation test**

                #### Large samples (n >= 15/group):
                1. ✅ **Visualize and test variance differences**
                2. ✅ **Specify weights explicitly for WLS**
                3. ⚠️ **Compare with HC3 (sensitivity analysis)**

                #### Always:
                - 📊 **Diagnostic Plots required** (check residual patterns)
                - 📝 **Document the rationale for weights** (for publication)
                - 🔄 **Verify robustness of results using multiple methods**

                ---

                ### 📚 Reference: Comparison of Methods for Handling Heteroscedasticity

                | Method | Sample Size Requirement | Robustness | Recommendation (small n) | Recommendation (large n) |
                |--------|------------------------|------------|--------------------------|--------------------------|
                | **Welch t-test** | Small OK | High | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ |
                | **HC3 Robust SE** | Small OK | High | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ |
                | **Permutation Test** | Small OK | Very High | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
                | **WLS** | Large Required | Low (small n) | ⭐ | ⭐⭐⭐⭐⭐ |
                | **Standard OLS** | - | None | ❌ | ⭐ |

                **Conclusion**: If your design has small samples, avoid WLS and use Welch + HC3 + Permutation Test instead.
                """)


    if "LMM" in analysis_method:
        st.markdown("#### LMM Options")
        reml_method = st.radio(
            "Estimation Method:",
            ["REML (Recommended)", "ML"],
            index=0,
            help="REML: Unbiased variance estimation. ML: For model comparison (AIC/BIC)"
        )

    st.markdown("---")
    st.markdown("### Welch's t-test")
    run_welch = st.checkbox("Run Welch's t-test", value=False,
                           help="**Only valid when main variable has 2 groups**\n\n"
                                "Two-group comparison without assuming equal variance (no covariate adjustment)\n"
                                "Calculates reference values as a simple two-group comparison")

    if run_welch:
        st.info("**About Welch's t-test:**\n\n"
                "- Direct comparison between two groups of the main variable (no covariate adjustment)\n"
                "- Robust test that does not assume equal variance\n"
                "- Will not be performed if main variable has 3 or more groups")

    st.markdown("---")
    st.markdown("### Permutation Test")
    run_permutation = st.checkbox("Run Permutation Test", value=False,
                                  help="**Test target: Main Variable 1 only**\n\n"
                                       "Recommended for small samples or when distributional assumptions are questionable\n"
                                       "Main Variable 2 and other covariates are treated as control variables")

    if run_permutation:
        st.info("**Target variable for Permutation Test:** Only Main Variable 1 is tested.\n\n"
                "Main Variable 2 and other variables are controlled as covariates.")

        n_permutations = st.slider("Number of Permutations:",
                                    min_value=1000, max_value=50000,
                                    value=10000, step=1000,
                                    help="More permutations = more accurate, but takes longer")
        perm_method = st.radio(
            "Permutation Strategy:",
            ["Freedman-Lane Method (Recommended)", "Simple label permutation"],
            index=0,
            help="Freedman-Lane: Rigorous test controlling for covariates\n"
                 "Tests the effect of Main Variable 1 while controlling for other variables"
        )

        # One-sided test option
        perm_sided = st.radio(
            "Test Direction:",
            ["Two-sided (two-sided)", "One-sided: Main Variable 1 > Reference (greater)", "One-sided: Main Variable 1 < Reference (less)"],
            index=0,
            help="**Two-sided (Recommended)**: Tests for any difference regardless of direction\n\n"
                 "**One-sided (greater)**: Use only when you expect 'Main Variable 1 > Reference' a priori\n"
                 "Example: Hypothesis that Female > Male\n\n"
                 "**One-sided (less)**: Use only when you expect 'Main Variable 1 < Reference' a priori\n\n"
                 "Use one-sided tests only when you have a prior hypothesis"
        )

    st.markdown("---")
    st.markdown("### Visualization Options")
    show_diagnostics = st.checkbox("Show Diagnostic Plots", value=True,
                                   help="Residual plots, Q-Q plots, etc.")
    show_emm = st.checkbox("Calculate Estimated Marginal Means (EMM)", value=True,
                          help="Calculate adjusted means by group")

# =============================================================================
# File Upload
# =============================================================================

st.markdown("## 1. Data Upload")

# Select data input method
input_method = st.radio(
    "Select Data Input Method:",
    ("PCA scores file + Enter metadata manually", "Upload PCA scores file (with metadata)"),
    help="When using a PCA scores file downloaded from pca.py, you can enter metadata manually."
)

if input_method == "PCA scores file + Enter metadata manually":
    # Enter metadata manually
    st.markdown("""
    **Enter PCA scores file and metadata separately**

    1. Upload PCA scores file downloaded from pca.py (containing sample names and PC1, PC2, ...)
    2. Enter metadata in the form below
    """)

    uploaded_file = st.file_uploader("Select PCA scores file", type=["tsv", "csv", "txt"],
                                     key="pca_scores_file")

else:  # Upload file with metadata
    st.markdown("""
    **Required data format**: TSV or CSV file

    - **1 row = 1 sample** (pseudobulk level)
    - **PC score columns**: PC1, PC2, PC3, ... etc.
    - **Categorical variables**: sex, subtype, condition, etc.
    - **Optional: Donor/Subject ID** (if repeated measures exist)
    - **Optional: Numeric covariates**: age, batch, depth, etc.
    """)

    uploaded_file = st.file_uploader("Select file", type=["tsv", "csv", "txt"])

if uploaded_file is not None:
    # Read file
    try:
        df = pd.read_csv(uploaded_file, sep="\t", index_col=0)
    except:
        uploaded_file.seek(0)
        df = pd.read_csv(uploaded_file, sep=",", index_col=0)

    # Reset index to convert to regular column
    df = df.reset_index()

    # Set appropriate index column name
    if df.columns[0] == 'index':
        df.rename(columns={'index': 'sample_id'}, inplace=True)

    st.success(f"File loaded successfully: {df.shape[0]} rows x {df.shape[1]} columns")

    with st.expander("Data Preview (first 10 rows)", expanded=True):
        st.dataframe(df.head(10))

    # Initialize session_state
    if 'metadata_submitted' not in st.session_state:
        st.session_state.metadata_submitted = False
    if 'metadata_df' not in st.session_state:
        st.session_state.metadata_df = None
    if 'combined_df' not in st.session_state:
        st.session_state.combined_df = None

    # Manual metadata entry
    if input_method == "PCA scores file + Enter metadata manually":
        st.markdown("## 2. Metadata Entry")
        st.markdown("Enter metadata for each sample")

        # Get sample names (first column or row names)
        if df.index.name or (df.index[0] != 0):  # If index is sample names
            sample_names = df.index.tolist()
        else:  # If first column is sample names
            sample_names = df.iloc[:, 0].tolist()

        # Metadata column name entry
        st.markdown("### Metadata Column Settings")
        col1, col2 = st.columns(2)

        with col1:
            meta_cols_input = st.text_area(
                "Metadata column names (comma or newline separated):",
                value="sex, subtype",
                help="Example: sex, subtype, donor, age, batch"
            )

        # Parse metadata column names
        if ',' in meta_cols_input:
            meta_col_names = [x.strip() for x in meta_cols_input.split(',') if x.strip()]
        else:
            meta_col_names = [x.strip() for x in meta_cols_input.split('\n') if x.strip()]

        with col2:
            st.markdown("**Columns to enter:**")
            for col in meta_col_names:
                st.write(f"- {col}")

        # Metadata entry form
        st.markdown("### Metadata Entry")

        with st.form("metadata_input"):
            # Enter metadata for each sample using data editor
            meta_df = pd.DataFrame(index=sample_names, columns=meta_col_names)

            # Set default values using same algorithm as DESeq2-LRT
            # 1. Remove common suffix
            sample_names_str = [str(s) for s in sample_names]
            group_base = remove_common_suffix(sample_names_str)
            # 2. Remove trailing numbers to infer group names
            group_names = [remove_sample_num(s) for s in group_base]

            # Set default values for each metadata column
            for col in meta_col_names:
                col_lower = col.lower()

                if col_lower in ['sex', 'gender']:
                    # Infer sex from group name (if male/female is included)
                    def infer_sex(name):
                        name_lower = name.lower()
                        if 'male' in name_lower and 'female' not in name_lower:
                            return 'Male'
                        elif 'female' in name_lower or 'fem' in name_lower:
                            return 'Female'
                        else:
                            return 'Male'  # Default
                    meta_df[col] = [infer_sex(g) for g in group_names]

                elif col_lower in ['donor', 'subject', 'patient', 'individual']:
                    # Infer donor ID from sample name
                    # Pattern 1: sample_D01_A -> D01
                    # Pattern 2: D01_typeA -> D01
                    # Pattern 3: Assign sequential numbers
                    donor_ids = []
                    for s in sample_names_str:
                        parts = s.split('_')
                        # Look for patterns like D01
                        donor_found = False
                        for part in parts:
                            if part.startswith('D') and len(part) >= 2 and part[1:].isdigit():
                                donor_ids.append(part)
                                donor_found = True
                                break
                        if not donor_found:
                            # Use group name
                            donor_ids.append(group_names[len(donor_ids)])
                    meta_df[col] = donor_ids

                elif col_lower in ['subtype', 'celltype', 'cell_type', 'type', 'tissue']:
                    # Use group name as is (often represents subtype)
                    meta_df[col] = group_names

                elif col_lower in ['batch', 'library', 'plate', 'run']:
                    # Infer batch number
                    # Pattern 1: sample_B1_typeA -> B1
                    # Pattern 2: If same for all samples, use Batch1
                    batch_ids = []
                    for s in sample_names_str:
                        parts = s.split('_')
                        batch_found = False
                        for part in parts:
                            if part.startswith('B') and len(part) >= 2 and part[1:].isdigit():
                                batch_ids.append(part)
                                batch_found = True
                                break
                            elif 'batch' in part.lower():
                                batch_ids.append(part)
                                batch_found = True
                                break
                        if not batch_found:
                            batch_ids.append('Batch1')
                    meta_df[col] = batch_ids

                elif col_lower in ['age', 'depth', 'weight', 'height', 'score']:
                    # Initialize numeric types with 0
                    meta_df[col] = 0

                elif col_lower in ['condition', 'treatment', 'group', 'status']:
                    # Use group name
                    meta_df[col] = group_names

                else:
                    # Use group name for others
                    meta_df[col] = group_names

            st.write("Edit metadata for each sample:")
            st.info("Default values have been auto-inferred from sample names. Edit as needed.")
            edited_meta_df = st.data_editor(meta_df, use_container_width=True)

            submitted = st.form_submit_button("Confirm Metadata", type="primary")

        if submitted:
            # Save to session_state
            st.session_state.metadata_submitted = True
            st.session_state.metadata_df = edited_meta_df.copy()
            st.session_state.meta_col_names = meta_col_names

            # Combine PCA scores and metadata
            if df.index.name or (df.index[0] != 0):
                df_combined = pd.concat([df, edited_meta_df], axis=1)
            else:
                # If first column is sample names
                df_temp = df.set_index(df.columns[0])
                df_combined = pd.concat([df_temp, edited_meta_df], axis=1)
                df_combined.reset_index(inplace=True)
                df_combined.rename(columns={'index': df.columns[0]}, inplace=True)

            st.session_state.combined_df = df_combined.copy()

        # Use combined dataframe if metadata has been confirmed
        if st.session_state.metadata_submitted and st.session_state.combined_df is not None:
            df = st.session_state.combined_df.copy()
            st.success("Metadata entry complete")

            st.write("Metadata confirmation:")
            for col in st.session_state.meta_col_names:
                st.write(f"**{col}**: " + ', '.join(df[col].unique().astype(str)))

    # =============================================================================
    # Column Mapping
    # =============================================================================

    step_number = "3" if input_method == "PCA scores file + Enter metadata manually" else "2"
    st.markdown(f"## {step_number}. Column Mapping")
    st.markdown("Select the columns corresponding to each variable")

    with st.form("column_mapping"):
        col1, col2, col3 = st.columns(3)

        with col1:
            st.markdown("### Required Columns")

            # Sample ID
            sample_col = st.selectbox("Sample ID Column:", df.columns,
                                       index=0 if 'sample' in df.columns[0].lower() else 0)

            # PC score columns
            numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
            pc_cols = [col for col in numeric_cols if col.upper().startswith('PC')]

            if len(pc_cols) == 0:
                st.warning("No columns starting with 'PC' found. Showing all numeric columns.")
                pc_cols = numeric_cols

            pc_col = st.selectbox("PC Score Column to Analyze:", pc_cols,
                                  index=0,
                                  help="Example: PC3, PC4, etc. For multiple PCs, run analysis one at a time")

        with col2:
            st.markdown("### Main Variables")

            # Get categorical columns
            cat_cols = df.select_dtypes(include=['object', 'category']).columns.tolist()

            # Exclude sample ID columns (columns with too many unique values)
            cat_cols_filtered = []
            for col in cat_cols:
                n_unique = df[col].nunique()
                n_samples = len(df)
                unique_ratio = n_unique / n_samples

                # Filtering criteria:
                # 1. Exclude sample IDs (all unique: ratio >= 0.9)
                # 2. Exclude columns with too many unique values (>50)
                # 3. Must have at least 2 unique values (exclude constant columns)
                if unique_ratio < 0.9 and n_unique >= 2 and n_unique <= 50:
                    cat_cols_filtered.append(col)

            if not cat_cols_filtered:
                st.warning("No suitable categorical variables found. Showing all categorical columns.")
                cat_cols_filtered = cat_cols

            # Main variable (e.g., sex)
            main_var = st.selectbox("Main Variable 1 (Required):",
                                    cat_cols_filtered,
                                    index=0 if len(cat_cols_filtered) > 0 else None,
                                    help="The main variable to test (e.g., sex, treatment, etc.).\n"
                                         "This variable is always included as a main effect in the model.\n"
                                         "Unique columns like sample ID are automatically excluded.")

            # Blocking variable (e.g., subtype)
            other_cats = [col for col in cat_cols_filtered if col != main_var]
            blocking_var = st.selectbox("Main Variable 2 / Blocking Variable (Optional):",
                                        ["None"] + other_cats,
                                        index=1 if len(other_cats) > 0 else 0,
                                        help="**Important:** This variable is also included as a main effect in the model.\n\n"
                                             "Example: If sex is Main Variable 1 and subtype is Main Variable 2,\n"
                                             "the model formula becomes `PC ~ C(sex) + C(subtype)`.\n\n"
                                             "**Use cases:**\n"
                                             "- When testing two main variables simultaneously (e.g., sex and subtype)\n"
                                             "- For blocking designs (stratification by donor, batch, etc.)\n\n"
                                             "Adding interaction terms allows testing the interaction between the two variables.")
            blocking_var = None if blocking_var == "None" else blocking_var

        with col3:
            st.markdown("### Hierarchical Structure")

            # Donor ID for random effects
            donor_col = st.selectbox("Donor/Subject ID (for Random Effects):",
                                     ["None", "Same as Sample ID"] + df.columns.tolist(),
                                     index=0,
                                     help="Specify when repeated measures exist. Used as (1|donor) random effect in LMM")

            if donor_col == "None":
                donor_col = None
            elif donor_col == "Same as Sample ID":
                donor_col = sample_col

            # Additional random effect
            additional_re = st.selectbox("Additional Random Effect (e.g., batch, library):",
                                         ["None"] + df.columns.tolist(),
                                         index=0,
                                         help="Use when treating technical batch effects as random effects")
            additional_re = None if additional_re == "None" else additional_re

        # Additional covariates
        st.markdown("### Additional Covariates (Fixed Effects)")

        # Exclude already used columns
        used_cols = [sample_col, pc_col, main_var, blocking_var, donor_col, additional_re]
        remaining_cols = [col for col in df.columns if col not in used_cols and col is not None]

        additional_covars = st.multiselect(
            "Select Additional Covariates (numeric or categorical):",
            remaining_cols,
            help="Variables to adjust for such as age, batch, depth. These are included as fixed effects in the model"
        )

        # Interaction terms
        st.markdown("### Interaction Terms")
        st.markdown("_Use interactions when the effects of two variables are not independent (e.g., sex differences vary by subtype)_")

        potential_interactions = []
        if blocking_var:
            potential_interactions.append(f"{main_var} : {blocking_var}")

        for covar in additional_covars:
            if df[covar].dtype in ['object', 'category']:
                potential_interactions.append(f"{main_var} : {covar}")

        interactions = st.multiselect(
            "Select Interaction Terms (Optional):",
            potential_interactions,
            help="Example: sex:subtype tests whether sex differences vary between subtypes. Use only with sufficient data"
        )

        # Form submit button
        mapping_submitted = st.form_submit_button("Confirm Column Mapping", type="primary")

    # =============================================================================
    # Run Analysis Button
    # =============================================================================

    # Show analysis button only when column mapping is confirmed
    if mapping_submitted or 'mapping_confirmed' in st.session_state:
        if mapping_submitted:
            # Clear all previous analysis results when column mapping is changed
            st.session_state.mapping_confirmed = True
            st.session_state.sample_col = sample_col
            st.session_state.pc_col = pc_col
            st.session_state.main_var = main_var
            st.session_state.blocking_var = blocking_var
            st.session_state.donor_col = donor_col
            st.session_state.additional_re = additional_re
            st.session_state.additional_covars = additional_covars
            st.session_state.interactions = interactions

            # Clear previous analysis results
            results_to_clear = [
                'ols_model', 'ols_coef', 'ols_anova', 'ols_emm',
                'lmm_result', 'lmm_coef',
                'perm_result',
                'fig_forest_ols', 'fig_diagnostic', 'fig_emm',
                'fig_forest_lmm', 'fig_null_dist'
            ]
            for key in results_to_clear:
                if key in st.session_state:
                    del st.session_state[key]

        st.success("Column mapping complete")

        # Get values from session_state
        sample_col = st.session_state.sample_col
        pc_col = st.session_state.pc_col
        main_var = st.session_state.main_var
        blocking_var = st.session_state.blocking_var
        donor_col = st.session_state.donor_col
        additional_re = st.session_state.additional_re
        additional_covars = st.session_state.additional_covars
        interactions = st.session_state.interactions

        if st.button("Run Analysis", type="primary"):

            step_number_check = "4" if input_method == "PCA scores file + Enter metadata manually" else "3"
            st.markdown(f"## {step_number_check}. Data Quality Check")

            # Prepare dataframe for analysis
            analysis_cols = [sample_col, pc_col, main_var]
            if blocking_var:
                analysis_cols.append(blocking_var)
            if donor_col and donor_col != sample_col:
                analysis_cols.append(donor_col)
            if additional_re:
                analysis_cols.append(additional_re)

            # Add additional_covars only if not empty
            if additional_covars:
                # Confirm it's a list before adding
                if isinstance(additional_covars, list):
                    analysis_cols.extend(additional_covars)
                else:
                    analysis_cols.append(additional_covars)

            # Exclude None and remove duplicates
            analysis_cols = [col for col in analysis_cols if col is not None]
            analysis_cols = list(dict.fromkeys(analysis_cols))  # Remove duplicates while preserving order

            # Debug information
            st.write("**Selected columns:**", analysis_cols)

            analysis_df = df[analysis_cols].copy()

            # Remove missing values
            n_before = len(analysis_df)
            analysis_df = analysis_df.dropna()
            n_after = len(analysis_df)

            if n_before > n_after:
                st.warning(f"Removed {n_before - n_after} rows containing missing values")

            # Standardize categorical variables
            for col in analysis_df.columns:
                # Verify col is a string
                if isinstance(col, str):
                    try:
                        if analysis_df[col].dtype == 'object':
                            analysis_df[col] = analysis_df[col].astype('category')
                    except (KeyError, AttributeError):
                        # Skip if column does not exist or has no type
                        continue

            # Check sample size (including missing data information)
            st.markdown("### Sample Size Summary")

            # Missing value check
            original_n = len(df)  # Original data count
            analysis_n = len(analysis_df)  # Data count used for analysis
            excluded_n = original_n - analysis_n

            if excluded_n > 0:
                st.warning(f"**Missing Value Exclusion**: {excluded_n} samples were excluded due to missing values (Original data: {original_n} -> Used for analysis: {analysis_n})")

            col1, col2 = st.columns(2)

            with col1:
                st.metric("Sample Size Used for Analysis", len(analysis_df))

                # Distribution of main variable
                st.write(f"**Distribution of {main_var}:**")
                main_var_counts = analysis_df[main_var].value_counts()
                st.dataframe(main_var_counts)

                # Imbalance warning
                if len(main_var_counts) > 1:
                    max_count = main_var_counts.max()
                    min_count = main_var_counts.min()
                    if max_count / min_count > 3:
                        st.warning(f"⚠️ **Group Imbalance**: More than 3-fold difference between the largest group ({max_count}) and smallest group ({min_count})")

            with col2:
                if blocking_var:
                    st.write(f"**Cross-tabulation: {main_var} x {blocking_var}**")
                    cross_tab = pd.crosstab(analysis_df[main_var], analysis_df[blocking_var])
                    st.dataframe(cross_tab)

                    # Small sample size warning
                    min_cell_size = cross_tab.min().min()
                    if min_cell_size < 3:
                        st.warning(f"⚠️ Some cells have very small sample sizes (minimum: {min_cell_size}). Estimates may be unstable.")

                if donor_col and donor_col != sample_col:
                    n_donors = analysis_df[donor_col].nunique()
                    st.metric("Number of Donors", n_donors)

                    if n_donors < 5:
                        st.warning("⚠️ Number of donors is less than 5. Variance estimation for random effects may be unstable.")

            # Collinearity check
            st.markdown("### Collinearity Diagnostics")

            formula_terms = [main_var]
            if blocking_var:
                formula_terms.append(blocking_var)
            formula_terms.extend(additional_covars)

            collinearity_issues = detect_collinearity(analysis_df, formula_terms)

            if collinearity_issues:
                st.error("**Perfect collinearity detected:**")
                for issue in collinearity_issues:
                    st.markdown(issue)
                st.warning("""
                **Perfectly collinear variables cannot have their effects estimated separately.**

                Example: If a certain subtype is all Female, then the sex effect and subtype effect are indistinguishable.

                **Solutions:**
                - Remove one of the variables
                - Collect more balanced data
                - Avoid using interaction terms and limit to descriptive comparisons
                """)
            else:
                st.success("✅ No perfect collinearity detected")

            # VIF for numeric variables
            numeric_covars = [col for col in additional_covars
                             if pd.api.types.is_numeric_dtype(analysis_df[col])]

            if len(numeric_covars) > 1:
                st.markdown("**Multicollinearity of Numeric Covariates (VIF):**")
                st.caption("VIF > 10: Strong multicollinearity. VIF > 5: Moderate multicollinearity.")

                try:
                    vif_df = calculate_vif(analysis_df, numeric_covars)
                    vif_df_display = vif_df.copy(); vif_df_display['VIF'] = vif_df_display['VIF'].round(2); st.dataframe(vif_df_display, use_container_width=True)

                    if (vif_df['VIF'] > 10).any():
                        st.warning("⚠️ Variables with VIF > 10 detected. Strong multicollinearity is suggested.")
                    elif (vif_df['VIF'] > 5).any():
                        st.info("ℹ️ Variables with VIF > 5 detected. Moderate multicollinearity is suggested.")
                except Exception as e:
                    st.info(f"VIF calculation skipped: {str(e)}")

            # =============================================================================
            # Build Formula
            # =============================================================================

            st.markdown("## 4. Model Formula Construction")

            # Debug information
            with st.expander("Model Construction Info (Selected Variables Check)", expanded=True):
                st.write(f"**Main Variable 1 (main_var):** {main_var}")
                st.write(f"**Main Variable 2 / Blocking Variable (blocking_var):** {blocking_var}")
                st.write(f"**Additional Covariates (additional_covars):** {additional_covars}")
                st.write(f"**Interaction Terms (interactions):** {interactions}")
                st.info("To include variables in the model formula, select them above under 'Main Variable 1' or 'Main Variable 2'.\n\n"
                        "'Additional Covariates' are used as adjustment variables and are not the main test targets.")

            # Fixed effects formula
            fixed_terms = [f"C({main_var})"]
            if blocking_var:
                fixed_terms.append(f"C({blocking_var})")

            for col in additional_covars:
                if analysis_df[col].dtype in ['object', 'category']:
                    fixed_terms.append(f"C({col})")
                else:
                    fixed_terms.append(col)

            # Add interaction terms
            if interactions:
                for interaction in interactions:
                    vars_in_interaction = [v.strip() for v in interaction.split(":")]
                    interaction_term = " * ".join([
                        f"C({v})" if analysis_df[v].dtype in ['object', 'category'] else v
                        for v in vars_in_interaction
                    ])
                    # Add only if not already included
                    if interaction_term not in fixed_terms:
                        fixed_terms.append(interaction_term)

            formula = f"{pc_col} ~ " + " + ".join(fixed_terms)

            # Explicitly display variables included in the model
            st.markdown("**Variables Included in Model:**")
            col_a, col_b = st.columns(2)
            with col_a:
                st.success(f"Main Variable 1: **{main_var}**")
                if blocking_var:
                    st.success(f"Main Variable 2: **{blocking_var}**")
                else:
                    st.warning("Main Variable 2: Not selected")
            with col_b:
                if additional_covars:
                    st.info(f"Adjustment Variables: {', '.join(additional_covars)}")
                if interactions:
                    st.info(f"Interactions: {', '.join(interactions)}")

            st.markdown("**Model Formula:**")
            st.code(formula, language="r")

            if donor_col:
                st.markdown(f"**Random Effects (LMM only):** `(1|{donor_col})`")

            # =============================================================================
            # Run OLS
            # =============================================================================

            if "OLS" in analysis_method:
                st.markdown("---")
                st.markdown("## 5. OLS Analysis Results")

                try:
                    # WLS vs OLS
                    if use_wls:
                        # First fit OLS to get residuals for weights
                        ols_temp = smf.ols(formula, data=analysis_df).fit()
                        # Calculate weights as inverse of squared residuals
                        weights = 1.0 / (ols_temp.resid ** 2)
                        # Normalize weights
                        weights = weights / weights.mean()
                        # Fit WLS
                        ols_model = smf.wls(formula, data=analysis_df, weights=weights).fit()
                        st.info(f"ℹ️ Using WLS (Weighted Least Squares). Weights are calculated as the inverse variance of residuals.")
                    else:
                        # Fit model (OLS)
                        ols_model = smf.ols(formula, data=analysis_df).fit()

                    # Get appropriate covariance matrix
                    if se_type == "HC3 (Robust, recommended)":
                        ols_model = ols_model.get_robustcov_results(cov_type='HC3')
                        st.info("ℹ️ Using HC3 robust standard errors (handles heteroscedasticity)")
                    elif se_type == "Cluster-robust":
                        cluster_var = donor_col if donor_col else additional_re
                        if cluster_var:
                            ols_model = ols_model.get_robustcov_results(
                                cov_type='cluster',
                                groups=analysis_df[cluster_var]
                            )
                            st.info(f"ℹ️ Using cluster-robust standard errors (Cluster variable: {cluster_var})")
                        else:
                            st.warning("Cluster variable not specified. Using standard HC3.")
                            ols_model = ols_model.get_robustcov_results(cov_type='HC3')

                    # Explicitly show model formula and reference levels
                    st.markdown("### Coefficient Table")
                    st.caption(f"**Model Formula**: `{formula}`")

                    # Extract and display reference levels
                    reference_levels = {}
                    for col in analysis_df.columns:
                        if analysis_df[col].dtype == 'object' or analysis_df[col].dtype.name == 'category':
                            if col in formula:
                                ref_level = analysis_df[col].iloc[0] if len(analysis_df) > 0 else "N/A"
                                # Patsy uses the first value alphabetically as reference
                                unique_vals = sorted(analysis_df[col].unique())
                                if len(unique_vals) > 0:
                                    reference_levels[col] = unique_vals[0]

                    if reference_levels:
                        ref_info = ", ".join([f"`{var}` reference={ref}" for var, ref in reference_levels.items()])
                        st.caption(f"**Reference Levels**: {ref_info}")
                    conf_int = ols_model.conf_int()

                    # Determine if conf_int is DataFrame or ndarray
                    if isinstance(conf_int, pd.DataFrame):
                        ci_lower = conf_int.iloc[:, 0]
                        ci_upper = conf_int.iloc[:, 1]
                    else:
                        # If numpy array
                        ci_lower = conf_int[:, 0]
                        ci_upper = conf_int[:, 1]

                    coef_table = pd.DataFrame({
                        'Coef': ols_model.params,
                        'Std_Error': ols_model.bse,
                        't_value': ols_model.tvalues,
                        'p_value': ols_model.pvalues,
                        'CI_lower': ci_lower,
                        'CI_upper': ci_upper
                    })

                    # Add significance symbols
                    def sig_symbol(p):
                        if p < 0.001:
                            return '***'
                        elif p < 0.01:
                            return '**'
                        elif p < 0.05:
                            return '*'
                        elif p < 0.1:
                            return '.'
                        else:
                            return ''

                    coef_table['Sig'] = coef_table['p_value'].apply(sig_symbol)

                    # Round numbers and create display DataFrame
                    display_table = coef_table.copy()
                    display_table['Coef'] = display_table['Coef'].round(4)
                    display_table['Std_Error'] = display_table['Std_Error'].round(4)
                    display_table['t_value'] = display_table['t_value'].round(3)
                    display_table['p_value'] = display_table['p_value'].apply(lambda x: f"{x:.4g}")
                    display_table['CI_lower'] = display_table['CI_lower'].round(4)
                    display_table['CI_upper'] = display_table['CI_upper'].round(4)

                    st.dataframe(display_table, use_container_width=True)

                    st.caption("Significance levels: *** p<0.001, ** p<0.01, * p<0.05, . p<0.1")

                    # Coefficient table interpretation guide
                    with st.expander("How to Read the Coefficient Table (Meaning of p_value Column)"):
                        st.markdown("""
                        ### Meaning of p_value Column

                        **Shows the significance of each coefficient (after adjusting for other variables)**

                        - **Covariate-adjusted**: Results considering all variables simultaneously
                        - This p-value tests the "pure effect" after removing the influence of other variables

                        ---

                        ### Notation Meaning
                        """)

                        # Generate explanations from variables included in the actual model
                        example_vars = []
                        for idx in display_table.index:
                            idx_str = str(idx)  # Ensure idx is a string
                            if idx_str != 'Intercept' and 'C(' in idx_str:
                                example_vars.append(idx_str)

                        if example_vars:
                            st.markdown("**Examples from this model:**")
                            # Use the first 1-2 variables as examples
                            for var_name in example_vars[:2]:
                                # Parse format like C(sex)[T.M]
                                if '[T.' in var_name:
                                    var_part = var_name.split('[')[0]  # C(sex)
                                    level_part = var_name.split('[T.')[1].rstrip(']')  # M
                                    base_var = var_part.replace('C(', '').replace(')', '')  # sex

                                    st.markdown(f"""
                                    - **`{var_name}`**
                                      - `C({base_var})`: Treat {base_var} as a categorical variable
                                      - `[T.{level_part}]`: Represents the {level_part} level in Treatment coding
                                      - **Meaning**: "How much does {level_part} differ from the reference category"
                                      - **p-value interpretation**: Is there a significant difference between {level_part} and the reference category (**adjusted for other variables**)
                                    """)

                        st.markdown("""
                        ---

                        ### General Notation

                        | Notation | Meaning |
                        |-----|------|
                        | `C(variable_name)` | Treat as a categorical variable |
                        | `[T.level_name]` | Represents the specified level in Treatment coding |
                        | `C(variable_name)[T.level_name]` | "How much does this level differ from the reference category" |

                        ### Reference Category (Reference level)

                        - Default: **First value in alphabetical order**
                        - The reference category row is not displayed in the coefficient table (coefficient = 0)
                        - Other levels are interpreted as "difference from the reference"

                        ### Difference from Welch t-test

                        | Test Method | Covariate Adjustment | Display Location |
                        |---------|----------|---------|
                        | Welch t-test | No (simple two-group comparison) | Info box |
                        | Coefficient table p_value | Yes (considering other variables) | This table |

                        **Recommended**: Use the p-values from this coefficient table for accurate statistical analysis.
                        """)

                    # Model fit
                    st.markdown("### 📈 Model Fit")
                    col1, col2, col3, col4 = st.columns(4)
                    with col1:
                        st.metric("R²", f"{ols_model.rsquared:.4f}")
                    with col2:
                        st.metric("Adjusted R²", f"{ols_model.rsquared_adj:.4f}")
                    with col3:
                        st.metric("AIC", f"{ols_model.aic:.2f}")
                    with col4:
                        st.metric("BIC", f"{ols_model.bic:.2f}")

                    # ANOVA
                    try:
                        if anova_type == "Type II (recommended)":
                            st.markdown("### 📊 Type II ANOVA")
                            st.caption("Tests overall significance of each variable (contribution in a model including other variables)")
                            anova_table = anova_lm(ols_model, typ=2)
                        else:  # Type III
                            st.markdown("### 📊 Type III ANOVA")
                            st.caption("Tests significance of each variable including all interactions")
                            anova_table = anova_lm(ols_model, typ=3)

                        # Display important notes (depending on SE type)
                        if se_type == "HC3 (Robust, recommended)":
                            st.warning("**ANOVA F-test assumes homoscedasticity and normality**. The p-values have different meanings from the coefficient table above (HC3 robust SE). If heteroscedasticity is a concern, prioritize the p-values from the coefficient table.")
                        elif se_type == "Cluster-robust":
                            st.warning("**ANOVA F-test assumes homoscedasticity and normality**. The p-values have different meanings from the coefficient table above (Cluster robust SE). If within-cluster correlation or heteroscedasticity is a concern, prioritize the p-values from the coefficient table.")
                        elif se_type == "Classical":
                            st.info("**ANOVA F-test and coefficient table use the same assumptions** (homoscedasticity and normality). If heteroscedasticity is a concern, consider using HC3 robust SE.")

                        # Round values for display
                        anova_display = anova_table.copy()
                        if 'sum_sq' in anova_display.columns:
                            anova_display['sum_sq'] = anova_display['sum_sq'].round(4)
                        if 'F' in anova_display.columns:
                            anova_display['F'] = anova_display['F'].round(3)
                        if 'PR(>F)' in anova_display.columns:
                            anova_display['PR(>F)'] = anova_display['PR(>F)'].apply(lambda x: f"{x:.4g}" if pd.notna(x) else x)

                        st.dataframe(anova_display, use_container_width=True)
                        st.session_state.ols_anova = anova_table  # Save for zip download

                        # ANOVA table interpretation guide
                        with st.expander("How to read the ANOVA table (Meaning of PR(>F) column)"):
                            st.markdown("""
                            ### Meaning of PR(>F) Column

                            **Shows the overall significance of each variable (adjusted for other variables)**

                            - **With covariate adjustment**: Results considering all variables simultaneously
                            - This p-value tests whether the entire variable significantly contributes after removing the influence of other variables

                            ---

                            ### Difference from Coefficient Table

                            | Table | Test Target | Covariate Adjustment |
                            |----|---------|----------|
                            | **Coefficient table** | Each level (e.g., Male vs Female) | Yes |
                            | **ANOVA table** | Entire variable (e.g., overall effect of sex) | Yes |
                            """)

                            # Generate explanations from variables included in the actual model
                            anova_vars = []
                            for idx in anova_display.index:
                                idx_str = str(idx)  # Ensure idx is a string
                                if idx_str != 'Residual' and 'C(' in idx_str:
                                    anova_vars.append(idx_str)

                            if anova_vars:
                                st.markdown("**Examples from this model:**")
                                for var_name in anova_vars[:2]:
                                    # Extract variable name from format like C(sex)
                                    base_var = var_name.replace('C(', '').replace(')', '')

                                    st.markdown(f"""
                                    - PR(>F) for **`{var_name}`**
                                      - **Meaning**: Does the entire {base_var} variable significantly affect the PC value
                                      - **Test content**: Whether all coefficients for {base_var} levels are simultaneously zero
                                      - **Adjustment**: Test after removing the influence of other variables (secondary variable, covariates, etc.)
                                    """)

                            st.markdown("""
                            ---

                            ### Understanding with a Concrete Example

                            **Model**: `PC4 ~ C(sex) + C(subtype)`

                            | Variable | PR(>F) | Meaning |
                            |------|--------|------|
                            | C(sex) | 0.0234 | After adjusting for subtype, sex significantly affects PC4 |
                            | C(subtype) | 0.0567 | After adjusting for sex, subtype significantly affects PC4 |

                            ### Difference from Welch t-test

                            | Test Method | Covariate Adjustment | Test Target | Display Location |
                            |---------|----------|---------|---------|
                            | Welch t-test | No | Primary variable 1 only (two-group comparison) | Info box |
                            | ANOVA table PR(>F) | Yes | Each entire variable | This table |

                            ### Notes

                            **ANOVA F-test assumes homoscedasticity and normality**

                            If heteroscedasticity is a concern, prioritize the p-values from the coefficient table (when using HC3/Cluster robust SE).

                            ### Type II vs Type III

                            - **Type II ANOVA (Recommended)**: Tests main effects of each variable in a model without interactions
                            - **Type III ANOVA**: Tests effects of each variable including all interactions

                            **Recommended**: Use PR(>F) from this ANOVA table for accurate overall variable significance.
                            """)

                    except Exception as e:
                        st.warning(f"Could not calculate ANOVA table: {str(e)}")

                    # Forest plot
                    st.markdown("### 🌲 Forest Plot (Coefficient Visualization)")
                    st.caption("Estimated values and 95% confidence intervals for each coefficient. Significant if not overlapping zero")

                    # Exclude intercept
                    coef_for_plot = coef_table.iloc[1:].copy()
                    if len(coef_for_plot) > 0:
                        fig_forest = plot_forest(coef_for_plot, title=f"OLS Coefficient Estimates ({se_type})")
                        st.pyplot(fig_forest)
                        st.session_state.fig_forest_ols = fig_forest  # Save for zip download
                    else:
                        st.info("No coefficients available for plotting (intercept only)")

                    # Diagnostic plots
                    if show_diagnostics:
                        st.markdown("### 🔬 Diagnostic Plots")
                        st.caption("""
                        - **Residuals vs Fitted**: Linearity and homoscedasticity are satisfied if no pattern is present
                        - **Normal Q-Q**: Normality is satisfied if points fall on a straight line
                        - **Scale-Location**: Homoscedasticity is satisfied if it shows a horizontal band
                        - **Residual Distribution**: Check if close to normal distribution
                        """)
                        fig_diag = plot_diagnostic(ols_model, title="OLS Model Diagnostics")
                        st.pyplot(fig_diag)
                        st.session_state.fig_diagnostic = fig_diag  # Save for zip download

                    # EMM
                    if show_emm and blocking_var:
                        st.markdown("### 📊 Estimated Marginal Means (EMM)")
                        st.caption("Predicted mean values for each group after adjusting for covariates")

                        try:
                            emm_df = calculate_emm(ols_model, analysis_df,
                                                  [main_var, blocking_var] if blocking_var else [main_var])

                            # Round values for display
                            emm_display = emm_df.copy()
                            for col in ['mean', 'se', 'ci_lower', 'ci_upper']:
                                if col in emm_display.columns:
                                    emm_display[col] = emm_display[col].round(4)

                            st.dataframe(emm_display, use_container_width=True)

                            # EMM plot
                            fig, ax = plt.subplots(figsize=(10, 6))

                            unique_blocks = emm_df[blocking_var].unique()
                            colors = plt.cm.Set2(np.linspace(0, 1, len(unique_blocks)))

                            for idx, block in enumerate(unique_blocks):
                                block_data = emm_df[emm_df[blocking_var] == block]
                                x_labels = block_data[main_var].values
                                x_pos = np.arange(len(block_data)) + idx * 0.2

                                ax.errorbar(x_pos, block_data['mean'],
                                           yerr=[block_data['mean'] - block_data['ci_lower'],
                                                 block_data['ci_upper'] - block_data['mean']],
                                           label=str(block), fmt='o-', capsize=5,
                                           markersize=8, linewidth=2, color=colors[idx])

                            ax.set_xticks(np.arange(len(x_labels)) + (len(unique_blocks)-1) * 0.1)
                            ax.set_xticklabels(x_labels)
                            ax.set_xlabel(main_var, fontsize=12)
                            ax.set_ylabel(f'Estimated Mean {pc_col}', fontsize=12)
                            ax.set_title('Estimated Marginal Means (95% CI)', fontsize=13, fontweight='bold')
                            ax.legend(title=blocking_var, fontsize=10)
                            ax.grid(True, alpha=0.3)
                            plt.tight_layout()
                            st.pyplot(fig)
                            st.session_state.fig_emm = fig  # Save for zip download
                            st.session_state.ols_emm = emm_df  # Save data as well

                        except Exception as e:
                            st.warning(f"Error occurred during EMM calculation: {str(e)}")

                    # Save results
                    st.session_state.ols_model = ols_model
                    st.session_state.ols_coef = coef_table

                except Exception as e:
                    st.error(f"OLS analysis failed")
                    st.exception(e)

            # =============================================================================
            # Run LMM
            # =============================================================================

            if "LMM" in analysis_method and donor_col:
                st.markdown("---")
                st.markdown("## 6️⃣ Linear Mixed Model (LMM) Analysis Results")

                try:
                    # Prepare random effects formula
                    re_formula = "1"  # Random intercept
                    groups = analysis_df[donor_col]

                    # Fit LMM
                    lmm_model = smf.mixedlm(formula, data=analysis_df,
                                            groups=groups,
                                            re_formula=re_formula)

                    lmm_result = lmm_model.fit(reml=(reml_method == "REML (Recommended)"))

                    st.info(f"Estimation method: {reml_method}")

                    # Coefficient table
                    st.markdown("### 📋 Fixed Effects Coefficient Table")
                    coef_table_lmm = pd.DataFrame({
                        'Coef': lmm_result.params,
                        'Std_Error': lmm_result.bse,
                        'z_value': lmm_result.tvalues,
                        'p_value': lmm_result.pvalues,
                        'CI_lower': lmm_result.conf_int()[0],
                        'CI_upper': lmm_result.conf_int()[1]
                    })

                    coef_table_lmm['Sig'] = coef_table_lmm['p_value'].apply(
                        lambda p: '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else '.' if p < 0.1 else ''
                    )

                    # Round values for display
                    lmm_display = coef_table_lmm.copy()
                    lmm_display['Coef'] = lmm_display['Coef'].round(4)
                    lmm_display['Std_Error'] = lmm_display['Std_Error'].round(4)
                    lmm_display['z_value'] = lmm_display['z_value'].round(3)
                    lmm_display['p_value'] = lmm_display['p_value'].apply(lambda x: f"{x:.4g}")
                    lmm_display['CI_lower'] = lmm_display['CI_lower'].round(4)
                    lmm_display['CI_upper'] = lmm_display['CI_upper'].round(4)

                    st.dataframe(lmm_display, use_container_width=True)

                    st.caption("Significance levels: *** p<0.001, ** p<0.01, * p<0.05, . p<0.1")

                    # Variance components
                    st.markdown("### 📊 Variance Components")
                    st.caption("Shows where the variation in the data comes from")

                    var_random = lmm_result.cov_re.values[0][0]
                    var_residual = lmm_result.scale
                    var_total = var_random + var_residual

                    var_comp = pd.DataFrame({
                        'Component': [f'Random effect ({donor_col})', 'Residual (within-subject variation)'],
                        'Variance': [var_random, var_residual],
                        'Std_Dev': [np.sqrt(var_random), np.sqrt(var_residual)],
                        'Proportion': [var_random / var_total, var_residual / var_total]
                    })

                    # Round values for display
                    var_comp_display = var_comp.copy()
                    var_comp_display['Variance'] = var_comp_display['Variance'].round(4)
                    var_comp_display['Std_Dev'] = var_comp_display['Std_Dev'].round(4)
                    var_comp_display['Proportion'] = (var_comp_display['Proportion'] * 100).round(2).astype(str) + '%'

                    st.dataframe(var_comp_display, use_container_width=True)

                    # ICC
                    icc = var_random / var_total
                    st.metric("Intraclass Correlation Coefficient (ICC)", f"{icc:.4f}",
                             help="Correlation between samples within the same donor. Close to 0 means small random effect, close to 1 means large between-donor variation")

                    if icc < 0.05:
                        st.info("ICC < 0.05: Random effect is very small. OLS may be sufficient.")
                    elif icc > 0.5:
                        st.success("ICC > 0.5: Random effect is large. Using LMM is appropriate.")

                    # Model fit
                    st.markdown("### 📈 Model Fit")
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Log Likelihood", f"{lmm_result.llf:.2f}")
                    with col2:
                        st.metric("AIC", f"{lmm_result.aic:.2f}")
                    with col3:
                        st.metric("BIC", f"{lmm_result.bic:.2f}")

                    # Forest plot
                    st.markdown("### 🌲 Forest Plot (Fixed Effects)")
                    coef_for_plot_lmm = coef_table_lmm.iloc[1:].copy()
                    if len(coef_for_plot_lmm) > 0:
                        fig_forest_lmm = plot_forest(coef_for_plot_lmm,
                                                     title="LMM Fixed Effects Estimates")
                        st.pyplot(fig_forest_lmm)
                        st.session_state.fig_forest_lmm = fig_forest_lmm  # Save for zip download

                    # Save results
                    st.session_state.lmm_result = lmm_result
                    st.session_state.lmm_coef = coef_table_lmm

                except Exception as e:
                    st.error(f"LMM analysis failed")
                    st.exception(e)
                    st.info("Hint: LMM convergence may fail when the number of donors is small or data variation is small. Try OLS instead.")

            elif "LMM" in analysis_method and not donor_col:
                st.warning("To run LMM, please specify Donor/Subject ID.")

            # =============================================================================
            # Welch's t-test
            # =============================================================================

            if run_welch:
                st.markdown("---")
                st.markdown("## 7️⃣ Welch's t-test Results")
                st.caption("Two-group comparison without assuming equal variance (no covariate adjustment)")
                st.warning(f"**Test target variable:** Primary variable 1 = **{main_var}**")

                # Check number of levels in primary variable
                unique_main_levels = analysis_df[main_var].nunique()
                if unique_main_levels == 2:
                    from scipy import stats
                    groups = analysis_df[main_var].unique()
                    group1_data = analysis_df[analysis_df[main_var] == groups[0]][pc_col]
                    group2_data = analysis_df[analysis_df[main_var] == groups[1]][pc_col]
                    welch_t, welch_p = stats.ttest_ind(group1_data, group2_data, equal_var=False)

                    # Descriptive statistics
                    group1_mean = group1_data.mean()
                    group1_std = group1_data.std()
                    group1_n = len(group1_data)
                    group2_mean = group2_data.mean()
                    group2_std = group2_data.std()
                    group2_n = len(group2_data)

                    # Effect size (Cohen's d, Hedges' g)
                    pooled_std = np.sqrt(((group1_n - 1) * group1_std**2 + (group2_n - 1) * group2_std**2) / (group1_n + group2_n - 2))
                    cohens_d = (group1_mean - group2_mean) / pooled_std
                    correction_factor = 1 - (3 / (4 * (group1_n + group2_n) - 9))
                    hedges_g = cohens_d * correction_factor

                    # Display results
                    st.success("Welch's t-test completed")

                    col1, col2 = st.columns(2)

                    with col1:
                        st.markdown("### 📊 Descriptive Statistics")
                        desc_df = pd.DataFrame({
                            'Group': [groups[0], groups[1]],
                            'n': [group1_n, group2_n],
                            'Mean': [f"{group1_mean:.4f}", f"{group2_mean:.4f}"],
                            'Std_Dev': [f"{group1_std:.4f}", f"{group2_std:.4f}"]
                        })
                        st.dataframe(desc_df, use_container_width=True)

                    with col2:
                        st.markdown("### 📈 Test Results")
                        st.metric("t-value", f"{welch_t:.4f}")
                        st.metric("p-value", f"{welch_p:.4g}")
                        if welch_p < 0.001:
                            st.success("***  p < 0.001")
                        elif welch_p < 0.01:
                            st.success("**  p < 0.01")
                        elif welch_p < 0.05:
                            st.success("*  p < 0.05")
                        elif welch_p < 0.1:
                            st.info(".  p < 0.1")
                        else:
                            st.info("n.s.  p ≥ 0.1")

                    st.markdown("### 📏 Effect Size")
                    effect_col1, effect_col2 = st.columns(2)
                    with effect_col1:
                        st.metric("Cohen's d", f"{cohens_d:.4f}")
                    with effect_col2:
                        st.metric("Hedges' g (corrected)", f"{hedges_g:.4f}")

                    st.caption("Effect size guidelines: |d| < 0.2 (small), 0.2-0.5 (medium), 0.5-0.8 (large), >= 0.8 (very large)")

                    # Important notes
                    st.warning(
                        f"**Important Note**\n\n"
                        f"This test is a **simple two-group comparison**:\n"
                        f"- Effects of other variables (secondary variable, additional covariates, etc.) are **NOT adjusted**\n"
                        f"- If covariate adjustment is needed, please check the OLS/LMM 'Coefficient table' or 'ANOVA table'\n\n"
                        f"**Usage guidelines**:\n"
                        f"- Welch's t-test: Simple comparison of primary variable only (reference value)\n"
                        f"- OLS coefficient table: Effect estimation after covariate adjustment (main analysis)\n"
                        f"- Permutation test: Robust test without distribution assumptions (supplementary)"
                    )

                    # Save to session_state
                    st.session_state.welch_result = {
                        't': welch_t,
                        'p': welch_p,
                        'cohens_d': cohens_d,
                        'hedges_g': hedges_g,
                        'group1': groups[0],
                        'group2': groups[1],
                        'group1_mean': group1_mean,
                        'group2_mean': group2_mean
                    }

                else:
                    st.error(f"Welch's t-test only supports two-group comparison. Primary variable 1 ({main_var}) has {unique_main_levels} levels.")
                    st.info("For 3 or more groups, please use the OLS ANOVA table or permutation test.")

            # =============================================================================
            # Permutation Test
            # =============================================================================

            if run_permutation:
                st.markdown("---")
                section_number = "8️⃣" if run_welch else "7️⃣"
                st.markdown(f"## {section_number} Permutation Test Results")
                st.caption("A robust testing method that does not rely on parametric assumptions")
                st.warning(f"🎯 **Variable being tested:** Primary variable 1 = **{main_var}**")

                with st.spinner(f"Running {n_permutations:,} permutations..."):
                    try:
                        # Prepare design matrix (get design_info)
                        y = analysis_df[pc_col].values
                        X_full = dmatrix(formula.split('~')[1], data=analysis_df, return_type='dataframe')
                        design_info = X_full.design_info  # Get patsy design_info

                        # Prepare stratification variable
                        stratify_var = None
                        if blocking_var is not None:
                            stratify_var = analysis_df[blocking_var].values
                            st.warning(f"🔹 **Using stratified permutation**: Residuals are permuted only within {blocking_var}\n\n"
                                      f"This generates a null distribution for the {main_var} effect while preserving the {blocking_var} structure")
                        else:
                            st.info("ℹ️ Using non-stratified permutation because primary variable 2 is not specified")

                        # One-sided test settings
                        one_sided_param = None
                        if perm_sided == "One-sided: Main Variable 1 > Reference (greater)":
                            one_sided_param = 'greater'
                            st.info(f"📊 **One-sided test (greater)**: Testing if {main_var} level 1 > reference (level 0)")
                        elif perm_sided == "One-sided: Main Variable 1 < Reference (less)":
                            one_sided_param = 'less'
                            st.info(f"📊 **One-sided test (less)**: Testing if {main_var} level 1 < reference (level 0)")

                        # Run permutation test
                        if perm_method == "Freedman-Lane (Recommended)":
                            perm_result = freedman_lane_permutation(
                                y, X_full, main_var,
                                n_perm=n_permutations,
                                stratify_by=stratify_var,
                                one_sided=one_sided_param,
                                design_info=design_info  # Added for robust column identification
                            )

                            info_msg = f"ℹ️ **Freedman-Lane method:** Testing the effect of primary variable 1 (**{main_var}**)\n\n"
                            info_msg += f"Other variables (primary variable 2, additional covariates, etc.) are treated as control variables\n\n"
                            if perm_result.get('stratified'):
                                info_msg += f"✅ Stratified permutation: Preserving {blocking_var} structure\n\n"
                            if perm_result.get('one_sided'):
                                info_msg += f"📊 Test direction: {perm_result['one_sided']}\n\n"
                            st.success(info_msg)
                        else:
                            # Simple label permutation (y-value permutation version: fast)
                            np.random.seed(42)  # Use the same seed as Freedman-Lane

                            # Observed statistic
                            obs_model = sm.OLS(y, X_full).fit()
                            interest_cols = [col for col in X_full.columns
                                           if f'C({main_var})' in col and col != 'Intercept']
                            main_effect_col = [col for col in interest_cols if ':' not in col][0]
                            obs_stat = float(obs_model.tvalues[main_effect_col])

                            null_dist = []
                            for _ in range(n_permutations):
                                # Permute y values (with stratification support)
                                if stratify_var is not None:
                                    # Stratified permutation: shuffle y values only within each stratum
                                    y_perm = np.zeros_like(y)
                                    for stratum in np.unique(stratify_var):
                                        stratum_mask = (stratify_var == stratum)
                                        stratum_indices = np.where(stratum_mask)[0]
                                        perm_indices = np.random.permutation(stratum_indices)
                                        y_perm[stratum_indices] = y[perm_indices]
                                else:
                                    # Non-stratified permutation: shuffle all
                                    y_perm = np.random.permutation(y)

                                # Fit (design matrix is fixed)
                                perm_model = sm.OLS(y_perm, X_full).fit()
                                null_dist.append(float(perm_model.tvalues[main_effect_col]))

                            null_dist = np.array(null_dist)

                            # Calculate p-value (with continuity correction)
                            if one_sided_param == 'greater':
                                p_value = (np.sum(null_dist >= obs_stat) + 1) / (n_permutations + 1)
                            elif one_sided_param == 'less':
                                p_value = (np.sum(null_dist <= obs_stat) + 1) / (n_permutations + 1)
                            else:
                                p_value = (np.sum(np.abs(null_dist) >= np.abs(obs_stat)) + 1) / (n_permutations + 1)

                            perm_result = {
                                'observed': obs_stat,
                                'null_distribution': null_dist,
                                'p_value': p_value,
                                'stratified': stratify_var is not None,
                                'one_sided': one_sided_param
                            }

                            info_msg = f"ℹ️ **Simple permutation:** Permute y values (fast version)\n\n"
                            info_msg += "📊 **Difference in null hypotheses:**\n"
                            info_msg += "- **Freedman-Lane**: Null hypothesis conditioned on covariates (permute residuals)\n"
                            info_msg += "- **Simple**: Null hypothesis ignoring covariate structure (permute y values)\n\n"
                            info_msg += "⚠️ **p-value differences**: It is normal for p-values to differ because the null hypotheses are different.\n"
                            if stratify_var is not None:
                                info_msg += f"✅ Stratified permutation: y values are swapped only within {blocking_var}\n\n"
                            else:
                                info_msg += "y values are swapped across all samples\n\n"
                            info_msg += "💡 **Recommendation**: Use the **Freedman-Lane method** if you want to control for covariates."
                            st.warning(info_msg)

                        # Display results
                        col1, col2, col3, col4 = st.columns(4)
                        with col1:
                            st.metric("Observed statistic (t-value)", f"{perm_result['observed']:.4f}")
                        with col2:
                            test_type = "One-sided" if perm_result.get('one_sided') else "Two-sided"
                            st.metric("Test type", test_type)
                        with col3:
                            st.metric("Permutation p-value", f"{perm_result['p_value']:.4g}")
                        with col4:
                            sig_level = "***" if perm_result['p_value'] < 0.001 else \
                                        "**" if perm_result['p_value'] < 0.01 else \
                                        "*" if perm_result['p_value'] < 0.05 else "n.s."
                            st.metric("Significance", sig_level)

                        # Display stratification information
                        if perm_result.get('stratified'):
                            st.success(f"✅ **Stratified permutation**: Residuals were permuted only within each level of {blocking_var}\n\n"
                                      f"Generated null distribution by nullifying {main_var} effect within each {blocking_var}")
                        else:
                            st.info("ℹ️ **Non-stratified permutation**: Residuals were permuted across all samples")

                        # Display reproducibility information (method_info)
                        if 'method_info' in perm_result:
                            with st.expander("🔬 Analysis details (information for reproducibility)", expanded=False):
                                info = perm_result['method_info']
                                st.markdown(f"""
**Model Formula:**
- **Null model**: `{info['reduced_model']}`
- **Alternative model**: `{info['full_model']}`

**Test settings:**
- **Number of permutations**: {info['n_permutations']:,}
- **Random seed**: {info['random_seed']}
- **Test statistic**: {info['test_statistic']}
- **Continuity correction**: {'Yes (+1 correction)' if info['continuity_correction'] else 'No'}
- **Stratification**: {info['stratification_var']}

**Variable being tested**: Effect of `{main_var}` (other variables are controlled)
""")

                        # Histogram
                        st.markdown("### 📊 Null Distribution")
                        st.caption("Gray: Distribution of the statistic under the null hypothesis. Red line: Observed statistic")

                        fig, ax = plt.subplots(figsize=(10, 6))
                        ax.hist(perm_result['null_distribution'], bins=50,
                               alpha=0.7, edgecolor='black', density=True, color='lightgray',
                               label='Null distribution')
                        ax.axvline(perm_result['observed'], color='red',
                                  linestyle='--', linewidth=2.5, label=f'Observed (t={perm_result["observed"]:.3f})')
                        ax.axvline(-perm_result['observed'], color='red',
                                  linestyle='--', linewidth=2.5)

                        # Fill the p-value region
                        null_dist = perm_result['null_distribution']
                        obs = perm_result['observed']
                        extreme_vals = null_dist[np.abs(null_dist) >= np.abs(obs)]
                        if len(extreme_vals) > 0:
                            ax.hist(extreme_vals, bins=20, alpha=0.5, color='red',
                                   density=True, label=f'p={perm_result["p_value"]:.4g}')

                        ax.set_xlabel('Test Statistic (t-value)', fontsize=12)
                        ax.set_ylabel('Density', fontsize=12)
                        ax.set_title(f'Permutation Null Distribution ({n_permutations:,} permutations)',
                                    fontsize=13, fontweight='bold')
                        ax.legend(fontsize=10)
                        ax.grid(alpha=0.3)
                        plt.tight_layout()
                        st.pyplot(fig)
                        st.session_state.fig_null_dist = fig  # Save for zip download

                        # Save results
                        st.session_state.perm_result = perm_result

                    except Exception as e:
                        st.error(f"❌ Permutation test failed")
                        st.exception(e)

            # =============================================================================
            # Download Results
            # =============================================================================

            st.markdown("---")
            section_number_download = "9️⃣" if run_welch else "8️⃣"
            st.markdown(f"## {section_number_download} Download Results")
            st.caption("You can download analysis results and plots as a ZIP file")

            # Center the ZIP download button
            col1, col2, col3 = st.columns([1, 2, 1])
            with col2:
                # Show download button only if at least one result exists
                has_results = any([
                    'ols_coef' in st.session_state,
                    'lmm_coef' in st.session_state,
                    'perm_result' in st.session_state,
                    'welch_result' in st.session_state
                ])

                if has_results:
                    # Generate filename
                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    zip_filename = f"pca_stats_{pc_col}_{timestamp}.zip"

                    # Create ZIP
                    zip_data = create_results_zip(pc_col)

                    st.download_button(
                        label="📦 Download all results (ZIP)",
                        data=zip_data,
                        file_name=zip_filename,
                        mime="application/zip",
                        help="Download statistical tables (CSV) and plots (PNG 300dpi + PDF) as a ZIP file",
                        type="primary",
                        use_container_width=True
                    )

                    # Display contents
                    st.markdown("---")
                    st.markdown("#### 📋 ZIP file contents:")

                    contents = []
                    if 'welch_result' in st.session_state:
                        contents.append("✅ Welch's t-test results (CSV)")
                    if 'ols_coef' in st.session_state:
                        contents.append("✅ OLS coefficient table (CSV)")
                    if 'ols_anova' in st.session_state:
                        contents.append("✅ ANOVA table (CSV)")
                    if 'ols_emm' in st.session_state:
                        contents.append("✅ Estimated marginal means (CSV)")
                    if 'lmm_coef' in st.session_state:
                        contents.append("✅ LMM coefficient table (CSV)")
                    if 'perm_result' in st.session_state:
                        contents.append("✅ Permutation test results (CSV)")
                    if 'fig_forest_ols' in st.session_state:
                        contents.append("✅ OLS Forest plot (PNG + PDF)")
                    if 'fig_diagnostic' in st.session_state:
                        contents.append("✅ Diagnostic plots (PNG + PDF)")
                    if 'fig_emm' in st.session_state:
                        contents.append("✅ EMM plot (PNG + PDF)")
                    if 'fig_forest_lmm' in st.session_state:
                        contents.append("✅ LMM Forest plot (PNG + PDF)")
                    if 'fig_null_dist' in st.session_state:
                        contents.append("✅ Null distribution plot (PNG + PDF)")

                    # Display in 2 columns
                    col_a, col_b = st.columns(2)
                    half = len(contents) // 2 + len(contents) % 2
                    with col_a:
                        for item in contents[:half]:
                            st.markdown(f"- {item}")
                    with col_b:
                        for item in contents[half:]:
                            st.markdown(f"- {item}")

                else:
                    st.info("💡 Results will be available for download after running the analysis")

else:
    st.info("👆 Please upload a data file to start the analysis")

    st.markdown("---")
    st.markdown("### 📋 Data format example")
    st.markdown("""
    Data should be in **TSV** or **CSV** format with the following structure:

    - **1 row = 1 sample** (pseudobulk level)
    - **PC score columns**: PC1, PC2, PC3, PC4, etc.
    - **Categorical variables**: sex, subtype, condition, etc.
    - **Optional: Donor/Subject ID** (for repeated measurements)
    - **Optional: Numeric covariates**: age, batch, depth, etc.

    **Sample data structure:**
    ```
    sample_id    sex       subtype    donor_id    age    batch    PC1      PC2      PC3      PC4
    sample_001   Female    TypeA      donor_01    45     1        -2.3     1.2      0.5      -0.8
    sample_002   Male      TypeA      donor_02    38     1         1.5    -0.8     -1.2      0.3
    sample_003   Female    TypeB      donor_03    52     1         0.2     2.1      1.8     -1.5
    sample_004   Male      TypeB      donor_04    41     2        -1.1    -1.5      0.9      1.2
    ...
    ```

    Test data can be generated from any PCA result with sample metadata.
    """)

st.markdown("---")

with st.expander("📚 Statistical Methods Details", expanded=False):
    st.markdown("""
#### **OLS (Ordinary Least Squares)**
Fixed effects model. Assumes all samples are independent.

- **Classical SE**: Standard standard errors (assumes homoscedasticity)
- **HC3 Robust SE**: Robust to heteroscedasticity (recommended for small sample sizes)
- **Cluster-robust SE**: Accounts for within-cluster correlation (e.g., correlation between samples from the same donor)
- **Welch type**: Does not assume equal variances (two-group comparison only)
- **WLS (Weighted Least Squares)**: Handles heteroscedasticity

**Use cases**: Only one sample per donor, or when between-donor variation is small

---

#### **LMM (Linear Mixed Model)**
Handles hierarchically structured data. Accounts for between-subject variation using random effects.

- **Random intercept `(1|donor)`**: Allows different intercepts for each donor
- **REML**: Unbiased variance estimation (recommended)
- **ML**: For model comparison (AIC/BIC comparison)
- **ICC (Intraclass Correlation Coefficient)**: Indicates the proportion of between-donor variation

**Use cases**: Multiple samples per donor, or technical replicates

**Note**: Variance estimation becomes unstable with fewer than 5 donors.

---

#### **Permutation Test**
A nonparametric testing method that does not rely on distributional assumptions.

**🎯 Variable being tested:** Only primary variable 1 is tested. Primary variable 2 and other covariates are treated as control variables.

- **Freedman-Lane method**: Rigorous permutation test controlling for covariates (recommended)
  - Tests the effect of primary variable 1
  - Fits null and alternative models
  - Generates pseudo-data by permuting residuals from the null model
  - Refits the alternative model to compute the test statistic

- **Simple permutation**: Simply swaps labels
  - Swaps labels of primary variable 1
  - Use when there are no covariates or few covariates

**Use cases**:
- Small sample size (<30)
- Questionable normality of residuals
- Presence of outliers

---

#### **ANOVA Type II vs Type III**
- **Type II** (recommended): Tests main effects of each variable adjusted for other variables (without interactions)
- **Type III**: Tests including all interactions (when interactions are present)

---

#### **EMM (Estimated Marginal Means)**
Predicted mean values for each group after adjusting for covariates.

- Fixes numeric covariates at their mean values
- Computes predicted values at each level of categorical variables
- Compares with confidence intervals

**Use cases**:
- Visualizing sex differences by subtype
- Group comparisons adjusted for age or batch

---

### ⚠️ Important Notes

#### **1. Perfect Collinearity**
When two variables are perfectly aligned, their effects cannot be separated.

**Example**: If a certain subtype is entirely Female, sex effect and subtype effect are not identifiable

**Solutions**:
- Remove one of the variables
- Collect more balanced data
- Limit to descriptive comparisons

#### **2. Sample Size**
- **Number of donors for LMM**: 5 or more recommended (variance estimation unstable with fewer)
- **Minimum samples per cell**: 3 or more recommended
- **Permutation test**: Particularly useful for small sample sizes (<30)

#### **3. Interaction Terms**
Recommended only when sufficient data is available. Each cell requires sufficient samples (>=5).

#### **4. Multiple Testing**
When analyzing multiple PCs, apply **Benjamini-Hochberg FDR correction**.

#### **5. Standardization of Numeric Covariates**
For numeric covariates with widely different scales (e.g., age), pre-standardization is recommended.

---

### 📖 References

1. **Robust standard errors**:
   - Long, J. S., & Ervin, L. H. (2000). Using heteroscedasticity consistent standard errors in the linear regression model. *The American Statistician*, 54(3), 217-224.

2. **Linear Mixed Models**:
   - Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). Fitting linear mixed-effects models using lme4. *Journal of Statistical Software*, 67(1), 1-48.

3. **Freedman-Lane permutation**:
   - Freedman, D., & Lane, D. (1983). A nonstochastic interpretation of reported significance levels. *Journal of Business & Economic Statistics*, 1(4), 292-298.
   - Winkler, A. M., et al. (2014). Permutation inference for the general linear model. *NeuroImage*, 92, 381-397.

4. **Estimated Marginal Means**:
   - Searle, S. R., Speed, F. M., & Milliken, G. A. (1980). Population marginal means in the linear model: An alternative to least squares means. *The American Statistician*, 34(4), 216-221.

---

**Developer**: Claude Code
**Version**: 1.0
**Last updated**: 2025-01
""")
