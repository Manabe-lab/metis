"""
Distribution Comparison Statistical Analysis Tool
Kolmogorov-Smirnov Test & Anderson-Darling Test

Usage:
streamlit run distribution_analysis.py

Required packages:
pip install streamlit plotly scipy pandas numpy
"""

import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from scipy import stats
from scipy.stats import anderson_ksamp
import io
import warnings
warnings.filterwarnings('ignore')

# Page configuration
st.set_page_config(
    page_title="Distribution Comparison Analysis Tool",
    page_icon="📊",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
<style>
    .main {
        padding: 0rem 1rem;
    }
    .stat-box {
        background: white;
        padding: 20px;
        border-radius: 10px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        margin: 10px 0;
        border-left: 4px solid #3498db;
    }
    .p-significant {
        color: #e74c3c;
        font-weight: bold;
        font-size: 1.2em;
    }
    .p-not-significant {
        color: #27ae60;
        font-weight: bold;
        font-size: 1.2em;
    }
    .metric-card {
        background: #f8f9fa;
        padding: 15px;
        border-radius: 8px;
        text-align: center;
    }
</style>
""", unsafe_allow_html=True)

# Title and subtitle
st.title("📊 Distribution Comparison Statistical Analysis Tool")
st.markdown("**Two-group distribution comparison using Kolmogorov-Smirnov Test & Anderson-Darling Test**")

# Sidebar
with st.sidebar:
    st.header("⚙️ Settings")

    st.subheader("Data Input Method")
    input_method = st.radio(
        "Select input method",
        ["Text Input", "CSV Upload", "Sample Data"]
    )

    st.subheader("Data Transformation")
    use_log = st.checkbox("Apply Log transformation", value=False)
    if use_log:
        log_base = st.selectbox("Logarithm base", ["Natural log (e)", "Common log (10)", "Binary log (2)"])

    st.subheader("Display Settings")
    show_raw_data = st.checkbox("Show raw data", value=False)
    show_qq = st.checkbox("Show Q-Q plot", value=True)
    show_violin = st.checkbox("Show violin plot", value=True)
    bin_size = st.slider("Number of histogram bins", 10, 50, 30)

    st.subheader("Statistical Settings")
    alpha = st.slider("Significance level (α)", 0.01, 0.10, 0.05, 0.01)

    # Information box
    with st.expander("ℹ️ About the tests"):
        st.markdown("""
        **Kolmogorov-Smirnov Test**
        - Evaluates the maximum difference between two cumulative distribution functions
        - Detects differences in location, scale, and shape of distributions
        - Sensitive to differences in the central part

        **Anderson-Darling Test**
        - Evaluates the entire distribution with weighting
        - Particularly sensitive to differences in the tails (both ends)
        - Generally has higher detection power than KS test
        """)

# Main content
# Data input section
st.header("1. Data Input")

col1, col2 = st.columns(2)

# Function definitions
def parse_data(text):
    """Parse text data"""
    if not text or text.strip() == "":
        return np.array([])
    try:
        # Split by comma, space, newline, or tab
        values = text.replace(',', ' ').replace('\n', ' ').replace('\t', ' ').split()
        data = np.array([float(v) for v in values if v.strip()])
        return data[data > 0]  # Only positive values
    except:
        return np.array([])

def generate_sample_data():
    """Generate sample data (log-normal distribution)"""
    np.random.seed(42)
    # Group 1: Log-normal distribution around mean 1000
    group1 = np.random.lognormal(mean=6.9, sigma=0.5, size=500)
    # Group 2: Slightly larger (around mean 1300)
    group2 = np.random.lognormal(mean=7.1, sigma=0.5, size=480)
    return group1, group2

def calculate_stats(data):
    """Calculate descriptive statistics"""
    if len(data) == 0:
        return None
    return {
        'n': len(data),
        'mean': np.mean(data),
        'median': np.median(data),
        'std': np.std(data, ddof=1),
        'cv': (np.std(data, ddof=1) / np.mean(data)) * 100,
        'min': np.min(data),
        'max': np.max(data),
        'q1': np.percentile(data, 25),
        'q3': np.percentile(data, 75),
        'iqr': np.percentile(data, 75) - np.percentile(data, 25),
        'skew': stats.skew(data),
        'kurtosis': stats.kurtosis(data)
    }

def ks_test_manual(data1, data2):
    """Implementation of Kolmogorov-Smirnov test"""
    result = stats.ks_2samp(data1, data2)
    # Effect size interpretation
    d = result.statistic
    if d < 0.15:
        effect = "Negligible"
    elif d < 0.33:
        effect = "Small"
    elif d < 0.47:
        effect = "Medium"
    else:
        effect = "Large"
    return {
        'statistic': d,
        'p_value': result.pvalue,
        'effect_size': effect
    }

def anderson_darling_test(data1, data2):
    """Anderson-Darling test"""
    try:
        result = anderson_ksamp([data1, data2])
        # Comparison with critical value
        critical_5 = 1.961  # 5% significance level
        is_significant = result.statistic > critical_5
        return {
            'statistic': result.statistic,
            'p_value': result.pvalue if hasattr(result, 'pvalue') else None,
            'critical_value': critical_5,
            'significant': is_significant
        }
    except:
        return None

def apply_log_transform(data, base="e"):
    """Log transformation"""
    if base == "Natural log (e)":
        return np.log(data)
    elif base == "Common log (10)":
        return np.log10(data)
    else:  # Binary log
        return np.log2(data)

# Data input processing
group1_data = np.array([])
group2_data = np.array([])

if input_method == "Text Input":
    with col1:
        st.subheader("Group 1")
        group1_text = st.text_area(
            "Enter data (space, comma, or newline separated)",
            height=150,
            placeholder="Example: 123.4 234.5 345.6\nor: 123.4, 234.5, 345.6"
        )
        group1_data = parse_data(group1_text)
        if len(group1_data) > 0:
            st.success(f"✅ {len(group1_data)} data points")

    with col2:
        st.subheader("Group 2")
        group2_text = st.text_area(
            "Enter data (space, comma, or newline separated)",
            height=150,
            placeholder="Example: 234.5 345.6 456.7"
        )
        group2_data = parse_data(group2_text)
        if len(group2_data) > 0:
            st.success(f"✅ {len(group2_data)} data points")

elif input_method == "CSV Upload":
    st.info("CSV file supports 2 columns (Group1, Group2) or single column format (Group and Value columns)")
    uploaded_file = st.file_uploader("Select CSV file", type=['csv'])

    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)
        st.write("Data preview:", df.head())

        # Get column names
        columns = df.columns.tolist()

        if len(columns) >= 2:
            col1_name = st.selectbox("Group 1 column", columns, index=0)
            col2_name = st.selectbox("Group 2 column", columns, index=1 if len(columns) > 1 else 0)

            group1_data = df[col1_name].dropna().values
            group2_data = df[col2_name].dropna().values

else:  # Sample Data
    if st.button("Generate sample data"):
        group1_data, group2_data = generate_sample_data()
        st.success("Sample data generated (log-normal distribution)")

# Run analysis only if data is available
if len(group1_data) > 0 and len(group2_data) > 0:

    # Apply Log transformation
    if use_log:
        group1_transformed = apply_log_transform(group1_data, log_base if 'log_base' in locals() else "e")
        group2_transformed = apply_log_transform(group2_data, log_base if 'log_base' in locals() else "e")
        st.warning(f"⚠️ Applying {log_base if 'log_base' in locals() else 'natural log'} transformation")
    else:
        group1_transformed = group1_data
        group2_transformed = group2_data

    # Descriptive statistics
    st.header("2. Descriptive Statistics")
    
    stats1 = calculate_stats(group1_data)
    stats2 = calculate_stats(group2_data)
    
    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Group 1**")
        met_cols = st.columns(3)
        met_cols[0].metric("n", f"{stats1['n']}")
        met_cols[1].metric("Mean", f"{stats1['mean']:.2f}")
        met_cols[2].metric("Median", f"{stats1['median']:.2f}")

        met_cols2 = st.columns(3)
        met_cols2[0].metric("Std Dev", f"{stats1['std']:.2f}")
        met_cols2[1].metric("CV%", f"{stats1['cv']:.1f}%")
        met_cols2[2].metric("Skewness", f"{stats1['skew']:.2f}")

    with col2:
        st.markdown("**Group 2**")
        met_cols = st.columns(3)
        met_cols[0].metric("n", f"{stats2['n']}")
        met_cols[1].metric("Mean", f"{stats2['mean']:.2f}")
        met_cols[2].metric("Median", f"{stats2['median']:.2f}")

        met_cols2 = st.columns(3)
        met_cols2[0].metric("Std Dev", f"{stats2['std']:.2f}")
        met_cols2[1].metric("CV%", f"{stats2['cv']:.1f}%")
        met_cols2[2].metric("Skewness", f"{stats2['skew']:.2f}")

    # Raw data display
    if show_raw_data:
        st.header("3. Raw Data")
        col1, col2 = st.columns(2)
        with col1:
            st.write("Group 1 (first 100)")
            st.dataframe(pd.DataFrame(group1_data[:100], columns=["Value"]))
        with col2:
            st.write("Group 2 (first 100)")
            st.dataframe(pd.DataFrame(group2_data[:100], columns=["Value"]))

    # Statistical tests
    st.header("3. Statistical Test Results")
    
    # KS test
    ks_result = ks_test_manual(group1_transformed, group2_transformed)

    # AD test
    ad_result = anderson_darling_test(group1_transformed, group2_transformed)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("### Kolmogorov-Smirnov Test")
        st.markdown('<div class="stat-box">', unsafe_allow_html=True)

        # D statistic
        st.write(f"**D Statistic**: {ks_result['statistic']:.4f}")
        st.write(f"**Effect Size**: {ks_result['effect_size']}")

        # p-value
        if ks_result['p_value'] < alpha:
            st.markdown(f"**p-value**: <span class='p-significant'>{ks_result['p_value']:.4f} ✗</span>",
                       unsafe_allow_html=True)
            st.error(f"Significant difference (α = {alpha})")
        else:
            st.markdown(f"**p-value**: <span class='p-not-significant'>{ks_result['p_value']:.4f} ✓</span>",
                       unsafe_allow_html=True)
            st.success(f"No significant difference (α = {alpha})")

        st.markdown('</div>', unsafe_allow_html=True)

        # Interpretation
        with st.expander("Interpretation"):
            if ks_result['p_value'] < alpha:
                st.write(f"The distributions of the two groups are statistically significantly different.")
                st.write(f"The maximum cumulative distribution difference (D = {ks_result['statistic']:.3f}) indicates a {ks_result['effect_size'].lower()} effect.")
            else:
                st.write(f"No statistically significant difference was found between the distributions of the two groups.")
    
    with col2:
        st.markdown("### Anderson-Darling Test")
        st.markdown('<div class="stat-box">', unsafe_allow_html=True)

        if ad_result:
            st.write(f"**A² Statistic**: {ad_result['statistic']:.4f}")
            st.write(f"**Critical Value (5%)**: {ad_result['critical_value']:.3f}")

            if ad_result['significant']:
                st.markdown(f"**Result**: <span class='p-significant'>Significant difference ✗</span>",
                           unsafe_allow_html=True)
                st.error("Distributions are significantly different")
            else:
                st.markdown(f"**Result**: <span class='p-not-significant'>No significant difference ✓</span>",
                           unsafe_allow_html=True)
                st.success("Distributions are similar")
        else:
            st.warning("Failed to run AD test")

        st.markdown('</div>', unsafe_allow_html=True)

        # Interpretation
        with st.expander("Interpretation"):
            if ad_result and ad_result['significant']:
                st.write("The A² statistic exceeds the critical value, indicating the distributions of the two groups are significantly different.")
                st.write("The Anderson-Darling test is particularly sensitive to differences in the tails of distributions.")
            else:
                st.write("The A² statistic is below the critical value, indicating no significant difference between the distributions of the two groups.")
    
    # Graph display
    st.header("4. Visualization")

    # Histogram
    st.subheader("Histogram")
    
    fig_hist = go.Figure()

    # Group 1
    fig_hist.add_trace(go.Histogram(
        x=group1_transformed,
        name='Group 1',
        nbinsx=bin_size,
        opacity=0.6,
        marker_color='#3498db'
    ))

    # Group 2
    fig_hist.add_trace(go.Histogram(
        x=group2_transformed,
        name='Group 2',
        nbinsx=bin_size,
        opacity=0.6,
        marker_color='#e74c3c'
    ))

    fig_hist.update_layout(
        barmode='overlay',
        title='Distribution Comparison' + (' (After Log Transformation)' if use_log else ''),
        xaxis_title='Value' + (' (Log Transformed)' if use_log else ''),
        yaxis_title='Frequency',
        height=400,
        showlegend=True
    )
    
    st.plotly_chart(fig_hist, use_container_width=True)

    # Cumulative distribution function
    st.subheader("Cumulative Distribution Function (CDF)")

    sorted1 = np.sort(group1_transformed)
    sorted2 = np.sort(group2_transformed)

    fig_cdf = go.Figure()

    fig_cdf.add_trace(go.Scatter(
        x=sorted1,
        y=np.arange(1, len(sorted1) + 1) / len(sorted1),
        mode='lines',
        name='Group 1',
        line=dict(color='#3498db', width=2)
    ))

    fig_cdf.add_trace(go.Scatter(
        x=sorted2,
        y=np.arange(1, len(sorted2) + 1) / len(sorted2),
        mode='lines',
        name='Group 2',
        line=dict(color='#e74c3c', width=2)
    ))

    # Visualization of KS statistic (maximum difference location)
    # Omitted for simplification

    fig_cdf.update_layout(
        title='CDF Comparison' + (' (After Log Transformation)' if use_log else ''),
        xaxis_title='Value' + (' (Log Transformed)' if use_log else ''),
        yaxis_title='Cumulative Probability',
        height=400,
        showlegend=True
    )
    
    st.plotly_chart(fig_cdf, use_container_width=True)

    # Q-Q plot
    if show_qq:
        st.subheader("Q-Q Plot")

        # Compare at common quantiles
        n_quantiles = min(len(group1_transformed), len(group2_transformed))
        quantiles = np.linspace(0, 1, min(100, n_quantiles))

        q1 = np.quantile(group1_transformed, quantiles)
        q2 = np.quantile(group2_transformed, quantiles)

        fig_qq = go.Figure()

        # Data points
        fig_qq.add_trace(go.Scatter(
            x=q1,
            y=q2,
            mode='markers',
            name='Data',
            marker=dict(size=8, color='#764ba2', opacity=0.6)
        ))

        # 45-degree line
        min_val = min(q1.min(), q2.min())
        max_val = max(q1.max(), q2.max())
        fig_qq.add_trace(go.Scatter(
            x=[min_val, max_val],
            y=[min_val, max_val],
            mode='lines',
            name='45-degree line (identical distribution)',
            line=dict(dash='dash', color='gray')
        ))

        fig_qq.update_layout(
            title='Q-Q Plot' + (' (After Log Transformation)' if use_log else ''),
            xaxis_title='Group 1 Quantiles',
            yaxis_title='Group 2 Quantiles',
            height=400
        )

        st.plotly_chart(fig_qq, use_container_width=True)
    
    # Violin plot
    if show_violin:
        st.subheader("Violin Plot")

        # Create dataframe
        df_violin = pd.DataFrame({
            'Value': np.concatenate([group1_transformed, group2_transformed]),
            'Group': ['Group 1'] * len(group1_transformed) + ['Group 2'] * len(group2_transformed)
        })

        fig_violin = px.violin(
            df_violin,
            y='Value',
            x='Group',
            color='Group',
            box=True,
            title='Distribution Shape Comparison' + (' (After Log Transformation)' if use_log else ''),
            color_discrete_map={'Group 1': '#3498db', 'Group 2': '#e74c3c'}
        )

        fig_violin.update_layout(height=400)
        st.plotly_chart(fig_violin, use_container_width=True)
    
    # Export results
    st.header("5. Export Results")

    # Summarize results in dataframe
    results_df = pd.DataFrame({
        'Test': ['Kolmogorov-Smirnov', 'Anderson-Darling'],
        'Statistic': [ks_result['statistic'], ad_result['statistic'] if ad_result else None],
        'p-value': [ks_result['p_value'], None],
        'Significance': [
            'Significant' if ks_result['p_value'] < alpha else 'Not significant',
            'Significant' if ad_result and ad_result['significant'] else 'Not significant'
        ]
    })

    # CSV download
    csv = results_df.to_csv(index=False)
    st.download_button(
        label="📥 Download results as CSV",
        data=csv,
        file_name="distribution_test_results.csv",
        mime="text/csv"
    )
    
    # Generate report
    if st.button("📄 Generate Report"):
        report = f"""
# Distribution Comparison Analysis Report

## Data Overview
- Group 1: n = {stats1['n']}, Mean = {stats1['mean']:.2f}, Median = {stats1['median']:.2f}
- Group 2: n = {stats2['n']}, Mean = {stats2['mean']:.2f}, Median = {stats2['median']:.2f}
- Transformation: {'Log transformation (' + log_base + ')' if use_log else 'None'}

## Statistical Test Results

### Kolmogorov-Smirnov Test
- D Statistic: {ks_result['statistic']:.4f}
- p-value: {ks_result['p_value']:.4f}
- Effect Size: {ks_result['effect_size']}
- Conclusion: {'Significant difference' if ks_result['p_value'] < alpha else 'No significant difference'} (α = {alpha})

### Anderson-Darling Test
- A² Statistic: {ad_result['statistic'] if ad_result else 'N/A':.4f}
- Critical Value: {ad_result['critical_value'] if ad_result else 'N/A':.3f}
- Conclusion: {'Significant difference' if ad_result and ad_result['significant'] else 'No significant difference'}

## Interpretation
{'The distributions of the two groups are statistically significantly different.' if ks_result['p_value'] < alpha else 'No statistically significant difference was found between the distributions of the two groups.'}
"""
        st.text_area("Report", report, height=400)

else:
    st.warning("⚠️ Please enter data for both groups")

    # Usage instructions
    with st.expander("📖 Usage Instructions"):
        st.markdown("""
        1. **Data Input**: Select input method in the sidebar and enter data
        2. **Transformation Options**: Apply log transformation if needed (effective for right-skewed distributions)
        3. **Review Results**: Check statistical test results and graphs
        4. **Export**: Save results as CSV or report

        **Data Format**:
        - Numbers separated by spaces, commas, or newlines
        - For CSV files, select the columns
        - Sample data available for testing
        """)

# Footer
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666;'>
    <p>Please interpret results considering both statistical significance and practical significance</p>
    <p>Created with Streamlit 📊</p>
</div>
""", unsafe_allow_html=True)