"""
divdistComparisonStatisticalAnalysistsu-ru
Kolmogorov-Smirnovcheckset & Anderson-Darlingcheckset

useusewaymethod:
streamlit run distribution_analysis.py

requirednapake-ji:
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

# pe-jiSettings
st.set_page_config(
    page_title="divdistComparisonAnalysistsu-ru",
    page_icon="📊",
    layout="wide",
    initial_sidebar_state="expanded"
)

# kasutamuCSS
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

# taitoruandsabutaitoru
st.title("📊 divdistComparisonStatisticalAnalysistsu-ru")
st.markdown("**Kolmogorov-Smirnovcheckset & Anderson-Darlingchecksettoyoru2groupbetweendivdistComparison**")

# saidoba-
with st.sidebar:
    st.header("⚙️ Settings")
    
    st.subheader("DataInputwaymethod")
    input_method = st.radio(
        "InputwaymethodtheSelect",
        ["tekisutoInput", "CSVUpload", "SampleData"]
    )
    
    st.subheader("Datachangechange")
    use_log = st.checkbox("Logchangechangethefituse", value=False)
    if use_log:
        log_base = st.selectbox("pairnumof底", ["selfnaturalpairnum (e)", "normalusepairnum (10)", "2progpairnum (2)"])
    
    st.subheader("DisplaySettings")
    show_raw_data = st.checkbox("genDatatheDisplay", value=False)
    show_qq = st.checkbox("Q-QpurototheDisplay", value=True)
    show_violin = st.checkbox("baiorinpurototheDisplay", value=True)
    bin_size = st.slider("hisutoguramuofbinnum", 10, 50, 30)
    
    st.subheader("StatisticalSettings")
    alpha = st.slider("Significant水level (α)", 0.01, 0.10, 0.05, 0.01)
    
    # infoinfobokusu
    with st.expander("ℹ️ checksettotsuite"):
        st.markdown("""
        **Kolmogorov-Smirnovcheckset**
        - 2tsuofcumaccumulatedivdistrelnumofmaximumdifftheevalval
        - divdistofrankplace、suke-ru、shapestateofdiffithecheckout
        - midcenterpartofdiffitosenssense
        
        **Anderson-Darlingcheckset**
        - divdistallbodytheweightmiattachkeshiteevalval
        - 裾part（bothedge）ofdiffitospectosenssense
        - onegenaltoKSchecksetthancheckOutputishighi
        """)

# meinkontentsu
# DataInputsekushiyon
st.header("1. DataInput")

col1, col2 = st.columns(2)

# relnumsetdef
def parse_data(text):
    """tekisutoDatathepa-su"""
    if not text or text.strip() == "":
        return np.array([])
    try:
        # kanma、supe-su、reformrow、tabuwithdivratio
        values = text.replace(',', ' ').replace('\n', ' ').replace('\t', ' ').split()
        data = np.array([float(v) for v in values if v.strip()])
        return data[data > 0]  # correctofvalofmi
    except:
        return np.array([])

def generate_sample_data():
    """SampleDatagenbecome（pairnumcorrectruledivdist）"""
    np.random.seed(42)
    # Group1: flatavg1000attachnearofpairnumcorrectruledivdist
    group1 = np.random.lognormal(mean=6.9, sigma=0.5, size=500)
    # Group2: yayabigkime（flatavg1300attachnear）
    group2 = np.random.lognormal(mean=7.1, sigma=0.5, size=480)
    return group1, group2

def calculate_stats(data):
    """recdescribeStatisticalamountofCalculation"""
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
    """Kolmogorov-Smirnovchecksetofrealinstall"""
    result = stats.ks_2samp(data1, data2)
    # effectresultamountofsolveinterp
    d = result.statistic
    if d < 0.15:
        effect = "noviewwithkiru"
    elif d < 0.33:
        effect = "small"
    elif d < 0.47:
        effect = "mid"
    else:
        effect = "big"
    return {
        'statistic': d,
        'p_value': result.pvalue,
        'effect_size': effect
    }

def anderson_darling_test(data1, data2):
    """Anderson-Darlingcheckset"""
    try:
        result = anderson_ksamp([data1, data2])
        # 临worldvalandofComparison
        critical_5 = 1.961  # 5%Significant水level
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
    """pairnumchangechange"""
    if base == "selfnaturalpairnum (e)":
        return np.log(data)
    elif base == "normalusepairnum (10)":
        return np.log10(data)
    else:  # 2progpairnum
        return np.log2(data)

# DataInputprocproc
group1_data = np.array([])
group2_data = np.array([])

if input_method == "tekisutoInput":
    with col1:
        st.subheader("Group1")
        group1_text = st.text_area(
            "DatatheInput（supe-su、kanma、reformrowareacutri）",
            height=150,
            placeholder="Example: 123.4 234.5 345.6\nalsois: 123.4, 234.5, 345.6"
        )
        group1_data = parse_data(group1_text)
        if len(group1_data) > 0:
            st.success(f"✅ {len(group1_data)}pieceofData")
    
    with col2:
        st.subheader("Group2")
        group2_text = st.text_area(
            "DatatheInput（supe-su、kanma、reformrowareacutri）",
            height=150,
            placeholder="Example: 234.5 345.6 456.7"
        )
        group2_data = parse_data(group2_text)
        if len(group2_data) > 0:
            st.success(f"✅ {len(group2_data)}pieceofData")

elif input_method == "CSVUpload":
    st.info("CSVFileis2col（Group1, Group2）alsois1col（GroupcolandValuecol）ofshapeformattopairrespond")
    uploaded_file = st.file_uploader("CSVFiletheSelect", type=['csv'])
    
    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)
        st.write("Datapurebiyu-:", df.head())
        
        # colnamethegetget
        columns = df.columns.tolist()
        
        if len(columns) >= 2:
            col1_name = st.selectbox("Group1ofcol", columns, index=0)
            col2_name = st.selectbox("Group2ofcol", columns, index=1 if len(columns) > 1 else 0)
            
            group1_data = df[col1_name].dropna().values
            group2_data = df[col2_name].dropna().values

else:  # SampleData
    if st.button("SampleDatathegenbecome"):
        group1_data, group2_data = generate_sample_data()
        st.success("SampleDatathegenbecomeshimashita（pairnumcorrectruledivdist）")

# DataisexistplacematchofmiAnalysisRun
if len(group1_data) > 0 and len(group2_data) > 0:
    
    # Logchangechangeoffituse
    if use_log:
        group1_transformed = apply_log_transform(group1_data, log_base if 'log_base' in locals() else "e")
        group2_transformed = apply_log_transform(group2_data, log_base if 'log_base' in locals() else "e")
        st.warning(f"⚠️ {log_base if 'log_base' in locals() else 'selfnaturalpairnum'}changechangethefitusemid")
    else:
        group1_transformed = group1_data
        group2_transformed = group2_data
    
    # recdescribeStatistical
    st.header("2. recdescribeStatisticalamount")
    
    stats1 = calculate_stats(group1_data)
    stats2 = calculate_stats(group2_data)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("**Group1**")
        met_cols = st.columns(3)
        met_cols[0].metric("n", f"{stats1['n']}")
        met_cols[1].metric("flatavg", f"{stats1['mean']:.2f}")
        met_cols[2].metric("midcenterval", f"{stats1['median']:.2f}")
        
        met_cols2 = st.columns(3)
        met_cols2[0].metric("marklevelbiasdiff", f"{stats1['std']:.2f}")
        met_cols2[1].metric("CV%", f"{stats1['cv']:.1f}%")
        met_cols2[2].metric("skewdegree", f"{stats1['skew']:.2f}")
    
    with col2:
        st.markdown("**Group2**")
        met_cols = st.columns(3)
        met_cols[0].metric("n", f"{stats2['n']}")
        met_cols[1].metric("flatavg", f"{stats2['mean']:.2f}")
        met_cols[2].metric("midcenterval", f"{stats2['median']:.2f}")
        
        met_cols2 = st.columns(3)
        met_cols2[0].metric("marklevelbiasdiff", f"{stats2['std']:.2f}")
        met_cols2[1].metric("CV%", f"{stats2['cv']:.1f}%")
        met_cols2[2].metric("skewdegree", f"{stats2['skew']:.2f}")
    
    # genDataDisplay
    if show_raw_data:
        st.header("3. genData")
        col1, col2 = st.columns(2)
        with col1:
            st.write("Group1 (mostfirstof100piece)")
            st.dataframe(pd.DataFrame(group1_data[:100], columns=["Value"]))
        with col2:
            st.write("Group2 (mostfirstof100piece)")
            st.dataframe(pd.DataFrame(group2_data[:100], columns=["Value"]))
    
    # Statisticalcheckset
    st.header("3. StatisticalchecksetResult")
    
    # KScheckset
    ks_result = ks_test_manual(group1_transformed, group2_transformed)
    
    # ADcheckset
    ad_result = anderson_darling_test(group1_transformed, group2_transformed)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("### Kolmogorov-Smirnovcheckset")
        st.markdown('<div class="stat-box">', unsafe_allow_html=True)
        
        # DStatisticalamount
        st.write(f"**DStatisticalamount**: {ks_result['statistic']:.4f}")
        st.write(f"**effectresultamount**: {ks_result['effect_size']}")
        
        # pval
        if ks_result['p_value'] < alpha:
            st.markdown(f"**pval**: <span class='p-significant'>{ks_result['p_value']:.4f} ✗</span>", 
                       unsafe_allow_html=True)
            st.error(f"Significantdiffari（α = {alpha}）")
        else:
            st.markdown(f"**pval**: <span class='p-not-significant'>{ks_result['p_value']:.4f} ✓</span>", 
                       unsafe_allow_html=True)
            st.success(f"Significantdiffnashi（α = {alpha}）")
        
        st.markdown('</div>', unsafe_allow_html=True)
        
        # solveinterp
        with st.expander("solveinterp"):
            if ks_result['p_value'] < alpha:
                st.write(f"2groupofdivdistisStatisticalaltoSignificanttodiffnarimasu。")
                st.write(f"maximumcumaccumulatedivdistdiff（D = {ks_result['statistic']:.3f}）is{ks_result['effect_size']}effectresulttheshowshiteimasu。")
            else:
                st.write(f"2groupofdivdisttoStatisticalalSignificantdiffisadmitmeraremasen。")
    
    with col2:
        st.markdown("### Anderson-Darlingcheckset")
        st.markdown('<div class="stat-box">', unsafe_allow_html=True)
        
        if ad_result:
            st.write(f"**A²Statisticalamount**: {ad_result['statistic']:.4f}")
            st.write(f"**临worldval (5%)**: {ad_result['critical_value']:.3f}")
            
            if ad_result['significant']:
                st.markdown(f"**Result**: <span class='p-significant'>Significantdiffari ✗</span>", 
                           unsafe_allow_html=True)
                st.error("divdistisSignificanttodiffbecome")
            else:
                st.markdown(f"**Result**: <span class='p-not-significant'>Significantdiffnashi ✓</span>", 
                           unsafe_allow_html=True)
                st.success("divdististypelikeshiteexist")
        else:
            st.warning("ADchecksetofRuntofailfailshimashita")
        
        st.markdown('</div>', unsafe_allow_html=True)
        
        # solveinterp
        with st.expander("solveinterp"):
            if ad_result and ad_result['significant']:
                st.write("A²Statisticalamountis临worldvalthesupereteori、2groupofdivdistisSignificanttodiffnarimasu。")
                st.write("Anderson-Darlingchecksetisdivdistof裾partofdiffitospectosenssensewithsu。")
            else:
                st.write("A²Statisticalamountis临worldvalorbelowwithari、2groupofdivdisttoSignificantdiffisadmitmeraremasen。")
    
    # gurafuDisplay
    st.header("4. Visualization")
    
    # hisutoguramu
    st.subheader("hisutoguramu")
    
    fig_hist = go.Figure()
    
    # Group1
    fig_hist.add_trace(go.Histogram(
        x=group1_transformed,
        name='Group1',
        nbinsx=bin_size,
        opacity=0.6,
        marker_color='#3498db'
    ))
    
    # Group2
    fig_hist.add_trace(go.Histogram(
        x=group2_transformed,
        name='Group2',
        nbinsx=bin_size,
        opacity=0.6,
        marker_color='#e74c3c'
    ))
    
    fig_hist.update_layout(
        barmode='overlay',
        title='divdistofComparison' + (' (Logchangechangeafter)' if use_log else ''),
        xaxis_title='val' + (' (Logchangechange)' if use_log else ''),
        yaxis_title='freqdegree',
        height=400,
        showlegend=True
    )
    
    st.plotly_chart(fig_hist, use_container_width=True)
    
    # cumaccumulatedivdistrelnum
    st.subheader("cumaccumulatedivdistrelnum (CDF)")
    
    sorted1 = np.sort(group1_transformed)
    sorted2 = np.sort(group2_transformed)
    
    fig_cdf = go.Figure()
    
    fig_cdf.add_trace(go.Scatter(
        x=sorted1,
        y=np.arange(1, len(sorted1) + 1) / len(sorted1),
        mode='lines',
        name='Group1',
        line=dict(color='#3498db', width=2)
    ))
    
    fig_cdf.add_trace(go.Scatter(
        x=sorted2,
        y=np.arange(1, len(sorted2) + 1) / len(sorted2),
        mode='lines',
        name='Group2',
        line=dict(color='#e74c3c', width=2)
    ))
    
    # KSStatisticalamountofVisualization（maximumdiffofrankplace）
    # simpleabbrevizeoffor省abbrev
    
    fig_cdf.update_layout(
        title='cumaccumulatedivdistrelnumofComparison' + (' (Logchangechangeafter)' if use_log else ''),
        xaxis_title='val' + (' (Logchangechange)' if use_log else ''),
        yaxis_title='cumaccumulatecertainrate',
        height=400,
        showlegend=True
    )
    
    st.plotly_chart(fig_cdf, use_container_width=True)
    
    # Q-Qpuroto
    if show_qq:
        st.subheader("Q-Qpuroto")
        
        # copassofdivrankpointwithComparison
        n_quantiles = min(len(group1_transformed), len(group2_transformed))
        quantiles = np.linspace(0, 1, min(100, n_quantiles))
        
        q1 = np.quantile(group1_transformed, quantiles)
        q2 = np.quantile(group2_transformed, quantiles)
        
        fig_qq = go.Figure()
        
        # Datapoint
        fig_qq.add_trace(go.Scatter(
            x=q1,
            y=q2,
            mode='markers',
            name='Data',
            marker=dict(size=8, color='#764ba2', opacity=0.6)
        ))
        
        # 45degreeline
        min_val = min(q1.min(), q2.min())
        max_val = max(q1.max(), q2.max())
        fig_qq.add_trace(go.Scatter(
            x=[min_val, max_val],
            y=[min_val, max_val],
            mode='lines',
            name='45degreeline（sameonedivdist）',
            line=dict(dash='dash', color='gray')
        ))
        
        fig_qq.update_layout(
            title='Q-Qpuroto' + (' (Logchangechangeafter)' if use_log else ''),
            xaxis_title='Group1ofdivrankpoint',
            yaxis_title='Group2ofdivrankpoint',
            height=400
        )
        
        st.plotly_chart(fig_qq, use_container_width=True)
    
    # baiorinpuroto
    if show_violin:
        st.subheader("baiorinpuroto")
        
        # Datafure-mumakebecome
        df_violin = pd.DataFrame({
            'Value': np.concatenate([group1_transformed, group2_transformed]),
            'Group': ['Group1'] * len(group1_transformed) + ['Group2'] * len(group2_transformed)
        })
        
        fig_violin = px.violin(
            df_violin, 
            y='Value', 
            x='Group', 
            color='Group',
            box=True,
            title='divdistofshapestateComparison' + (' (Logchangechangeafter)' if use_log else ''),
            color_discrete_map={'Group1': '#3498db', 'Group2': '#e74c3c'}
        )
        
        fig_violin.update_layout(height=400)
        st.plotly_chart(fig_violin, use_container_width=True)
    
    # Resultofekusupo-to
    st.header("5. Resultofekusupo-to")
    
    # ResulttheDatafure-mutomaandmeru
    results_df = pd.DataFrame({
        'checkset': ['Kolmogorov-Smirnov', 'Anderson-Darling'],
        'Statisticalamount': [ks_result['statistic'], ad_result['statistic'] if ad_result else None],
        'pval': [ks_result['p_value'], None],
        'Significantnature': [
            'Significant' if ks_result['p_value'] < alpha else 'nonSignificant',
            'Significant' if ad_result and ad_result['significant'] else 'nonSignificant'
        ]
    })
    
    # CSVDownload
    csv = results_df.to_csv(index=False)
    st.download_button(
        label="📥 ResulttheCSVwithDownload",
        data=csv,
        file_name="distribution_test_results.csv",
        mime="text/csv"
    )
    
    # repo-togenbecome
    if st.button("📄 repo-tothegenbecome"):
        report = f"""
# divdistComparisonAnalysisrepo-to

## Dataoverviewneed
- Group1: n = {stats1['n']}, flatavg = {stats1['mean']:.2f}, midcenterval = {stats1['median']:.2f}
- Group2: n = {stats2['n']}, flatavg = {stats2['mean']:.2f}, midcenterval = {stats2['median']:.2f}
- changechange: {'Logchangechange (' + log_base + ')' if use_log else 'nashi'}

## StatisticalchecksetResult

### Kolmogorov-Smirnovcheckset
- DStatisticalamount: {ks_result['statistic']:.4f}
- pval: {ks_result['p_value']:.4f}
- effectresultamount: {ks_result['effect_size']}
- resultlogic: {'Significantdiffari' if ks_result['p_value'] < alpha else 'Significantdiffnashi'} (α = {alpha})

### Anderson-Darlingcheckset
- A²Statisticalamount: {ad_result['statistic'] if ad_result else 'N/A':.4f}
- 临worldval: {ad_result['critical_value'] if ad_result else 'N/A':.3f}
- resultlogic: {'Significantdiffari' if ad_result and ad_result['significant'] else 'Significantdiffnashi'}

## solveinterp
{'2groupofdivdistisStatisticalaltoSignificanttodiffbecomekoandisshowsaremashita。' if ks_result['p_value'] < alpha else '2groupofdivdisttoStatisticalalSignificantdiffisadmitmeraremasenwithshita。'}
"""
        st.text_area("repo-to", report, height=400)

else:
    st.warning("⚠️ bothGrouptoDatatheInputshitekudasai")
    
    # useusewaymethod
    with st.expander("📖 useusewaymethod"):
        st.markdown("""
        1. **DataInput**: saidoba-withInputwaymethodtheSelectshi、DatatheInput
        2. **changechangeOption**: requiredtorespondjiteLogchangechangethefituse（righttoskewndadivdisttovalid）
        3. **ResultConfirm**: StatisticalchecksetResultandgurafutheConfirm
        4. **ekusupo-to**: ResulttheCSVyarepo-toandshiteSave
        
        **Datashapeformat**:
        - supe-su、kanma、reformrowwithareacutraretanumval
        - CSVFileofplacematchiscoltheSelect
        - SampleDatawithmovemakeConfirmpossible
        """)

# futa-
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666;'>
    <p>StatisticalalSignificantnatureandrealqualityalmeandefofbothwaythethinkconsidershiteResultthesolveinterpshitekudasai</p>
    <p>Created with Streamlit 📊</p>
</div>
""", unsafe_allow_html=True)