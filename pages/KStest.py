"""
分布比較統計解析ツール
Kolmogorov-Smirnov検定 & Anderson-Darling検定

使用方法:
streamlit run distribution_analysis.py

必要なパッケージ:
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

# ページ設定
st.set_page_config(
    page_title="分布比較解析ツール",
    page_icon="📊",
    layout="wide",
    initial_sidebar_state="expanded"
)

# カスタムCSS
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

# タイトルとサブタイトル
st.title("📊 分布比較統計解析ツール")
st.markdown("**Kolmogorov-Smirnov検定 & Anderson-Darling検定による2群間分布比較**")

# サイドバー
with st.sidebar:
    st.header("⚙️ 設定")
    
    st.subheader("データ入力方法")
    input_method = st.radio(
        "入力方法を選択",
        ["テキスト入力", "CSVアップロード", "サンプルデータ"]
    )
    
    st.subheader("データ変換")
    use_log = st.checkbox("Log変換を適用", value=False)
    if use_log:
        log_base = st.selectbox("対数の底", ["自然対数 (e)", "常用対数 (10)", "二進対数 (2)"])
    
    st.subheader("表示設定")
    show_raw_data = st.checkbox("生データを表示", value=False)
    show_qq = st.checkbox("Q-Qプロットを表示", value=True)
    show_violin = st.checkbox("バイオリンプロットを表示", value=True)
    bin_size = st.slider("ヒストグラムのビン数", 10, 50, 30)
    
    st.subheader("統計設定")
    alpha = st.slider("有意水準 (α)", 0.01, 0.10, 0.05, 0.01)
    
    # 情報ボックス
    with st.expander("ℹ️ 検定について"):
        st.markdown("""
        **Kolmogorov-Smirnov検定**
        - 2つの累積分布関数の最大差を評価
        - 分布の位置、スケール、形状の違いを検出
        - 中央部の違いに敏感
        
        **Anderson-Darling検定**
        - 分布全体を重み付けして評価
        - 裾部（両端）の違いに特に敏感
        - 一般的にKS検定より検出力が高い
        """)

# メインコンテンツ
# データ入力セクション
st.header("1. データ入力")

col1, col2 = st.columns(2)

# 関数定義
def parse_data(text):
    """テキストデータをパース"""
    if not text or text.strip() == "":
        return np.array([])
    try:
        # カンマ、スペース、改行、タブで分割
        values = text.replace(',', ' ').replace('\n', ' ').replace('\t', ' ').split()
        data = np.array([float(v) for v in values if v.strip()])
        return data[data > 0]  # 正の値のみ
    except:
        return np.array([])

def generate_sample_data():
    """サンプルデータ生成（対数正規分布）"""
    np.random.seed(42)
    # グループ1: 平均1000付近の対数正規分布
    group1 = np.random.lognormal(mean=6.9, sigma=0.5, size=500)
    # グループ2: やや大きめ（平均1300付近）
    group2 = np.random.lognormal(mean=7.1, sigma=0.5, size=480)
    return group1, group2

def calculate_stats(data):
    """記述統計量の計算"""
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
    """Kolmogorov-Smirnov検定の実装"""
    result = stats.ks_2samp(data1, data2)
    # 効果量の解釈
    d = result.statistic
    if d < 0.15:
        effect = "無視できる"
    elif d < 0.33:
        effect = "小"
    elif d < 0.47:
        effect = "中"
    else:
        effect = "大"
    return {
        'statistic': d,
        'p_value': result.pvalue,
        'effect_size': effect
    }

def anderson_darling_test(data1, data2):
    """Anderson-Darling検定"""
    try:
        result = anderson_ksamp([data1, data2])
        # 臨界値との比較
        critical_5 = 1.961  # 5%有意水準
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
    """対数変換"""
    if base == "自然対数 (e)":
        return np.log(data)
    elif base == "常用対数 (10)":
        return np.log10(data)
    else:  # 二進対数
        return np.log2(data)

# データ入力処理
group1_data = np.array([])
group2_data = np.array([])

if input_method == "テキスト入力":
    with col1:
        st.subheader("グループ1")
        group1_text = st.text_area(
            "データを入力（スペース、カンマ、改行区切り）",
            height=150,
            placeholder="例: 123.4 234.5 345.6\nまたは: 123.4, 234.5, 345.6"
        )
        group1_data = parse_data(group1_text)
        if len(group1_data) > 0:
            st.success(f"✅ {len(group1_data)}個のデータ")
    
    with col2:
        st.subheader("グループ2")
        group2_text = st.text_area(
            "データを入力（スペース、カンマ、改行区切り）",
            height=150,
            placeholder="例: 234.5 345.6 456.7"
        )
        group2_data = parse_data(group2_text)
        if len(group2_data) > 0:
            st.success(f"✅ {len(group2_data)}個のデータ")

elif input_method == "CSVアップロード":
    st.info("CSVファイルは2列（Group1, Group2）または1列（Group列とValue列）の形式に対応")
    uploaded_file = st.file_uploader("CSVファイルを選択", type=['csv'])
    
    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)
        st.write("データプレビュー:", df.head())
        
        # 列名を取得
        columns = df.columns.tolist()
        
        if len(columns) >= 2:
            col1_name = st.selectbox("グループ1の列", columns, index=0)
            col2_name = st.selectbox("グループ2の列", columns, index=1 if len(columns) > 1 else 0)
            
            group1_data = df[col1_name].dropna().values
            group2_data = df[col2_name].dropna().values

else:  # サンプルデータ
    if st.button("サンプルデータを生成"):
        group1_data, group2_data = generate_sample_data()
        st.success("サンプルデータを生成しました（対数正規分布）")

# データがある場合のみ解析実行
if len(group1_data) > 0 and len(group2_data) > 0:
    
    # Log変換の適用
    if use_log:
        group1_transformed = apply_log_transform(group1_data, log_base if 'log_base' in locals() else "e")
        group2_transformed = apply_log_transform(group2_data, log_base if 'log_base' in locals() else "e")
        st.warning(f"⚠️ {log_base if 'log_base' in locals() else '自然対数'}変換を適用中")
    else:
        group1_transformed = group1_data
        group2_transformed = group2_data
    
    # 記述統計
    st.header("2. 記述統計量")
    
    stats1 = calculate_stats(group1_data)
    stats2 = calculate_stats(group2_data)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("**グループ1**")
        met_cols = st.columns(3)
        met_cols[0].metric("n", f"{stats1['n']}")
        met_cols[1].metric("平均", f"{stats1['mean']:.2f}")
        met_cols[2].metric("中央値", f"{stats1['median']:.2f}")
        
        met_cols2 = st.columns(3)
        met_cols2[0].metric("標準偏差", f"{stats1['std']:.2f}")
        met_cols2[1].metric("CV%", f"{stats1['cv']:.1f}%")
        met_cols2[2].metric("歪度", f"{stats1['skew']:.2f}")
    
    with col2:
        st.markdown("**グループ2**")
        met_cols = st.columns(3)
        met_cols[0].metric("n", f"{stats2['n']}")
        met_cols[1].metric("平均", f"{stats2['mean']:.2f}")
        met_cols[2].metric("中央値", f"{stats2['median']:.2f}")
        
        met_cols2 = st.columns(3)
        met_cols2[0].metric("標準偏差", f"{stats2['std']:.2f}")
        met_cols2[1].metric("CV%", f"{stats2['cv']:.1f}%")
        met_cols2[2].metric("歪度", f"{stats2['skew']:.2f}")
    
    # 生データ表示
    if show_raw_data:
        st.header("3. 生データ")
        col1, col2 = st.columns(2)
        with col1:
            st.write("グループ1 (最初の100個)")
            st.dataframe(pd.DataFrame(group1_data[:100], columns=["Value"]))
        with col2:
            st.write("グループ2 (最初の100個)")
            st.dataframe(pd.DataFrame(group2_data[:100], columns=["Value"]))
    
    # 統計検定
    st.header("3. 統計検定結果")
    
    # KS検定
    ks_result = ks_test_manual(group1_transformed, group2_transformed)
    
    # AD検定
    ad_result = anderson_darling_test(group1_transformed, group2_transformed)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("### Kolmogorov-Smirnov検定")
        st.markdown('<div class="stat-box">', unsafe_allow_html=True)
        
        # D統計量
        st.write(f"**D統計量**: {ks_result['statistic']:.4f}")
        st.write(f"**効果量**: {ks_result['effect_size']}")
        
        # p値
        if ks_result['p_value'] < alpha:
            st.markdown(f"**p値**: <span class='p-significant'>{ks_result['p_value']:.4f} ✗</span>", 
                       unsafe_allow_html=True)
            st.error(f"有意差あり（α = {alpha}）")
        else:
            st.markdown(f"**p値**: <span class='p-not-significant'>{ks_result['p_value']:.4f} ✓</span>", 
                       unsafe_allow_html=True)
            st.success(f"有意差なし（α = {alpha}）")
        
        st.markdown('</div>', unsafe_allow_html=True)
        
        # 解釈
        with st.expander("解釈"):
            if ks_result['p_value'] < alpha:
                st.write(f"2群の分布は統計的に有意に異なります。")
                st.write(f"最大累積分布差（D = {ks_result['statistic']:.3f}）は{ks_result['effect_size']}効果を示しています。")
            else:
                st.write(f"2群の分布に統計的有意差は認められません。")
    
    with col2:
        st.markdown("### Anderson-Darling検定")
        st.markdown('<div class="stat-box">', unsafe_allow_html=True)
        
        if ad_result:
            st.write(f"**A²統計量**: {ad_result['statistic']:.4f}")
            st.write(f"**臨界値 (5%)**: {ad_result['critical_value']:.3f}")
            
            if ad_result['significant']:
                st.markdown(f"**結果**: <span class='p-significant'>有意差あり ✗</span>", 
                           unsafe_allow_html=True)
                st.error("分布は有意に異なる")
            else:
                st.markdown(f"**結果**: <span class='p-not-significant'>有意差なし ✓</span>", 
                           unsafe_allow_html=True)
                st.success("分布は類似している")
        else:
            st.warning("AD検定の実行に失敗しました")
        
        st.markdown('</div>', unsafe_allow_html=True)
        
        # 解釈
        with st.expander("解釈"):
            if ad_result and ad_result['significant']:
                st.write("A²統計量が臨界値を超えており、2群の分布は有意に異なります。")
                st.write("Anderson-Darling検定は分布の裾部の違いに特に敏感です。")
            else:
                st.write("A²統計量が臨界値以下であり、2群の分布に有意差は認められません。")
    
    # グラフ表示
    st.header("4. 可視化")
    
    # ヒストグラム
    st.subheader("ヒストグラム")
    
    fig_hist = go.Figure()
    
    # グループ1
    fig_hist.add_trace(go.Histogram(
        x=group1_transformed,
        name='グループ1',
        nbinsx=bin_size,
        opacity=0.6,
        marker_color='#3498db'
    ))
    
    # グループ2
    fig_hist.add_trace(go.Histogram(
        x=group2_transformed,
        name='グループ2',
        nbinsx=bin_size,
        opacity=0.6,
        marker_color='#e74c3c'
    ))
    
    fig_hist.update_layout(
        barmode='overlay',
        title='分布の比較' + (' (Log変換後)' if use_log else ''),
        xaxis_title='値' + (' (Log変換)' if use_log else ''),
        yaxis_title='頻度',
        height=400,
        showlegend=True
    )
    
    st.plotly_chart(fig_hist, use_container_width=True)
    
    # 累積分布関数
    st.subheader("累積分布関数 (CDF)")
    
    sorted1 = np.sort(group1_transformed)
    sorted2 = np.sort(group2_transformed)
    
    fig_cdf = go.Figure()
    
    fig_cdf.add_trace(go.Scatter(
        x=sorted1,
        y=np.arange(1, len(sorted1) + 1) / len(sorted1),
        mode='lines',
        name='グループ1',
        line=dict(color='#3498db', width=2)
    ))
    
    fig_cdf.add_trace(go.Scatter(
        x=sorted2,
        y=np.arange(1, len(sorted2) + 1) / len(sorted2),
        mode='lines',
        name='グループ2',
        line=dict(color='#e74c3c', width=2)
    ))
    
    # KS統計量の可視化（最大差の位置）
    # 簡略化のため省略
    
    fig_cdf.update_layout(
        title='累積分布関数の比較' + (' (Log変換後)' if use_log else ''),
        xaxis_title='値' + (' (Log変換)' if use_log else ''),
        yaxis_title='累積確率',
        height=400,
        showlegend=True
    )
    
    st.plotly_chart(fig_cdf, use_container_width=True)
    
    # Q-Qプロット
    if show_qq:
        st.subheader("Q-Qプロット")
        
        # 共通の分位点で比較
        n_quantiles = min(len(group1_transformed), len(group2_transformed))
        quantiles = np.linspace(0, 1, min(100, n_quantiles))
        
        q1 = np.quantile(group1_transformed, quantiles)
        q2 = np.quantile(group2_transformed, quantiles)
        
        fig_qq = go.Figure()
        
        # データ点
        fig_qq.add_trace(go.Scatter(
            x=q1,
            y=q2,
            mode='markers',
            name='データ',
            marker=dict(size=8, color='#764ba2', opacity=0.6)
        ))
        
        # 45度線
        min_val = min(q1.min(), q2.min())
        max_val = max(q1.max(), q2.max())
        fig_qq.add_trace(go.Scatter(
            x=[min_val, max_val],
            y=[min_val, max_val],
            mode='lines',
            name='45度線（同一分布）',
            line=dict(dash='dash', color='gray')
        ))
        
        fig_qq.update_layout(
            title='Q-Qプロット' + (' (Log変換後)' if use_log else ''),
            xaxis_title='グループ1の分位点',
            yaxis_title='グループ2の分位点',
            height=400
        )
        
        st.plotly_chart(fig_qq, use_container_width=True)
    
    # バイオリンプロット
    if show_violin:
        st.subheader("バイオリンプロット")
        
        # データフレーム作成
        df_violin = pd.DataFrame({
            'Value': np.concatenate([group1_transformed, group2_transformed]),
            'Group': ['グループ1'] * len(group1_transformed) + ['グループ2'] * len(group2_transformed)
        })
        
        fig_violin = px.violin(
            df_violin, 
            y='Value', 
            x='Group', 
            color='Group',
            box=True,
            title='分布の形状比較' + (' (Log変換後)' if use_log else ''),
            color_discrete_map={'グループ1': '#3498db', 'グループ2': '#e74c3c'}
        )
        
        fig_violin.update_layout(height=400)
        st.plotly_chart(fig_violin, use_container_width=True)
    
    # 結果のエクスポート
    st.header("5. 結果のエクスポート")
    
    # 結果をデータフレームにまとめる
    results_df = pd.DataFrame({
        '検定': ['Kolmogorov-Smirnov', 'Anderson-Darling'],
        '統計量': [ks_result['statistic'], ad_result['statistic'] if ad_result else None],
        'p値': [ks_result['p_value'], None],
        '有意性': [
            '有意' if ks_result['p_value'] < alpha else '非有意',
            '有意' if ad_result and ad_result['significant'] else '非有意'
        ]
    })
    
    # CSVダウンロード
    csv = results_df.to_csv(index=False)
    st.download_button(
        label="📥 結果をCSVでダウンロード",
        data=csv,
        file_name="distribution_test_results.csv",
        mime="text/csv"
    )
    
    # レポート生成
    if st.button("📄 レポートを生成"):
        report = f"""
# 分布比較解析レポート

## データ概要
- グループ1: n = {stats1['n']}, 平均 = {stats1['mean']:.2f}, 中央値 = {stats1['median']:.2f}
- グループ2: n = {stats2['n']}, 平均 = {stats2['mean']:.2f}, 中央値 = {stats2['median']:.2f}
- 変換: {'Log変換 (' + log_base + ')' if use_log else 'なし'}

## 統計検定結果

### Kolmogorov-Smirnov検定
- D統計量: {ks_result['statistic']:.4f}
- p値: {ks_result['p_value']:.4f}
- 効果量: {ks_result['effect_size']}
- 結論: {'有意差あり' if ks_result['p_value'] < alpha else '有意差なし'} (α = {alpha})

### Anderson-Darling検定
- A²統計量: {ad_result['statistic'] if ad_result else 'N/A':.4f}
- 臨界値: {ad_result['critical_value'] if ad_result else 'N/A':.3f}
- 結論: {'有意差あり' if ad_result and ad_result['significant'] else '有意差なし'}

## 解釈
{'2群の分布は統計的に有意に異なることが示されました。' if ks_result['p_value'] < alpha else '2群の分布に統計的有意差は認められませんでした。'}
"""
        st.text_area("レポート", report, height=400)

else:
    st.warning("⚠️ 両グループにデータを入力してください")
    
    # 使用方法
    with st.expander("📖 使用方法"):
        st.markdown("""
        1. **データ入力**: サイドバーで入力方法を選択し、データを入力
        2. **変換オプション**: 必要に応じてLog変換を適用（右に歪んだ分布に有効）
        3. **結果確認**: 統計検定結果とグラフを確認
        4. **エクスポート**: 結果をCSVやレポートとして保存
        
        **データ形式**:
        - スペース、カンマ、改行で区切られた数値
        - CSVファイルの場合は列を選択
        - サンプルデータで動作確認可能
        """)

# フッター
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666;'>
    <p>統計的有意性と実質的意義の両方を考慮して結果を解釈してください</p>
    <p>Created with Streamlit 📊</p>
</div>
""", unsafe_allow_html=True)