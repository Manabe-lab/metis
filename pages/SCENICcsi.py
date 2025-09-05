import streamlit as st
import pandas as pd
import numpy as np
import tempfile
import os
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
from concurrent.futures import ThreadPoolExecutor
import threading
import zipfile
import io

def install_r_packages():
    """必要なRパッケージをインストールする"""
    utils = importr('utils')
    required_packages = ['tidyverse', 'pheatmap', 'viridis']
    
    # CRANミラーの設定
    utils.chooseCRANmirror(ind=1)
    
    for package in required_packages:
        if not rpackages.isinstalled(package):
            st.write(f"Installing {package}...")
            try:
                utils.install_packages(package)
                st.success(f"Successfully installed {package}")
            except Exception as e:
                st.error(f"Failed to install {package}: {str(e)}")
                return False
    return True

def setup_r_environment():
    """R環境のセットアップと診断"""
    try:
        if not install_r_packages():
            return False
        
        for package in ['tidyverse', 'pheatmap', 'viridis']:
            ro.r(f'library({package})')
        
        return True
    except Exception as e:
        st.error(f"Error setting up R environment: {str(e)}")
        return False

def compare_pcc(vector_of_pcc, pcc):
    """PCC比較関数"""
    pcc_larger = np.sum(vector_of_pcc > pcc)
    return 0 if pcc_larger == len(vector_of_pcc) else len(vector_of_pcc)

def calc_csi(reg, reg2, pearson_cor):
    """CSI計算関数"""
    test_cor = pearson_cor[reg, reg2]
    total_n = pearson_cor.shape[1]
    pearson_cor_sub = pearson_cor[[reg, reg2], :]
    
    sums = np.apply_along_axis(compare_pcc, 0, pearson_cor_sub, pcc=test_cor)
    return np.sum(sums == pearson_cor_sub.shape[0]) / total_n

def calculate_csi_batch(regulon_pairs, pearson_cor, regulon_names):
    """バッチ処理でCSIを計算"""
    results = []
    for reg, reg2 in regulon_pairs:
        fraction_lower = calc_csi(reg, reg2, pearson_cor)
        results.append((regulon_names[reg], regulon_names[reg2], fraction_lower))
    return results

@st.cache_data(show_spinner=False, hash_funcs={pd.DataFrame: lambda x: hash(str(x.values.tobytes()))})
def calculate_csi(_regulonAUC, calc_extended=False, verbose=False):
    """CSI計算のメイン関数"""
    regulonAUC_sub = _regulonAUC.T
    if calc_extended:
        regulonAUC_sub = regulonAUC_sub.loc[:, regulonAUC_sub.columns.str.contains("extended")]
    else:
        regulonAUC_sub = regulonAUC_sub.loc[:, ~regulonAUC_sub.columns.str.contains("extended")]
    
    pearson_cor = regulonAUC_sub.corr().values
    regulon_names = regulonAUC_sub.columns.values
    n = len(regulon_names)
    
    # 全ての組み合わせを生成
    regulon_pairs = [(i, j) for i in range(n) for j in range(n)]
    
    # バッチサイズの計算（メモリ使用量を考慮）
    batch_size = min(1000, len(regulon_pairs))
    batches = [regulon_pairs[i:i + batch_size] for i in range(0, len(regulon_pairs), batch_size)]
    
    csi_data = []
    progress_text = st.empty()
    progress_bar = st.progress(0)
    
    for i, batch in enumerate(batches):
        if verbose:
            progress = (i + 1) / len(batches)
            progress_bar.progress(progress)
            progress_text.text(f"Processing batch {i+1}/{len(batches)}")
        
        results = calculate_csi_batch(batch, pearson_cor, regulon_names)
        csi_data.extend(results)
    
    progress_bar.empty()
    progress_text.empty()
    
    csi_regulons = pd.DataFrame(csi_data, columns=["regulon_1", "regulon_2", "CSI"])
    return csi_regulons.astype({'CSI': 'float'})

def plot_csi_modules_r(csi_df, nclust=10, font_size_regulons=6, plot_width=2400, plot_height=1800, file_format="png"):
    """CSIモジュールのプロット関数"""
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_csv = os.path.join(temp_dir, 'csi_data.csv')
            temp_plot = os.path.join(temp_dir, f'csi_plot.{file_format}')
            
            # CSVファイルの保存
            csi_df.to_csv(temp_csv, index=False)
            
            # Rパッケージのロード
            ro.r('library(tidyverse)')
            ro.r('library(pheatmap)')
            ro.r('library(viridis)')
            
            # データの読み込みと処理
            ro.r(f'csi_df <- read.csv("{temp_csv}")')
            ro.r('csi_test_mat <- csi_df %>% spread(regulon_2, CSI)')
            ro.r('future_rownames <- csi_test_mat$regulon_1')
            ro.r('csi_test_mat <- as.matrix(csi_test_mat[,2:ncol(csi_test_mat)])')
            ro.r('rownames(csi_test_mat) <- future_rownames')
            
            # プロットデバイスの設定
            if file_format == "png":
                ro.r(f'png("{temp_plot}", width={plot_width}, height={plot_height}, units="px", res=100)')
            else:  # pdf
                # PDFの場合、サイズをインチ単位に変換
                width_inch = plot_width / 100
                height_inch = plot_height / 100
                ro.r(f'pdf("{temp_plot}", width={width_inch}, height={height_inch})')
            
            # ヒートマップの生成
            ro.r(f'''
            pheatmap(csi_test_mat,
                    show_colnames = FALSE,
                    color = viridis(n = 10),
                    cutree_cols = {nclust},
                    cutree_rows = {nclust},
                    fontsize_row = {font_size_regulons},
                    cluster_cols = TRUE,
                    cluster_rows = TRUE,
                    treeheight_row = 20,
                    treeheight_col = 20,
                    clustering_distance_rows = "euclidean",
                    clustering_distance_cols = "euclidean")
            ''')
            
            # プロットデバイスを閉じる
            ro.r('dev.off()')
            
            # 生成されたプロットを読み込む
            with open(temp_plot, 'rb') as f:
                return f.read()
                
    except Exception as e:
        st.error(f"Error in plot generation: {str(e)}")
        return None

def create_zip_file(files):
    """Zipファイル作成関数"""
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, 'w') as zipf:
        for file_name, file_content in files:
            zipf.writestr(file_name, file_content)
    return zip_buffer.getvalue()


def main():
    """メイン関数"""
    st.title('Calculate connection specificity index (CSI) for all regulons')

    st.markdown("""CSIは2つのレギュロン間の活性パターンの類似性がどれだけ特異的かを示す

CSI値は0から1の範囲で表され：

　高いCSI値（1に近い）：2つのレギュロンが非常に似た活性パターンを特異的に共有していることを示す

　低いCSI値：特異的な関係性が弱いことを示す

単純な相関ではなく、他のレギュロンとの関係と比較して特異的に類似したパターンを持つレギュロンのペアを識別

高いCSIを共有するレギュロン：

　下流の遺伝子を共同で制御している可能性が高い 特定の細胞機能に一緒に関与している可能性
""")
 
     # キャッシュコントロール
    with st.expander("Cache Control Options"):
        if st.button('Clear All Cache', help="Click to clear all cached calculations"):
            st.cache_data.clear()
            st.success("Cache cleared successfully!")
            st.rerun()
    st.write("You may need to clear cache to calculate CIS for a new dataset.")
    st.markdown("---")

    if not setup_r_environment():
        st.error("Failed to set up R environment. Please check the R installation and required packages.")
        return
    
    uploaded_file = st.file_uploader("Choose AUC_per_cell.txt", type="txt", key='auc_file')
    
    if uploaded_file is not None:
        try:
            # ファイルの内容を読み込む
            file_content = uploaded_file.read()
            
            # ファイルの内容をデータフレームに変換
            try:
                df = pd.read_csv(io.StringIO(file_content.decode('utf-8')), index_col=0, sep='\t')
            except UnicodeDecodeError:
                df = pd.read_csv(io.StringIO(file_content.decode('latin1')), index_col=0, sep='\t')
                
            st.success("Data loaded successfully!")
            
            with st.form("plot_parameters"):
                st.header("Plot Parameters")
                cols = st.columns(2)
                with cols[0]:
                    plot_width = st.number_input("Plot Width", min_value=1000, max_value=5000, value=2400, step=100)
                    font_size = st.number_input("Font Size", min_value=4, max_value=20, value=6, step=1)
                    file_format = st.selectbox("File Format", ["png", "pdf"])
                with cols[1]:
                    plot_height = st.number_input("Plot Height", min_value=800, max_value=4000, value=2350, step=100)
                    n_clusters = st.number_input("Number of Clusters", min_value=5, max_value=20, value=10, step=1)
                
                submit_button = st.form_submit_button(label='Run CSI Analysis')
            
            if submit_button:
                with st.spinner('Calculating CSI...'):
                    # ファイルの内容に基づいてCSIを計算
                    result = calculate_csi(df, calc_extended=False, verbose=True)
                    
                    plot_content = plot_csi_modules_r(
                        result,
                        nclust=n_clusters,
                        font_size_regulons=font_size,
                        plot_width=plot_width,
                        plot_height=plot_height,
                        file_format=file_format
                    )
                
                if plot_content is not None:
                    st.success('Analysis completed successfully!')
                    
                    if file_format == "png":
                        st.image(plot_content, caption='CSI Heatmap', use_container_width=True)
                    else:
                        st.info("PDF file generated. You can download it using the button below.")
                    
                    csv_buffer = io.StringIO()
                    result.to_csv(csv_buffer, index=False)
                    
                    zip_content = create_zip_file([
                        ("csi_results.csv", csv_buffer.getvalue()),
                        (f"csi_heatmap.{file_format}", plot_content)
                    ])
                    
                    st.download_button(
                        label="Download Results and Plot",
                        data=zip_content,
                        file_name="csi_results.zip",
                        mime="application/zip"
                    )
                
        except Exception as e:
            st.error(f"Error processing file: {str(e)}")
            st.error(f"Error details: {str(e.__class__.__name__)}")
            import traceback
            st.error(f"Traceback: {traceback.format_exc()}")
    else:
        st.write("Please upload a TSV file to begin the analysis.")

if __name__ == "__main__":
    main()