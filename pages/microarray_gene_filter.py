#!/usr/bin/env python3
"""
Microarray Gene Name Filter
マイクロアレイデータから遺伝子名を抽出・フィルタリング・集約するツール
"""

import streamlit as st
import pandas as pd
import numpy as np
import io
import re
from collections import defaultdict, Counter

class GeneInfoExtractor:
    """遺伝子情報の位置パターンを学習して抽出するクラス"""

    def __init__(self):
        self.position_patterns = {}
        self.separator_pattern = None
        self.learned = False

    def detect_separator(self, sample_texts):
        """サンプルテキストから使用されている区切り文字を自動検出"""
        separators = {
            ' // ': r'\s*//\s*',
            '//': r'//',
            ' / ': r'\s*/\s*',
            '/': r'/',
            '\t': r'\t',
            ', ': r',\s*',
            ',': r',',
            ' ': r'\s+',
            '|': r'\|',
            ';': r';'
        }

        # 各区切り文字の出現回数をカウント
        counts = {}
        for name, pattern in separators.items():
            total_count = 0
            for text in sample_texts:
                if pd.notna(text) and text != "---":
                    count = len(re.findall(pattern, str(text)))
                    total_count += count
            if total_count > 0:
                counts[name] = total_count

        # 最も頻出する区切り文字を返す
        if counts:
            most_common = max(counts, key=counts.get)
            return most_common, separators[most_common]
        else:
            return ' // ', r'\s*//\s*'  # デフォルト

    def classify_element(self, element):
        """要素の種類を分類"""
        element = str(element).strip()

        if not element or element == "---":
            return "empty"

        # RefSeq ID
        if any(element.startswith(prefix) for prefix in ["NM_", "NR_", "XM_", "XR_", "NP_", "XP_"]):
            return "refseq"

        # Ensembl ID
        if element.startswith("ENSMUST") or element.startswith("ENSMUSG") or element.startswith("ENSMUS") or \
           element.startswith("ENST") or element.startswith("ENSG") or element.startswith("ENS"):
            return "ensembl"

        # GenBank clone ID
        if element.startswith("BC") and len(element) > 2 and element[2:].replace(".", "").isdigit():
            return "genbank"

        if element.startswith("AK") and len(element) > 2 and element[2:].isdigit():
            return "genbank"

        # Entrez Gene ID (数字のみ)
        if element.isdigit():
            return "entrez"

        # 染色体位置情報
        if "|" in element or "cM" in element or element.startswith("chr"):
            return "location"

        # 説明文の判定
        words = element.split()
        if len(words) >= 3:
            # 長い説明文
            if not re.match(r'^[A-Za-z0-9]+[-]?[A-Za-z0-9]+$', element):
                return "description"

        # 特定のキーワードを含む説明文
        skip_keywords = ["uncharacterized", "predicted", "hypothetical", "family", "domain",
                        "containing", "protein", "member", "homolog", "like", "antigen", "receptor"]
        if any(keyword in element.lower() for keyword in skip_keywords):
            if len(element.split()) > 2 or len(element) > 30:
                return "description"

        # 遺伝子シンボルの可能性が高い
        if re.match(r'^([A-Za-z][A-Za-z0-9_-]*\d*[A-Za-z0-9]*|[0-9]+[A-Za-z][A-Za-z0-9]*[A-Za-z]+)$', element):
            return "gene_symbol"

        # その他（おそらく説明文）
        return "description"

    def learn_pattern(self, sample_data, separator=None):
        """サンプルデータから位置パターンを学習"""
        # 区切り文字を検出
        if separator is None:
            sample_texts = [str(text) for text in sample_data if pd.notna(text) and text != "---"]
            sep_name, self.separator_pattern = self.detect_separator(sample_texts[:20])
            st.info(f"🔍 検出された区切り文字: '{sep_name}'")
        else:
            self.separator_pattern = separator

        # 各位置での要素タイプの出現頻度をカウント
        position_counts = defaultdict(Counter)
        max_positions = 0

        for text in sample_data:
            if pd.isna(text) or text == "---" or text == "":
                continue

            text = str(text)
            # 複数エントリー（///区切り）の場合は最初のエントリーのみ使用
            if " /// " in text:
                text = text.split(" /// ")[0]

            # 区切り文字で分割
            parts = re.split(self.separator_pattern, text)
            max_positions = max(max_positions, len(parts))

            for pos, part in enumerate(parts):
                element_type = self.classify_element(part)
                if element_type != "empty":
                    position_counts[pos][element_type] += 1

        # 各位置で最も頻出する要素タイプを決定
        self.position_patterns = {}

        st.write("📊 位置パターン分析結果:")
        pattern_summary = []

        for pos in range(max_positions):
            if pos in position_counts:
                # 最も頻出するタイプを取得
                most_common = position_counts[pos].most_common()
                if most_common:
                    top_type = most_common[0][0]
                    top_count = most_common[0][1]
                    total = sum(position_counts[pos].values())
                    confidence = top_count / total * 100

                    # 信頼度が40%以上の場合のみパターンとして採用
                    if confidence >= 40:
                        self.position_patterns[pos] = top_type
                        pattern_summary.append({
                            "位置": pos + 1,
                            "推定タイプ": self.get_type_label(top_type),
                            "信頼度": f"{confidence:.1f}%",
                            "サンプル数": total
                        })

        if pattern_summary:
            pattern_df = pd.DataFrame(pattern_summary)
            st.dataframe(pattern_df, use_container_width=True)

        self.learned = True
        return self.position_patterns

    def get_type_label(self, type_key):
        """タイプキーを日本語ラベルに変換"""
        labels = {
            "gene_symbol": "Gene Symbol",
            "refseq": "RefSeq ID",
            "ensembl": "Ensembl ID",
            "genbank": "GenBank ID",
            "entrez": "Entrez ID",
            "description": "Description",
            "location": "Location",
            "empty": "Empty"
        }
        return labels.get(type_key, type_key)

    def extract_with_pattern(self, text):
        """学習したパターンを使用して情報を抽出"""
        result = {
            'gene_symbol': '',
            'refseq_id': '',
            'ensembl_id': '',
            'genbank_id': '',
            'entrez_id': '',
            'description': '',
            'location': ''
        }

        if pd.isna(text) or text == "---" or text == "":
            return result

        text = str(text)

        # 複数エントリー（///区切り）の場合は最初のエントリーのみ使用
        if " /// " in text:
            text = text.split(" /// ")[0]

        # 区切り文字で分割
        parts = re.split(self.separator_pattern, text)

        # 各タイプごとに値を収集
        symbols = []
        refseqs = []
        ensembls = []
        genbanks = []
        entrezs = []
        descriptions = []
        locations = []

        for pos, part in enumerate(parts):
            part = part.strip()
            if not part or part == "---":
                continue

            # 学習したパターンがある場合は優先
            if pos in self.position_patterns:
                expected_type = self.position_patterns[pos]
                actual_type = self.classify_element(part)

                # 期待されるタイプと実際のタイプが一致する場合
                if expected_type == actual_type:
                    if expected_type == "gene_symbol" and not symbols:
                        symbols.append(part)
                    elif expected_type == "refseq":
                        refseqs.append(part)
                    elif expected_type == "ensembl":
                        ensembls.append(part)
                    elif expected_type == "genbank":
                        genbanks.append(part)
                    elif expected_type == "entrez":
                        entrezs.append(part)
                    elif expected_type == "description":
                        descriptions.append(part)
                    elif expected_type == "location":
                        locations.append(part)
                # タイプが一致しない場合でも、実際のタイプに基づいて保存
                else:
                    if actual_type == "gene_symbol" and not symbols:
                        symbols.append(part)
                    elif actual_type == "refseq":
                        refseqs.append(part)
                    elif actual_type == "ensembl":
                        ensembls.append(part)
                    elif actual_type == "genbank":
                        genbanks.append(part)
                    elif actual_type == "entrez":
                        entrezs.append(part)
                    elif actual_type == "description":
                        descriptions.append(part)
                    elif actual_type == "location":
                        locations.append(part)
            else:
                # パターンがない位置は実際のタイプで分類
                actual_type = self.classify_element(part)
                if actual_type == "gene_symbol" and not symbols:
                    symbols.append(part)
                elif actual_type == "refseq":
                    refseqs.append(part)
                elif actual_type == "ensembl":
                    ensembls.append(part)
                elif actual_type == "genbank":
                    genbanks.append(part)
                elif actual_type == "entrez":
                    entrezs.append(part)
                elif actual_type == "description":
                    descriptions.append(part)
                elif actual_type == "location":
                    locations.append(part)

        # 結果をまとめる
        if symbols:
            result['gene_symbol'] = symbols[0]
        if refseqs:
            result['refseq_id'] = "; ".join(refseqs)
        if ensembls:
            result['ensembl_id'] = "; ".join(ensembls)
        if genbanks:
            result['genbank_id'] = "; ".join(genbanks)
        if entrezs:
            result['entrez_id'] = "; ".join(entrezs)
        if descriptions:
            result['description'] = " | ".join(descriptions)
        if locations:
            result['location'] = "; ".join(locations)

        return result

@st.cache_data
def convert_df(df):
    """DataFrameをTSV形式に変換"""
    return df.to_csv(index=True, sep='\t').encode('utf-8')

@st.cache_data
def load_and_parse_file(file_content, file_name):
    """ファイルを読み込んでパースする（キャッシュ付き）"""
    delimiter = '\t' if '\t' in file_content.split('\n')[0] else ','
    df = pd.read_csv(
        io.StringIO(file_content),
        sep=delimiter,
        engine='python'
    )
    return df

@st.cache_data
def extract_gene_info(df_dict, gene_column, sample_size):
    """遺伝子情報を抽出する（キャッシュ付き）"""
    # dictからDataFrameを復元
    df = pd.DataFrame(df_dict)

    extractor = GeneInfoExtractor()

    # サンプルデータから学習
    sample_size_actual = min(sample_size, len(df))
    sample_data = df[gene_column].head(sample_size_actual).dropna()
    extractor.learn_pattern(sample_data)

    if not extractor.learned:
        return None, None

    # バッチ処理で効率化
    batch_size = 1000
    results = []

    for i in range(0, len(df), batch_size):
        batch = df[gene_column].iloc[i:i+batch_size]
        batch_results = batch.apply(extractor.extract_with_pattern)
        results.extend(batch_results.tolist())

    # 抽出した情報を個別の列に展開
    info_df = pd.DataFrame(results)

    # 空の列を削除
    cols_to_keep = []
    for col in info_df.columns:
        if info_df[col].str.len().gt(0).any():
            cols_to_keep.append(col)

    return info_df, cols_to_keep

def main():
    st.set_page_config(page_title="Microarray Gene Name Filter", page_icon="🧬", layout="wide")
    st.title("🧬 Microarray Gene Name Filter")
    st.markdown("---")

    st.markdown("""
    ### 機能
    - マイクロアレイデータから遺伝子情報を抽出
    - 遺伝子名カラムを選択して1列目に配置
    - 遺伝子名がNoneの行を削除可能
    - 重複遺伝子の検出と集約（Mean/Max）
    """)

    # ファイルアップロード
    uploaded_file = st.file_uploader(
        "TSV/CSVファイルを選択してください",
        type=['tsv', 'txt', 'csv'],
        help="マイクロアレイデータファイル"
    )

    if uploaded_file is not None:
        try:
            # ファイル読み込み（キャッシュ付き）
            content = uploaded_file.getvalue().decode('utf-8')
            df = load_and_parse_file(content, uploaded_file.name)

            st.success(f"✅ ファイル読み込み成功: {len(df):,}行 × {len(df.columns)}列")

            # データプレビュー
            st.subheader("📊 データプレビュー")

            preview_rows = st.slider(
                "プレビュー行数",
                min_value=5,
                max_value=min(100, len(df)),
                value=10
            )

            st.dataframe(df.head(preview_rows))

            # 遺伝子列の選択
            st.subheader("🔍 遺伝子情報列の選択")

            # デフォルト列を探す
            default_col = None
            for col in df.columns:
                if 'gene' in col.lower() or 'symbol' in col.lower() or 'description' in col.lower():
                    default_col = col
                    break

            gene_column = st.selectbox(
                "遺伝子情報が含まれる列を選択してください",
                options=df.columns.tolist(),
                index=df.columns.tolist().index(default_col) if default_col else 0,
                help="通常は'gene_assignment'列に遺伝子情報が含まれています"
            )

            # 選択した列のサンプル表示
            if gene_column:
                st.info("選択した列のサンプル（最初の5行）:")
                sample_data = df[gene_column].head(5).to_frame()
                st.dataframe(sample_data)

                # 学習サンプルサイズの設定
                st.subheader("⚙️ 学習設定")

                col1, col2 = st.columns(2)
                with col1:
                    sample_size = st.number_input(
                        "学習に使用するサンプル数",
                        min_value=20,
                        max_value=min(500, len(df)),
                        value=min(100, len(df)),
                        help="多いほど正確ですが処理時間が増えます"
                    )

                with col2:
                    st.info(f"""
                    📌 推奨設定:
                    - 小規模データ（<1000行）: 50-100
                    - 中規模データ（1000-10000行）: 100-200
                    - 大規模データ（>10000行）: 200-500
                    """)

                # 処理実行ボタン
                if st.button("🚀 遺伝子情報を抽出", type="primary"):
                    st.subheader("📈 位置パターンの学習と抽出")

                    with st.spinner("遺伝子情報を抽出中..."):
                        # DataFrameをdictに変換してキャッシュ可能にする
                        df_dict = df.to_dict('list')
                        info_df, cols_to_keep = extract_gene_info(df_dict, gene_column, sample_size)

                    if info_df is None:
                        st.error("パターンの学習に失敗しました")
                    else:
                        # セッションステートに保存
                        st.session_state['info_df'] = info_df
                        st.session_state['cols_to_keep'] = cols_to_keep
                        st.session_state['original_df'] = df.copy()
                        st.success("✨ 抽出完了！")

                # 抽出結果がある場合
                if 'info_df' in st.session_state and 'cols_to_keep' in st.session_state:
                    info_df = st.session_state['info_df']
                    cols_to_keep = st.session_state['cols_to_keep']
                    df = st.session_state['original_df'].copy()

                    # 統計情報表示
                    st.subheader("📈 抽出結果の統計")

                    cols = st.columns(len(cols_to_keep))
                    for i, col_name in enumerate(cols_to_keep):
                        with cols[i]:
                            count = info_df[col_name].str.len().gt(0).sum()
                            percentage = count / len(info_df) * 100

                            label = {
                                'gene_symbol': "Gene Symbol",
                                'refseq_id': "RefSeq ID",
                                'ensembl_id': "Ensembl ID",
                                'genbank_id': "GenBank ID",
                                'entrez_id': "Entrez ID",
                                'description': "Description",
                                'location': "Location"
                            }.get(col_name, col_name)

                            st.metric(
                                label,
                                f"{count:,}行",
                                f"{percentage:.1f}%"
                            )

                    # どの遺伝子名列を使用するか選択
                    st.subheader("🎯 1列目に配置する遺伝子名の選択")

                    # 利用可能な遺伝子名タイプを表示
                    available_types = []
                    for col_name in cols_to_keep:
                        if col_name in ['gene_symbol', 'refseq_id', 'ensembl_id', 'genbank_id', 'entrez_id']:
                            count = info_df[col_name].str.len().gt(0).sum()
                            if count > 0:
                                available_types.append(col_name)

                    if available_types:
                        selected_gene_col = st.selectbox(
                            "1列目に配置する遺伝子名タイプを選択してください",
                            options=available_types,
                            format_func=lambda x: {
                                'gene_symbol': "Gene Symbol",
                                'refseq_id': "RefSeq ID",
                                'ensembl_id': "Ensembl ID",
                                'genbank_id': "GenBank ID",
                                'entrez_id': "Entrez ID"
                            }.get(x, x),
                            help="選択した遺伝子名が行のインデックスになります"
                        )

                        # None行の削除オプション
                        remove_none = st.checkbox(
                            "遺伝子名が空の行を削除する",
                            value=True,
                            help="選択した遺伝子名タイプの値が空の行を削除します"
                        )

                        if st.button("✅ 遺伝子名を1列目に配置", type="primary"):
                            # 選択された遺伝子名列を取得
                            gene_names = info_df[selected_gene_col].copy()

                            # 他の抽出情報列を準備
                            other_cols = [col for col in cols_to_keep if col != selected_gene_col]
                            other_info = info_df[other_cols].copy()

                            # 元のデータと結合
                            result_df = pd.concat([other_info, df], axis=1)
                            result_df.index = gene_names
                            result_df.index.name = "Gene"

                            # None行の削除
                            if remove_none:
                                before_len = len(result_df)
                                result_df = result_df[result_df.index.str.len() > 0]
                                removed_count = before_len - len(result_df)
                                st.info(f"空の遺伝子名を持つ行を {removed_count} 行削除しました")

                            # セッションステートに保存
                            st.session_state['result_df'] = result_df
                            st.session_state['selected_gene_col'] = selected_gene_col

                            st.success("✨ 遺伝子名を1列目に配置しました！")

                        # 結果がある場合
                        if 'result_df' in st.session_state:
                            result_df = st.session_state['result_df']

                            st.subheader("📋 処理結果プレビュー")
                            st.write(f"データ長: {len(result_df):,}行")
                            st.dataframe(result_df.head(preview_rows))

                            # 重複遺伝子の検出
                            st.subheader("🔍 重複遺伝子の検出")

                            dup_d = result_df.loc[result_df.index.duplicated(keep=False), :].sort_index()
                            dup_count = len(set(dup_d.index))

                            if dup_count > 0:
                                st.warning(f"⚠️ 重複している遺伝子数: {dup_count}")
                                st.dataframe(dup_d)

                                # Aggregateオプション
                                st.subheader("📊 重複遺伝子の集約")

                                agg_method = st.radio(
                                    "集約方法を選択してください:",
                                    ('Mean', 'Max'),
                                    help="Mean: 平均値, Max: 最大値"
                                )

                                if st.button("🔄 重複遺伝子を集約", type="primary"):
                                    dup_gene = set(result_df.index[result_df.index.duplicated(keep=False)])
                                    df_nodup = result_df[~result_df.index.duplicated(keep='first')]
                                    grouping = result_df.groupby(level=0)

                                    if agg_method == "Mean":
                                        df_agg = grouping.mean(numeric_only=True)
                                    else:
                                        df_agg = grouping.max(numeric_only=True)

                                    for gene in dup_gene:
                                        # 数値列のみ更新
                                        df_nodup.loc[gene, df_agg.columns] = df_agg.loc[gene, df_agg.columns].values

                                    st.success(f"✨ {agg_method}方式で集約しました")
                                    st.write("集約後の重複遺伝子:")
                                    st.dataframe(df_nodup.loc[list(dup_gene), :].sort_index())

                                    # セッションステートを更新
                                    st.session_state['result_df'] = df_nodup
                                    st.session_state['aggregated'] = True

                                    st.rerun()
                            else:
                                st.success("✅ 重複遺伝子はありません")

                            # ダウンロード
                            st.markdown("---")
                            st.subheader("💾 ファイルダウンロード")

                            st.write(f"最終データ長: {len(result_df):,}行")

                            original_name = uploaded_file.name.rsplit('.', 1)[0]
                            output_filename = f"{original_name}_filtered.tsv"

                            csv_data = convert_df(result_df)
                            st.download_button(
                                label="📥 処理済みファイルをダウンロード (TSV)",
                                data=csv_data,
                                file_name=output_filename,
                                mime="text/tab-separated-values"
                            )
                    else:
                        st.warning("遺伝子名情報が抽出されませんでした")

        except Exception as e:
            st.error(f"❌ エラーが発生しました: {str(e)}")
            st.markdown("**デバッグ情報:**")
            st.code(str(e))

    # フッター
    st.markdown("---")
    st.markdown("""
    <div style='text-align: center; color: gray; font-size: 0.8em;'>
    Microarray Gene Name Filter v1.0 |
    マイクロアレイデータから遺伝子名を抽出・フィルタリング・集約
    </div>
    """, unsafe_allow_html=True)

if __name__ == "__main__":
    main()
