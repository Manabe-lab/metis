#!/usr/bin/env python3
"""
Microarray Gene Name Filter
maikuroareiDatafromGenenametheextractout,Filtering,gatheraboutdotsu-ru
"""

import streamlit as st
import pandas as pd
import numpy as np
import io
import re
from collections import defaultdict, Counter

class GeneInfoExtractor:
    """Geneinfoinfoofrankplacepata-nthelearnlearnshiteextractoutdokurasu"""

    def __init__(self):
        self.position_patterns = {}
        self.separator_pattern = None
        self.learned = False

    def detect_separator(self, sample_texts):
        """Sampletekisutofromuseusesareteexistareacutritextchartheselfmovecheckout"""
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

        # eachareacutritextcharofoutpresenttimenumthekaunto
        counts = {}
        for name, pattern in separators.items():
            total_count = 0
            for text in sample_texts:
                if pd.notna(text) and text != "---":
                    count = len(re.findall(pattern, str(text)))
                    total_count += count
            if total_count > 0:
                counts[name] = total_count

        # mostmofreqoutdoareacutritextcharthereturnsu
        if counts:
            most_common = max(counts, key=counts.get)
            return most_common, separators[most_common]
        else:
            return ' // ', r'\s*//\s*'  # defuoruto

    def classify_element(self, element):
        """needelemofkindtypethedivtype"""
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

        # Entrez Gene ID (numcharofmi)
        if element.isdigit():
            return "entrez"

        # 染colorbodyrankplaceinfoinfo
        if "|" in element or "cM" in element or element.startswith("chr"):
            return "location"

        # Descriptiontextofjudgeset
        words = element.split()
        if len(words) >= 3:
            # longiDescriptiontext
            if not re.match(r'^[A-Za-z0-9]+[-]?[A-Za-z0-9]+$', element):
                return "description"

        # specsetofki-wa-dotheincludemuDescriptiontext
        skip_keywords = ["uncharacterized", "predicted", "hypothetical", "family", "domain",
                        "containing", "protein", "member", "homolog", "like", "antigen", "receptor"]
        if any(keyword in element.lower() for keyword in skip_keywords):
            if len(element.split()) > 2 or len(element) > 30:
                return "description"

        # Geneshinboruofpossiblenatureishighi
        if re.match(r'^([A-Za-z][A-Za-z0-9_-]*\d*[A-Za-z0-9]*|[0-9]+[A-Za-z][A-Za-z0-9]*[A-Za-z]+)$', element):
            return "gene_symbol"

        # soofother（osorakuDescriptiontext）
        return "description"

    def learn_pattern(self, sample_data, separator=None):
        """SampleDatafromrankplacepata-nthelearnlearn"""
        # areacutritextcharthecheckout
        if separator is None:
            sample_texts = [str(text) for text in sample_data if pd.notna(text) and text != "---"]
            sep_name, self.separator_pattern = self.detect_separator(sample_texts[:20])
            st.info(f"🔍 checkoutsaretaareacutritextchar: '{sep_name}'")
        else:
            self.separator_pattern = separator

        # eachrankplacewithofneedelemtaipuofoutpresentfreqdegreethekaunto
        position_counts = defaultdict(Counter)
        max_positions = 0

        for text in sample_data:
            if pd.isna(text) or text == "---" or text == "":
                continue

            text = str(text)
            # multinumentori-（///areacutri）ofplacematchismostfirstofentori-ofmiuseuse
            if " /// " in text:
                text = text.split(" /// ")[0]

            # areacutritextcharwithdivratio
            parts = re.split(self.separator_pattern, text)
            max_positions = max(max_positions, len(parts))

            for pos, part in enumerate(parts):
                element_type = self.classify_element(part)
                if element_type != "empty":
                    position_counts[pos][element_type] += 1

        # eachrankplacewithmostmofreqoutdoneedelemtaiputhedetset
        self.position_patterns = {}

        st.write("📊 rankplacepata-ndivanalyzeResult:")
        pattern_summary = []

        for pos in range(max_positions):
            if pos in position_counts:
                # mostmofreqoutdotaiputhegetget
                most_common = position_counts[pos].most_common()
                if most_common:
                    top_type = most_common[0][0]
                    top_count = most_common[0][1]
                    total = sum(position_counts[pos].values())
                    confidence = top_count / total * 100

                    # trustrelydegreeis40%orupofplacematchofmipata-nandshite採use
                    if confidence >= 40:
                        self.position_patterns[pos] = top_type
                        pattern_summary.append({
                            "rankplace": pos + 1,
                            "Estimationtaipu": self.get_type_label(top_type),
                            "trustrelydegree": f"{confidence:.1f}%",
                            "Samplenum": total
                        })

        if pattern_summary:
            pattern_df = pd.DataFrame(pattern_summary)
            st.dataframe(pattern_df, use_container_width=True)

        self.learned = True
        return self.position_patterns

    def get_type_label(self, type_key):
        """taipuki-the日mainwordraberutochangechange"""
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
        """learnlearnshitapata-ntheuseuseshiteinfoinfotheextractout"""
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

        # multinumentori-（///areacutri）ofplacematchismostfirstofentori-ofmiuseuse
        if " /// " in text:
            text = text.split(" /// ")[0]

        # areacutritextcharwithdivratio
        parts = re.split(self.separator_pattern, text)

        # eachtaipugoandtovalthecollgather
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

            # learnlearnshitapata-nisexistplacematchispriorfirst
            if pos in self.position_patterns:
                expected_type = self.position_patterns[pos]
                actual_type = self.classify_element(part)

                # periodwaitsarerutaipuandrealoccoftaipuisonecausedoplacematch
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
                # taipuisonecauseshinotplacematchwithmo、realoccoftaiputobaseduiteSave
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
                # pata-nisnotrankplaceisrealoccoftaipuwithdivtype
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

        # Resultthemaandmeru
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
    """DataFrametheTSVshapeformattochangechange"""
    return df.to_csv(index=True, sep='\t').encode('utf-8')

@st.cache_data
def load_and_parse_file(file_content, file_name):
    """Filethereadmiintonwithpa-sudo（kiyashiyuattachki）"""
    delimiter = '\t' if '\t' in file_content.split('\n')[0] else ','
    df = pd.read_csv(
        io.StringIO(file_content),
        sep=delimiter,
        engine='python'
    )
    return df

@st.cache_data
def extract_gene_info(df_dict, gene_column, sample_size):
    """Geneinfoinfotheextractoutdo（kiyashiyuattachki）"""
    # dictfromDataFramethe復source
    df = pd.DataFrame(df_dict)

    extractor = GeneInfoExtractor()

    # SampleDatafromlearnlearn
    sample_size_actual = min(sample_size, len(df))
    sample_data = df[gene_column].head(sample_size_actual).dropna()
    extractor.learn_pattern(sample_data)

    if not extractor.learned:
        return None, None

    # bachiprocprocwitheffectrateize
    batch_size = 1000
    results = []

    for i in range(0, len(df), batch_size):
        batch = df[gene_column].iloc[i:i+batch_size]
        batch_results = batch.apply(extractor.extract_with_pattern)
        results.extend(batch_results.tolist())

    # extractoutshitainfoinfothepiecesepofcoltoexpandopen
    info_df = pd.DataFrame(results)

    # emptyofcolthedeleteremove
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
    ### Function
    - maikuroareiDatafromGeneinfoinfotheextractout
    - GenenamekaramutheSelectshite1colidxtodistplace
    - GenenameisNoneofrowthedeleteremovepossible
    - weightmultiGeneofcheckoutandgatherabout（Mean/Max）
    """)

    # FileUpload
    uploaded_file = st.file_uploader(
        "TSV/CSVFiletheSelectshitekudasai",
        type=['tsv', 'txt', 'csv'],
        help="maikuroareiDataFile"
    )

    if uploaded_file is not None:
        try:
            # FileLoading（kiyashiyuattachki）
            content = uploaded_file.getvalue().decode('utf-8')
            df = load_and_parse_file(content, uploaded_file.name)

            st.success(f"✅ FileLoadingSuccess: {len(df):,}row × {len(df.columns)}col")

            # Datapurebiyu-
            st.subheader("📊 Datapurebiyu-")

            preview_rows = st.slider(
                "purebiyu-rownum",
                min_value=5,
                max_value=min(100, len(df)),
                value=10
            )

            st.dataframe(df.head(preview_rows))

            # GenecolofSelect
            st.subheader("🔍 GeneinfoinfocolofSelect")

            # defuorutocolthesearchsu
            default_col = None
            for col in df.columns:
                if 'gene' in col.lower() or 'symbol' in col.lower() or 'description' in col.lower():
                    default_col = col
                    break

            gene_column = st.selectbox(
                "GeneinfoinfoisincludemarerucoltheSelectshitekudasai",
                options=df.columns.tolist(),
                index=df.columns.tolist().index(default_col) if default_col else 0,
                help="passnormalis'gene_assignment'coltoGeneinfoinfoisincludemareteimasu"
            )

            # SelectshitacolofSampleDisplay
            if gene_column:
                st.info("SelectshitacolofSample（mostfirstof5row）:")
                sample_data = df[gene_column].head(5).to_frame()
                st.dataframe(sample_data)

                # learnlearnSamplesaizuofSettings
                st.subheader("⚙️ learnlearnSettings")

                col1, col2 = st.columns(2)
                with col1:
                    sample_size = st.number_input(
                        "learnlearntouseusedoSamplenum",
                        min_value=20,
                        max_value=min(500, len(df)),
                        value=min(100, len(df)),
                        help="manyihodocorrectcertainwithsuisprocproctimebetweenisincreaseemasu"
                    )

                with col2:
                    st.info(f"""
                    📌 inferrecSettings:
                    - smallrulemodelData（<1000row）: 50-100
                    - midrulemodelData（1000-10000row）: 100-200
                    - bigrulemodelData（>10000row）: 200-500
                    """)

                # procprocRunbotan
                if st.button("🚀 Geneinfoinfotheextractout", type="primary"):
                    st.subheader("📈 rankplacepata-noflearnlearnandextractout")

                    with st.spinner("Geneinfoinfotheextractoutmid..."):
                        # DataFramethedicttochangechangeshitekiyashiyupossibletodo
                        df_dict = df.to_dict('list')
                        info_df, cols_to_keep = extract_gene_info(df_dict, gene_column, sample_size)

                    if info_df is None:
                        st.error("pata-noflearnlearntofailfailshimashita")
                    else:
                        # seshiyonsute-totoSave
                        st.session_state['info_df'] = info_df
                        st.session_state['cols_to_keep'] = cols_to_keep
                        st.session_state['original_df'] = df.copy()
                        st.success("✨ extractoutComplete！")

                # extractoutResultisexistplacematch
                if 'info_df' in st.session_state and 'cols_to_keep' in st.session_state:
                    info_df = st.session_state['info_df']
                    cols_to_keep = st.session_state['cols_to_keep']
                    df = st.session_state['original_df'].copy()

                    # StatisticalinfoinfoDisplay
                    st.subheader("📈 extractoutResultofStatistical")

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
                                f"{count:,}row",
                                f"{percentage:.1f}%"
                            )

                    # doofGenenamecoltheuseusedokaSelect
                    st.subheader("🎯 1colidxtodistplacedoGenenameofSelect")

                    # useusepossiblenaGenenametaiputheDisplay
                    available_types = []
                    for col_name in cols_to_keep:
                        if col_name in ['gene_symbol', 'refseq_id', 'ensembl_id', 'genbank_id', 'entrez_id']:
                            count = info_df[col_name].str.len().gt(0).sum()
                            if count > 0:
                                available_types.append(col_name)

                    if available_types:
                        selected_gene_col = st.selectbox(
                            "1colidxtodistplacedoGenenametaiputheSelectshitekudasai",
                            options=available_types,
                            format_func=lambda x: {
                                'gene_symbol': "Gene Symbol",
                                'refseq_id': "RefSeq ID",
                                'ensembl_id': "Ensembl ID",
                                'genbank_id': "GenBank ID",
                                'entrez_id': "Entrez ID"
                            }.get(x, x),
                            help="SelectshitaGenenameisrowofindekusutonarimasu"
                        )

                        # NonerowofdeleteremoveOption
                        remove_none = st.checkbox(
                            "Genenameisemptyofrowthedeleteremovedo",
                            value=True,
                            help="SelectshitaGenenametaipuofvalisemptyofrowthedeleteremoveshimasu"
                        )

                        if st.button("✅ Genenamethe1colidxtodistplace", type="primary"):
                            # SelectsaretaGenenamecolthegetget
                            gene_names = info_df[selected_gene_col].copy()

                            # otherofextractoutinfoinfocolthelevelprep
                            other_cols = [col for col in cols_to_keep if col != selected_gene_col]
                            other_info = info_df[other_cols].copy()

                            # sourceofDataandresultmatch
                            result_df = pd.concat([other_info, df], axis=1)
                            result_df.index = gene_names
                            result_df.index.name = "Gene"

                            # Nonerowofdeleteremove
                            if remove_none:
                                before_len = len(result_df)
                                result_df = result_df[result_df.index.str.len() > 0]
                                removed_count = before_len - len(result_df)
                                st.info(f"emptyofGenenametheholdtsurowthe {removed_count} rowdeleteremoveshimashita")

                            # seshiyonsute-totoSave
                            st.session_state['result_df'] = result_df
                            st.session_state['selected_gene_col'] = selected_gene_col

                            st.success("✨ Genenamethe1colidxtodistplaceshimashita！")

                        # Resultisexistplacematch
                        if 'result_df' in st.session_state:
                            result_df = st.session_state['result_df']

                            st.subheader("📋 procprocResultpurebiyu-")
                            st.write(f"Datalong: {len(result_df):,}row")
                            st.dataframe(result_df.head(preview_rows))

                            # weightmultiGeneofcheckout
                            st.subheader("🔍 weightmultiGeneofcheckout")

                            dup_d = result_df.loc[result_df.index.duplicated(keep=False), :].sort_index()
                            dup_count = len(set(dup_d.index))

                            if dup_count > 0:
                                st.warning(f"⚠️ weightmultishiteexistGenenum: {dup_count}")
                                st.dataframe(dup_d)

                                # AggregateOption
                                st.subheader("📊 weightmultiGeneofgatherabout")

                                agg_method = st.radio(
                                    "gatheraboutwaymethodtheSelectshitekudasai:",
                                    ('Mean', 'Max'),
                                    help="Mean: flatavgval, Max: maximumval"
                                )

                                if st.button("🔄 weightmultiGenethegatherabout", type="primary"):
                                    dup_gene = set(result_df.index[result_df.index.duplicated(keep=False)])
                                    df_nodup = result_df[~result_df.index.duplicated(keep='first')]
                                    grouping = result_df.groupby(level=0)

                                    if agg_method == "Mean":
                                        df_agg = grouping.mean(numeric_only=True)
                                    else:
                                        df_agg = grouping.max(numeric_only=True)

                                    for gene in dup_gene:
                                        # numvalcolofmifurthernew
                                        df_nodup.loc[gene, df_agg.columns] = df_agg.loc[gene, df_agg.columns].values

                                    st.success(f"✨ {agg_method}wayformatwithgatheraboutshimashita")
                                    st.write("gatheraboutafterofweightmultiGene:")
                                    st.dataframe(df_nodup.loc[list(dup_gene), :].sort_index())

                                    # seshiyonsute-tothefurthernew
                                    st.session_state['result_df'] = df_nodup
                                    st.session_state['aggregated'] = True

                                    st.rerun()
                            else:
                                st.success("✅ weightmultiGeneisarimasen")

                            # Download
                            st.markdown("---")
                            st.subheader("💾 FileDownload")

                            st.write(f"finalDatalong: {len(result_df):,}row")

                            original_name = uploaded_file.name.rsplit('.', 1)[0]
                            output_filename = f"{original_name}_filtered.tsv"

                            csv_data = convert_df(result_df)
                            st.download_button(
                                label="📥 procprocdonemiFiletheDownload (TSV)",
                                data=csv_data,
                                file_name=output_filename,
                                mime="text/tab-separated-values"
                            )
                    else:
                        st.warning("Genenameinfoinfoisextractoutsaremasenwithshita")

        except Exception as e:
            st.error(f"❌ Errorisoccurgenshimashita: {str(e)}")
            st.markdown("**debaguinfoinfo:**")
            st.code(str(e))

    # futa-
    st.markdown("---")
    st.markdown("""
    <div style='text-align: center; color: gray; font-size: 0.8em;'>
    Microarray Gene Name Filter v1.0 |
    maikuroareiDatafromGenenametheextractout,Filtering,gatherabout
    </div>
    """, unsafe_allow_html=True)

if __name__ == "__main__":
    main()
