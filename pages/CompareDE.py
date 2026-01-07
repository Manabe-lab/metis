import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from upsetplot import plot
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import io
from difflib import SequenceMatcher
import itertools

def process_filename(filenames):
    """
    Filenamethe2tsuzutsuComparisonshi、copassdopartdivtheremovermdo
    """
    # expandextendchildtheremoverm
    names = [name.rsplit('.', 1)[0] for name in filenames]
    
    def compare_two_names(name1, name2):
        """2tsuofnamebeforetheComparisonshi、sorezoreoffixhaveofpartdivtheviewtsukeru"""
        # namebeforethedivratioshi、emptyofneedelemtheremoverm
        parts1 = [p for p in name1.replace('-', '_').split('_') if p]
        parts2 = [p for p in name2.replace('-', '_').split('_') if p]
        
        # copasspartdivtheviewtsukeru
        common = set(parts1) & set(parts2)
        
        # eachnamebeforeoffixhavepartdivtheextractout
        unique1 = [p for p in parts1 if p not in common]
        unique2 = [p for p in parts2 if p not in common]
        
        return unique1, unique2

    # eachFileoffixhavepartdivthekeephold
    processed = {f: set() for f in filenames}
    
    # 2tsuzutsuComparison
    for i in range(len(filenames)):
        for j in range(i + 1, len(filenames)):
            unique1, unique2 = compare_two_names(names[i], names[j])
            # fixhavepartdivthe蓄accumulate
            processed[filenames[i]].update(unique1)
            processed[filenames[j]].update(unique2)
    
    # Resultthekuri-napu
    result = {}
    for f in filenames:
        unique_parts = list(processed[f])
        if unique_parts:
            # emptyofneedelemtheremovermshiteresultmatch
            result[f] = '_'.join(filter(None, unique_parts)).strip('_')
        else:
            # sourceofnamebeforefromemptyofneedelemtheremoverm
            parts = [p for p in names[filenames.index(f)].replace('-', '_').split('_') if p]
            result[f] = '_'.join(parts)
    
    return result
    

@st.cache_data
def read_excel(file, index_col=None, header=0):
    return pd.read_excel(file, index_col=index_col, header=header)

@st.cache_data
def read_csv(file, index_col=None, sep=',', header=0):
    return pd.read_csv(file, index_col=index_col, header=header, sep=sep, engine='python')

def get_regulated_genes(df, lfc_col, pval_col, direction='up', 
                       lfc_threshold=1, pval_threshold=0.05):
    """pointsetshitawaydirtoExpressionchangemovedoGenetheextractoutdorelnum"""
    if direction == 'up':
        return set(df.index[
            (df[lfc_col] > lfc_threshold) & 
            (df[pval_col] < pval_threshold)
        ])
    else:  # down
        return set(df.index[
            (df[lfc_col] < -lfc_threshold) & 
            (df[pval_col] < pval_threshold)
        ])

def create_venn_diagram(datasets, names, title):
    """benfigthemakebecomedorelnum"""
    fig, ax = plt.subplots(figsize=(10, 10))
    
    if len(datasets) == 2:
        venn2(
            subsets=[datasets[0], datasets[1]],
            set_labels=names,
            ax=ax
        )
    elif len(datasets) == 3:
        venn3(
            subsets=[datasets[0], datasets[1], datasets[2]],
            set_labels=names,
            ax=ax
        )
    
    plt.title(title)
    return fig

def create_upset_plot(datasets, names, title):
    """
    UpSetpurotowindofVisualizationthemakebecomedorelnum
    """
    # copasspartdivofsetmimatchwasetheviewtsukeru
    def get_combinations():
        n = len(datasets)
        combos = []
        sizes = []
        labels = []
        
        # eachsetmimatchwasetotsuite
        for i in range(1, n + 1):
            for combo in itertools.combinations(range(n), i):
                # soofsetmimatchwasetoincludemareruDatasetoofcopasspartdiv
                common = set.intersection(*[datasets[j] for j in combo])
                if len(common) > 0:  # copasspartdivisexistatdoplacematchofmiaddadd
                    combos.append(combo)
                    sizes.append(len(common))
                    labels.append(' & '.join([names[j] for j in combo]))
        
        return combos, sizes, labels

    combos, sizes, labels = get_combinations()
    
    # saizuwithso-to
    sorted_indices = np.argsort(sizes)[::-1]
    sizes = [sizes[i] for i in sorted_indices]
    labels = [labels[i] for i in sorted_indices]
    combos = [combos[i] for i in sorted_indices]

    # purotoofmakebecome
    fig, (ax_sets, ax_matrix) = plt.subplots(2, 1, 
                                            figsize=(12, 8),
                                            gridspec_kw={'height_ratios': [3, 1]})
    
    # ba-chiya-toofplotdraw
    bars = ax_sets.bar(range(len(sizes)), sizes)
    ax_sets.set_xticks([])
    ax_sets.set_ylabel('Intersection Size')
    
    # rowcolofplotdraw
    matrix = np.zeros((len(names), len(combos)))
    for i, combo in enumerate(combos):
        for j in combo:
            matrix[j, i] = 1
    
    ax_matrix.imshow(matrix, cmap='binary', aspect='auto')
    ax_matrix.set_yticks(range(len(names)))
    ax_matrix.set_yticklabels(names)
    ax_matrix.set_xticks(range(len(labels)))
    ax_matrix.set_xticklabels(labels, rotation=45, ha='right')
    
    # reiautoofadjustarrange
    plt.tight_layout()
    ax_sets.set_title(title)
    
    return fig


def create_2d_scatter(df1, df2, x_col, y_col, name1, name2, pval_col1, pval_col2,
                     highlight_significant=True, lfc_threshold=1, pval_threshold=0.05,
                     use_union=False):
    """2nextsourcescatterdistfigthemakebecome"""
    # indekusuofprocproc（copass or sumgathermatch）
    if use_union:
        # sumgathermatchtheuseuseshi、misslossvalthe0withfillmeru
        all_index = df1.index.union(df2.index)
        df1 = df1.reindex(all_index).fillna(0)
        df2 = df2.reindex(all_index).fillna(0)
    else:
        # copasspartdivofmiuseuse
        common_index = df1.index.intersection(df2.index)
        df1 = df1.loc[common_index]
        df2 = df2.loc[common_index]

    fig = go.Figure()
    
    # hoba-tekisutoofmakebecome
    hovertemplate = (
        "Gene: %{customdata}<br>" +
        f"{name1} {x_col}: %{{x:.3f}}<br>" +
        f"{name2} {y_col}: %{{y:.3f}}<br>" +
        "<extra></extra>"
    )

    if 'neg_log_padj' not in x_col:  # Log2FC ofscatterdistfigofplacematch
        if highlight_significant:
            # Dataseto1withSignificant
            sig1 = (abs(df1[x_col]) > lfc_threshold) & (df1[pval_col1] < pval_threshold)
            # Dataseto2withSignificant
            sig2 = (abs(df2[y_col]) > lfc_threshold) & (df2[pval_col2] < pval_threshold)
            
            # bothwaywithSignificant
            both_sig = sig1 & sig2
            # dochirakaofmiwithSignificant
            only_sig1 = sig1 & ~sig2
            only_sig2 = ~sig1 & sig2
            # nonSignificant
            not_sig = ~sig1 & ~sig2

            # eachkategori-ofpuroto
            fig.add_trace(go.Scatter(
                x=df1[x_col][both_sig],
                y=df2[y_col][both_sig],
                mode='markers',
                name='Both significant',
                customdata=df1.index[both_sig],
                hovertemplate=hovertemplate,
                marker=dict(color='red', size=4)
            ))
            fig.add_trace(go.Scatter(
                x=df1[x_col][only_sig1],
                y=df2[y_col][only_sig1],
                mode='markers',
                name=f'Only in {name1}',
                customdata=df1.index[only_sig1],
                hovertemplate=hovertemplate,
                marker=dict(color='blue', size=4)
            ))
            fig.add_trace(go.Scatter(
                x=df1[x_col][only_sig2],
                y=df2[y_col][only_sig2],
                mode='markers',
                name=f'Only in {name2}',
                customdata=df1.index[only_sig2],
                hovertemplate=hovertemplate,
                marker=dict(color='green', size=4)
            ))
            fig.add_trace(go.Scatter(
                x=df1[x_col][not_sig],
                y=df2[y_col][not_sig],
                mode='markers',
                name='Not significant',
                customdata=df1.index[not_sig],
                hovertemplate=hovertemplate,
                marker=dict(color='grey', size=2, opacity=0.5)
            ))
        else:
            fig.add_trace(go.Scatter(
                x=df1[x_col],
                y=df2[y_col],
                mode='markers',
                name='Genes',
                customdata=df1.index,
                hovertemplate=hovertemplate,
                marker=dict(
                    color='blue',
                    size=3,
                    opacity=0.6
                )
            ))
    else:  # Pvalofscatterdistfigofplacematch
        # Pvaltoyorucolordivke
        low_p1 = df1[pval_col1] < pval_threshold
        low_p2 = df2[pval_col2] < pval_threshold
        
        # bothwaywithSignificant
        both_sig = low_p1 & low_p2
        # dochirakaofmiwithSignificant
        only_sig1 = low_p1 & ~low_p2
        only_sig2 = ~low_p1 & low_p2
        # nonSignificant
        not_sig = ~low_p1 & ~low_p2
        
        # eachkategori-ofpuroto
        fig.add_trace(go.Scatter(
            x=df1[x_col][both_sig],
            y=df2[y_col][both_sig],
            mode='markers',
            name='Both significant',
            customdata=df1.index[both_sig],
            hovertemplate=hovertemplate,
            marker=dict(color='red', size=4)
        ))
        fig.add_trace(go.Scatter(
            x=df1[x_col][only_sig1],
            y=df2[y_col][only_sig1],
            mode='markers',
            name=f'Only in {name1}',
            customdata=df1.index[only_sig1],
            hovertemplate=hovertemplate,
            marker=dict(color='blue', size=4)
        ))
        fig.add_trace(go.Scatter(
            x=df1[x_col][only_sig2],
            y=df2[y_col][only_sig2],
            mode='markers',
            name=f'Only in {name2}',
            customdata=df1.index[only_sig2],
            hovertemplate=hovertemplate,
            marker=dict(color='green', size=4)
        ))
        fig.add_trace(go.Scatter(
            x=df1[x_col][not_sig],
            y=df2[y_col][not_sig],
            mode='markers',
            name='Not significant',
            customdata=df1.index[not_sig],
            hovertemplate=hovertemplate,
            marker=dict(color='grey', size=2, opacity=0.5)
        ))

    # supp助lineofaddadd
    fig.add_hline(y=0, line_dash="dash", line_color="gray")
    fig.add_vline(x=0, line_dash="dash", line_color="gray")

    fig.update_layout(
        title=f"Comparison of {name1} vs {name2}",
        xaxis_title=f"{name1} - {x_col}",
        yaxis_title=f"{name2} - {y_col}",
        hovermode='closest',
        showlegend=True
    )

    return fig


def create_3d_scatter(df1, df2, df3, x_col, y_col, z_col, name1, name2, name3,
                     pval_col1, pval_col2, pval_col3,
                     highlight_significant=True, lfc_threshold=1, pval_threshold=0.05,
                     use_union=False):
    """3nextsourcescatterdistfigthemakebecome"""
    # indekusuofprocproc（copass or sumgathermatch）
    if use_union:
        # sumgathermatchtheuseuseshi、misslossvalthe0withfillmeru
        all_index = df1.index.union(df2.index).union(df3.index)
        df1 = df1.reindex(all_index).fillna(0)
        df2 = df2.reindex(all_index).fillna(0)
        df3 = df3.reindex(all_index).fillna(0)
    else:
        # copasspartdivofmiuseuse
        common_index = df1.index.intersection(df2.index).intersection(df3.index)
        df1 = df1.loc[common_index]
        df2 = df2.loc[common_index]
        df3 = df3.loc[common_index]

    fig = go.Figure()
    
    # hoba-tekisutoofmakebecome
    hovertemplate = (
        "Gene: %{customdata}<br>" +
        f"{name1} {x_col}: %{{x:.3f}}<br>" +
        f"{name2} {y_col}: %{{y:.3f}}<br>" +
        f"{name3} {z_col}: %{{z:.3f}}<br>" +
        "<extra></extra>"
    )

    if highlight_significant:
        # eachDatasetowithofSignificantnature
        sig1 = (abs(df1[x_col]) > lfc_threshold) & (df1[pval_col1] < pval_threshold)
        sig2 = (abs(df2[y_col]) > lfc_threshold) & (df2[pval_col2] < pval_threshold)
        sig3 = (abs(df3[z_col]) > lfc_threshold) & (df3[pval_col3] < pval_threshold)

        # alltewithSignificant
        all_sig = sig1 & sig2 & sig3
        # soreorout
        other = ~all_sig

        # puroto
        fig.add_trace(go.Scatter3d(
            x=df1[x_col][all_sig],
            y=df2[y_col][all_sig],
            z=df3[z_col][all_sig],
            mode='markers',
            name='Significant in all',
            customdata=df1.index[all_sig],
            hovertemplate=hovertemplate,
            marker=dict(size=4, color='red')
        ))
        fig.add_trace(go.Scatter3d(
            x=df1[x_col][other],
            y=df2[y_col][other],
            z=df3[z_col][other],
            mode='markers',
            name='Other',
            customdata=df1.index[other],
            hovertemplate=hovertemplate,
            marker=dict(size=2, color='grey', opacity=0.5)
        ))
    else:
        fig.add_trace(go.Scatter3d(
            x=df1[x_col],
            y=df2[y_col],
            z=df3[z_col],
            mode='markers',
            name='Genes',
            customdata=df1.index,
            hovertemplate=hovertemplate,
            marker=dict(size=3, color='blue', opacity=0.6)
        ))

    fig.update_layout(
        title="3D Comparison",
        scene=dict(
            xaxis_title=f"{name1} - {x_col}",
            yaxis_title=f"{name2} - {y_col}",
            zaxis_title=f"{name3} - {z_col}"
        ),
        margin=dict(l=0, r=0, b=0, t=30)
    )

    return fig


def main():
    st.title("DE results comparison")
    
    # FileUpload
    uploaded_files = st.file_uploader(
        "DEAnalysisResultFiletheSelectshitekudasai（multinumcan）", 
        accept_multiple_files=True,
        type=['csv', 'tsv', 'txt', 'xlsx', 'xls']
    )
    st.write("P-value, FCtheholdtsuFilethe２tsu")
    
    if not uploaded_files:
        st.info("AnalysisResultFiletheUploadshitekudasai。")
        return
    
    # Filenameofprocproc
    filename_mapping = process_filename([f.name for f in uploaded_files])
    
    # ThresholdofSettings
    col1, col2 = st.columns(2)
    with col1:
        lfc_threshold = st.number_input(
            "Log2 Fold Change threshold",
            value=1.0,
            step=0.1,
            min_value = 0.0
        )
    with col2:
        pval_threshold = st.number_input(
            "FDR threshold",
            value=0.05,
            step=0.01,
            min_value = 0.00,
            max_value = 1.00,
            format="%.3f"
        )
    
    # DatasetoofLoadingandprocproc
    datasets = {}
    for uploaded_file in uploaded_files:
        try:
            try:
                df = read_csv(uploaded_file, index_col=0, sep=None)
            except:
                df = read_excel(uploaded_file, index_col=0)
            numeric_cols = df.select_dtypes(include=[np.number]).columns
            datasets[filename_mapping[uploaded_file.name]] = {
                'data': df,
                'columns': numeric_cols
            }
        except Exception as e:
            st.error(f"Error ({uploaded_file.name}): {e}")
    
    if len(datasets) < 2:
        st.warning("Comparisontoisfewnakuandmo2tsuofDatasetoisrequiredwithsu。")
        return
    
    # DatasetoSettings
    st.subheader("DatasetoandkaramuofSettings")
    
    dataset_settings = []
    for name, dataset in datasets.items():
        st.markdown(f"#### {name}")
        # karamunameofpata-nmachingu
        pvalue = [i for i in dataset['columns'] if ('adj' in i.lower()) or  ('fdr' in i.lower()) or
         ('pvalue' in i.lower()) or ('p-val' in i.lower()) or  ('p val' in i.lower()) 
                 ]




        fc = [i for i in dataset['columns'] if ('log2fc' in i.lower()) or 
              ('fold change' in i.lower()) or ('log2foldchange' in i.lower()) or 
              ('logfc' in i.lower()) or ('coef' in i.lower())]
        
        lfc_col = st.selectbox(
            "Log2 Fold Change",
            fc,
            key=f"lfc_{name}"
        )

        pval_col = st.selectbox(
            "adjusted P",
            pvalue,
            key=f"pval_{name}"
        )
        use_in_analysis = True
        invert_sign = st.checkbox(
            "Log2 Fold Changeof符numtheantiturn",
            key=f"invert_{name}",
            help="kontoro-ruofsetdefisinvofplacematchtochiekushitekudasai"
        )
        
        # DatasetoofSettingstheSave
        dataset_settings.append({
            'name': name,
            'data': dataset['data'].copy(),  # Create a copy to avoid modifying original
            'lfc_col': lfc_col,
            'pval_col': pval_col,
            'invert_sign': invert_sign,
            'use_in_analysis': use_in_analysis
        })
        
        # Invert sign if requested
        if invert_sign:
            dataset_settings[-1]['data'][lfc_col] = -dataset_settings[-1]['data'][lfc_col]

        st.markdown("---")
    
    # AnalysistaipuofSelect
    analysis_type = st.radio(
        "AnalysistaiputheSelect：",
        ["ExpressionchangemoveGeneofComparison", "scatterdistfigtoyoruComparison"],
        horizontal=True
    )
    
    active_datasets = [ds for ds in dataset_settings if ds['use_in_analysis']]
    
    if analysis_type == "ExpressionchangemoveGeneofComparison":
        # ExpressionwaydirnatureofSelect
        direction = st.radio(
            "ExpressionchangemoveofwaydirnaturetheSelect：",
            ["Up-regulated", "Down-regulated", "Both"],
            horizontal=True
        )
        
        # purototaipuofSelect
        plot_type = st.radio(
            "purototaiputheSelect：",
            ["benfig", "UpSetpuroto"],
            horizontal=True
        )
        
        if direction in ["Up-regulated", "Down-regulated"]:
            dir_type = 'up' if direction == "Up-regulated" else 'down'
            significant_sets = []
            set_names = []
            
            for ds in active_datasets:
                sig_genes = get_regulated_genes(
                    ds['data'], 
                    ds['lfc_col'], 
                    ds['pval_col'],
                    dir_type,
                    lfc_threshold,
                    pval_threshold
                )
                significant_sets.append(sig_genes)
                set_names.append(ds['name'])
            
            title = f"{direction}"
            
            if plot_type == "benfig" and len(significant_sets) <= 3:
                fig = create_venn_diagram(significant_sets, set_names, title)
                st.pyplot(fig)
            else:
                fig = create_upset_plot(significant_sets, set_names, title)
                st.pyplot(fig)
            
            # copassGeneofDisplay
            common_genes = set.intersection(*significant_sets)
            st.write(f"copassshite{direction}doGenenum: {len(common_genes)}")
            
            if st.checkbox(f"copass{direction}GenerisutotheDisplay"):
                st.write(", ".join(sorted(list(common_genes))))
            
            # copassGeneofekusupo-to
            if len(common_genes) > 0:
                csv = pd.DataFrame(sorted(list(common_genes)), columns=['Gene_ID'])
                csv_data = csv.to_csv(index=False)
                st.download_button(
                    label=f"copass{direction}GenerisutotheDownload",
                    data=csv_data,
                    file_name=f"common_{dir_type}_regulated_genes.csv",
                    mime="text/csv"
                )
        
        else:  # bothwaysep々toDisplay
            for dir_type, dir_name in [('up', 'Up-regulated'), ('down', 'Down-regulated')]:
                st.subheader(f"{dir_name}GeneofComparison")
                significant_sets = []
                set_names = []
                
                for ds in active_datasets:
                    sig_genes = get_regulated_genes(
                        ds['data'], 
                        ds['lfc_col'], 
                        ds['pval_col'],
                        dir_type,
                        lfc_threshold,
                        pval_threshold
                    )
                    significant_sets.append(sig_genes)
                    set_names.append(ds['name'])
                
                title = f"{dir_name}GeneofComparison"
                
                if plot_type == "benfig" and len(significant_sets) <= 3:
                    fig = create_venn_diagram(significant_sets, set_names, title)
                    st.pyplot(fig)
                else:
                    fig = create_upset_plot(significant_sets, set_names, title)
                    st.pyplot(fig)
                
                # copassGeneofDisplay
                common_genes = set.intersection(*significant_sets)
                st.write(f"copassshite{dir_name}doGenenum: {len(common_genes)}")
                
                if st.checkbox(f"copass{dir_name}GenerisutotheDisplay"):
                    st.write(", ".join(sorted(list(common_genes))))
                
                # copassGeneofekusupo-to
                if len(common_genes) > 0:
                    csv = pd.DataFrame(sorted(list(common_genes)), columns=['Gene_ID'])
                    csv_data = csv.to_csv(index=False)
                    st.download_button(
                        label=f"copass{dir_name}GenerisutotheDownload",
                        data=csv_data,
                        file_name=f"common_{dir_type}_regulated_genes.csv",
                        mime="text/csv"
                    )
    
    else:  # scatterdistfigtoyoruComparison
        st.subheader("scatterdistfigtoyoruComparison")
        
        plot_dimension = st.radio(
            "scatterdistfigofnextsourcetheSelect：",
            ["2Dscatterdistfig", "3Dscatterdistfig"],
            horizontal=True
        )
        
        highlight = st.checkbox("SignificantnaGenethehairaitoDisplay", value=True)
        
        if plot_dimension == "2Dscatterdistfig" and len(active_datasets) >= 2:
            col1, col2 = st.columns(2)
            with col1:
                dataset1 = st.selectbox("1tsuidxofDataseto", [ds['name'] for ds in active_datasets], key='scatter1', index=0)
            with col2:
                dataset2 = st.selectbox("2tsuidxofDataseto", [ds['name'] for ds in active_datasets], key='scatter2', index=1)
            
            if dataset1 != dataset2:
                ds1 = next(ds for ds in active_datasets if ds['name'] == dataset1)
                ds2 = next(ds for ds in active_datasets if ds['name'] == dataset2)
                
                # GenesetoofSelect
                use_union = st.checkbox("allteofGenetheDisplay（existatshinotGeneis0andshitehandleu）", 
                                      value=False,
                                      help="chiekutheoutsuandcopassofGeneofmitheDisplayshimasu")
                
                # GenenumofDisplay
                common_genes = len(ds1['data'].index.intersection(ds2['data'].index))
                total_genes = len(ds1['data'].index.union(ds2['data'].index))
                st.write(f"copassofGenenum: {common_genes}")
                st.write(f"totalGenenum: {total_genes}")
                
                # valofSelect
                value_type = st.radio(
                    "ComparisondovaltheSelect：",
                    ["Log2 Fold Change", "adjusted P-value (-log10)"],
                    horizontal=True
                )
                
                if value_type == "Log2 Fold Change":
                    fig = create_2d_scatter(
                        ds1['data'], ds2['data'],
                        ds1['lfc_col'], ds2['lfc_col'],
                        dataset1, dataset2,
                        ds1['pval_col'], ds2['pval_col'],
                        highlight,
                        lfc_threshold,
                        pval_threshold,
                        use_union
                    )
                else:
                    # adjusted P-valueofplacematchis-log10tochangechange
                    ds1['data']['neg_log_padj'] = -np.log10(ds1['data'][ds1['pval_col']])
                    ds2['data']['neg_log_padj'] = -np.log10(ds2['data'][ds2['pval_col']])
                    fig = create_2d_scatter(
                        ds1['data'], ds2['data'],
                        'neg_log_padj', 'neg_log_padj',
                        dataset1, dataset2,
                        ds1['pval_col'], ds2['pval_col'],
                        False,
                        lfc_threshold,
                        pval_threshold,
                        use_union
                    )
                
                st.plotly_chart(fig, use_container_width=True)

                # mutualrelrelatenumofCalculationandDisplay
                if value_type == "Log2 Fold Change":
                    if use_union:
                        all_index = ds1['data'].index.union(ds2['data'].index)
                        df1 = ds1['data'].reindex(all_index).fillna(0)
                        df2 = ds2['data'].reindex(all_index).fillna(0)
                    else:
                        common_index = ds1['data'].index.intersection(ds2['data'].index)
                        df1 = ds1['data'].loc[common_index]
                        df2 = ds2['data'].loc[common_index]
                    
                    # NaNtheremoveoutshitemutualrelrelatenumtheCalculation
                    mask = ~np.isnan(df1[ds1['lfc_col']]) & ~np.isnan(df2[ds2['lfc_col']])
                    correlation = np.corrcoef(
                        df1[ds1['lfc_col']][mask].astype(float),
                        df2[ds2['lfc_col']][mask].astype(float)
                    )[0,1]
                    
                    # mutualrelrelatenumisCalculationwithkitaplacematchofmiDisplay
                    if not np.isnan(correlation):
                        st.write(f"Log2 Fold Change ofmutualrelrelatenum: {correlation:.3f}")
                    else:
                        st.warning("Log2 Fold Change ofmutualrelrelatenumtheCalculationwithkimasenwithshita")
                else:
                    if use_union:
                        all_index = ds1['data'].index.union(ds2['data'].index)
                        df1 = ds1['data'].reindex(all_index).fillna(1)  # Pvalofmisslossis1andshitehandleu
                        df2 = ds2['data'].reindex(all_index).fillna(1)
                    else:
                        common_index = ds1['data'].index.intersection(ds2['data'].index)
                        df1 = ds1['data'].loc[common_index]
                        df2 = ds2['data'].loc[common_index]
                    
                    df1['neg_log_padj'] = -np.log10(df1[ds1['pval_col']])
                    df2['neg_log_padj'] = -np.log10(df2[ds2['pval_col']])
                    
                    correlation = np.corrcoef(
                        df1['neg_log_padj'].values,
                        df2['neg_log_padj'].values
                    )[0,1]
                    st.write(f"-log10(adjusted P-value) ofmutualrelrelatenum: {correlation:.3f}")
                
        elif plot_dimension == "3Dscatterdistfig" and len(active_datasets) >= 3:
            col1, col2, col3 = st.columns(3)
            with col1:
                dataset1 = st.selectbox("XaxisofDataseto", [ds['name'] for ds in active_datasets], key='scatter3d1')
            with col2:
                dataset2 = st.selectbox("YaxisofDataseto", [ds['name'] for ds in active_datasets], key='scatter3d2')
            with col3:
                dataset3 = st.selectbox("ZaxisofDataseto", [ds['name'] for ds in active_datasets], key='scatter3d3')
            
            if len({dataset1, dataset2, dataset3}) == 3:
                ds1 = next(ds for ds in active_datasets if ds['name'] == dataset1)
                ds2 = next(ds for ds in active_datasets if ds['name'] == dataset2)
                ds3 = next(ds for ds in active_datasets if ds['name'] == dataset3)
                
                # GenesetoofSelect
                use_union = st.checkbox("allteofGenetheDisplay（existatshinotGeneis0andshitehandleu）", 
                                      value=False,
                                      help="chiekutheoutsuandcopassofGeneofmitheDisplayshimasu")
                
                # GenenumofDisplay
                common_genes = len(set.intersection(
                    set(ds1['data'].index),
                    set(ds2['data'].index),
                    set(ds3['data'].index)
                ))
                total_genes = len(set.union(
                    set(ds1['data'].index),
                    set(ds2['data'].index),
                    set(ds3['data'].index)
                ))
                st.write(f"copassofGenenum: {common_genes}")
                st.write(f"totalGenenum: {total_genes}")
                

                value_type = st.radio(
                    "ComparisondovaltheSelect：",
                    ["Log2 Fold Change", "adjusted P-value (-log10)"],
                    horizontal=True
                )
                
                if value_type == "Log2 Fold Change":
                    fig = create_3d_scatter(
                        ds1['data'], ds2['data'], ds3['data'],
                        ds1['lfc_col'], ds2['lfc_col'], ds3['lfc_col'],
                        dataset1, dataset2, dataset3,
                        ds1['pval_col'], ds2['pval_col'], ds3['pval_col'],
                        highlight,
                        lfc_threshold,
                        pval_threshold,
                        use_union
                    )
                else:
                    # adjusted P-valueofplacematchis-log10tochangechange
                    ds1['data']['neg_log_padj'] = -np.log10(ds1['data'][ds1['pval_col']])
                    ds2['data']['neg_log_padj'] = -np.log10(ds2['data'][ds2['pval_col']])
                    ds3['data']['neg_log_padj'] = -np.log10(ds3['data'][ds3['pval_col']])
                    fig = create_3d_scatter(
                        ds1['data'], ds2['data'], ds3['data'],
                        'neg_log_padj', 'neg_log_padj', 'neg_log_padj',
                        dataset1, dataset2, dataset3,
                        ds1['pval_col'], ds2['pval_col'], ds3['pval_col'],
                        False,
                        lfc_threshold,
                        pval_threshold,
                        use_union
                    )
                
                st.plotly_chart(fig, use_container_width=True)


if __name__ == "__main__":
    main()
