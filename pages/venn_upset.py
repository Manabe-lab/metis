import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3, venn2_circles, venn3_circles
from upsetplot import from_contents, UpSet, plot
import io
import re
import tempfile

def parse_group_data(file_path, direction='row'):
    if direction not in ['row', 'column']:
        raise ValueError("Direction must be either 'row' or 'column'")
    
    result = {}

        # ファイルの内容を行ごとに分割
    lines = file_content.decode('utf-8').splitlines()
    
    lines = [line.strip().split() for line in lines if line.strip()]
    
    if direction == 'row':
        for line in lines:
            if line:
                group = line[0]
                items = line[1:]
                if items:
                    result[group] = items
    else:  # column direction
        # Transpose the data
        transposed = list(map(list, zip(*lines)))
        for column in transposed:
            if column:
                group = column[0]
                items = column[1:]
                if items:
                    result[group] = items
    
    return result


def parse_input(input_text):
    elements = re.split(r'[,;\s\t]+', input_text)
    return list(set(filter(None, elements)))  # Remove duplicates and empty strings


def create_venn_diagram(data, title, font_size, colors, plot_width, plot_height):
    fig, ax = plt.subplots(figsize=(plot_width, plot_height))
    
    if len(data) == 2:
        set1 = set(data[list(data.keys())[0]])
        set2 = set(data[list(data.keys())[1]])
        intersection = set1 & set2
        only_set1 = set1 - intersection
        only_set2 = set2 - intersection
        
        if venn_line:
            v = venn2_circles(subsets=(len(only_set1), len(only_set2), len(intersection)), ax=ax)
        else:
            v = venn2(subsets=(len(only_set1), len(only_set2), len(intersection)),
                      set_labels=data.keys(), ax=ax, set_colors=colors[:2], alpha=0.5)
            
            # Hide labels for regions with zero elements
            if len(only_set1) == 0 and hasattr(v, 'get_label_by_id') and v.get_label_by_id('10') is not None:
                v.get_label_by_id('10').set_text('')
            if len(only_set2) == 0 and hasattr(v, 'get_label_by_id') and v.get_label_by_id('01') is not None:
                v.get_label_by_id('01').set_text('')
            if len(intersection) == 0 and hasattr(v, 'get_label_by_id') and v.get_label_by_id('11') is not None:
                v.get_label_by_id('11').set_text('')
        
        element_info = {
            f"Only in {list(data.keys())[0]}": only_set1,
            f"Only in {list(data.keys())[1]}": only_set2,
            "Intersection": intersection
        }
    
    elif len(data) == 3:
        set1 = set(data[list(data.keys())[0]])
        set2 = set(data[list(data.keys())[1]])
        set3 = set(data[list(data.keys())[2]])
        
        # Calculate proper subset sizes for each region
        only_set1 = set1 - (set2 | set3)  # In set1 only
        only_set2 = set2 - (set1 | set3)  # In set2 only 
        only_set3 = set3 - (set1 | set2)  # In set3 only
        only_set1_set2 = (set1 & set2) - set3  # In set1 and set2 only
        only_set1_set3 = (set1 & set3) - set2  # In set1 and set3 only
        only_set2_set3 = (set2 & set3) - set1  # In set2 and set3 only
        all_sets = set1 & set2 & set3  # In all three sets
        
        if venn_line:
            v = venn3_circles(subsets=(len(only_set1), len(only_set2), len(only_set1_set2),
                          len(only_set3), len(only_set1_set3),
                          len(only_set2_set3),
                          len(all_sets)),
                   ax=ax)
        else:
            v = venn3(subsets=(len(only_set1), len(only_set2), len(only_set1_set2),
                          len(only_set3), len(only_set1_set3),
                          len(only_set2_set3),
                          len(all_sets)),
                  set_labels=data.keys(), ax=ax, set_colors=colors, alpha=0.5)
            
            # Hide labels for regions with zero elements
            if hasattr(v, 'get_label_by_id'):
                region_ids = ['100', '010', '110', '001', '101', '011', '111']
                region_counts = [len(only_set1), len(only_set2), len(only_set1_set2), 
                                len(only_set3), len(only_set1_set3), len(only_set2_set3), 
                                len(all_sets)]
                
                for region_id, count in zip(region_ids, region_counts):
                    if count == 0 and v.get_label_by_id(region_id) is not None:
                        v.get_label_by_id(region_id).set_text('')
        
        element_info = {
            f"Only in {list(data.keys())[0]}": only_set1,
            f"Only in {list(data.keys())[1]}": only_set2,
            f"Only in {list(data.keys())[2]}": only_set3,
            f"In {list(data.keys())[0]} and {list(data.keys())[1]}": only_set1_set2,
            f"In {list(data.keys())[0]} and {list(data.keys())[2]}": only_set1_set3,
            f"In {list(data.keys())[1]} and {list(data.keys())[2]}": only_set2_set3,
            "In all three sets": all_sets
        }
    
    if v is not None and not venn_line:
        for text in v.set_labels:
            if text is not None:
                text.set_fontsize(font_size)
    
    plt.title(title, fontsize=font_size+2)
    return fig, element_info

def create_upsetplot(data, title, font_size, plot_width, plot_height):    
    try:
        data_upset = from_contents(data)
        
        fig = plt.figure(figsize=(plot_width, plot_height), dpi=100)
        
        plot_kwargs = {
            'sort_by': sort_by,
            'orientation': orientation,
            'sort_categories_by': sort_categories_by,
            'element_size': None
        }
        
        if show_counts:
            plot_kwargs['show_counts'] = True
        if show_percentages:
            plot_kwargs['show_percentages'] = True
        
        plot(data_upset, fig=fig, **plot_kwargs)
        
        plt.suptitle(title, fontsize=font_size+4, y=1.02)
            
        return fig
    except Exception as e:
        st.error(f"Error in create_upsetplot: {str(e)}")
        return None

st.sidebar.title("Options")
# Initialize session state
if 'data' not in st.session_state:
    st.session_state.data = {}
if 'num_groups' not in st.session_state:
    st.session_state.num_groups = 2
if 'input_method' not in st.session_state:
    st.session_state.input_method = "Direct Input"
if 'previous_num_groups' not in st.session_state:
    st.session_state.previous_num_groups = st.session_state.num_groups

st.markdown('## Upset plot and Venn diagram generator')

# サイドバーの設定
with st.sidebar:
    st.markdown("#### Plot Settings")
    plot_title = st.text_input("Plot Title", "Set Intersection")
    plot_width = st.number_input("Plot Width (inches)", min_value=2, max_value=20, value=8, step=1)
    plot_height = st.number_input("Plot Height (inches)", min_value=2, max_value=15, value=6, step=1)
    output_format = st.selectbox("Output Format", ["PNG", "PDF"], index=1)
    st.markdown("#### Venn diagram options:")
    font_size = st.number_input("Font Size", min_value=4, max_value=20, value=12, step=1)
    colors = st.color_picker("Color 1", "#3366cc"), st.color_picker("Color 2", "#dc3912"), st.color_picker("Color 3", "#ff9900")
    venn_line = st.checkbox("Line diagram?", value=False)
    st.markdown("#### Upset plot options:")
    sort_by = st.selectbox("How to sort:", ["cardinality", "degree", "input", "-cardinality", "-degree", "-input"], index=0) 
    sort_categories_by = st.selectbox("How to sort groups:", ['cardinality', 'input', '-cardinality', '-input'], index=0)
    orientation = st.selectbox("Orientation:", ['horizontal', 'vertical'], index=0)
    show_counts = st.checkbox("Show counts?", value=True)
    show_percentages = st.checkbox("Show percentages?", value=False)
    st.markdown("#### Clear cache")
    if st.button('Clear cache data and renew'):
      #  st.cache_data.clear()
        st.session_state.clear()
        st.success('キャッシュがクリアされました')
        st.rerun()

# 入力方法の選択
st.session_state.input_method = st.radio("Choose input method:", ("Upload CSV", "Direct Input"))

if st.session_state.input_method == "Upload CSV":
    direct = st.radio("Direction of data:", ["row","column"] )
    if direct == "row":
        st.markdown("Row-wise:")
        st.markdown("G1\ta\tb\tc\td\te")
        st.markdown("G2\td\te\tc")
    else:
        st.markdown("Column-wise:")
        st.markdown("""
| G1 | G2 | G3 |
| --- | --- | --- |
| a | c | e |
| b | d | f |
| d |  | g |""")
    st.write("Either tab or comma separated.")
    uploaded_file = st.file_uploader("Choose a CSV file", type=["csv",'tsv','txt'])
    if uploaded_file is not None:
        # ファイルの内容を読み込む
        file_content = uploaded_file.read()
        st.session_state.data = parse_group_data(uploaded_file, direction=direct)
else:
    new_num_groups = st.number_input("Number of groups", min_value=2, max_value=20, value=st.session_state.num_groups, step=1)
    
    # グループ数が変更された場合、データをリセット
    if new_num_groups != st.session_state.previous_num_groups:
        st.session_state.data = {}
        st.session_state.num_groups = new_num_groups
        st.session_state.previous_num_groups = new_num_groups
        st.rerun()  # 変更を反映するためにページを再読み込み

    with st.form(key='input_form'):
        new_data = {}
        for i in range(st.session_state.num_groups):
            col1, col2 = st.columns([1, 3])
            with col1:
                group_name = st.text_input(f"Group {i+1} name", key=f"group_name_{i}", 
                                           value=list(st.session_state.data.keys())[i] if i < len(st.session_state.data) else "")
            with col2:
                if group_name:
                    current_elements = ", ".join(map(str, st.session_state.data.get(group_name, [])))
                    group_elements = st.text_area(
                        f"Elements for {group_name}",
                        key=f"group_elements_{i}",
                        value=current_elements,
                        help="Enter elements separated by space, comma, semicolon, or tab"
                    )
                    new_data[group_name] = parse_input(group_elements)
        
        submit_button = st.form_submit_button(label='Update Data')
        
    if submit_button:
        st.session_state.data = new_data
        st.success("Data updated successfully!")

# データプレビューの表示
st.write("Current Data:")
for group, elements in st.session_state.data.items():
    st.write(f"{group}: {', '.join(elements)}")




# プロット生成ボタン
if st.button("Generate Plot"):

    if st.session_state.data:
        num_groups = len(st.session_state.data)

        if num_groups <= 3:
            st.write("Generating Scaled Venn Diagram...")
            try:
                fig, element_info = create_venn_diagram(st.session_state.data, plot_title, font_size, colors, plot_width, plot_height)
                st.pyplot(fig)
                
                # 要素情報の表示
                st.write("Set Elements:")
                for set_name, elements in element_info.items():
                    if elements:  # 空でない場合のみ表示
                        st.write(f"{set_name}: {', '.join(sorted(elements))}")
            except Exception as e:
                st.error(f"Error generating Venn diagram: {str(e)}")
                st.error(f"Traceback: {traceback.format_exc()}")
        else:
            st.write("Generating UpsetPlot...")
            
            fig = create_upsetplot(st.session_state.data, plot_title, font_size, plot_width, plot_height)
            if fig:
                st.pyplot(fig, use_container_width=True)

        if fig:
            # プロットのダウンロードボタンを追加
            buf = io.BytesIO()
            if output_format == "PNG":
                fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
                file_extension = "png"
                mime = "image/png"
            else:
                fig.savefig(buf, format='pdf', bbox_inches='tight')
                file_extension = "pdf"
                mime = "application/pdf"
            buf.seek(0)
            st.download_button(
                label=f"Download Figure as {output_format}",
                data=buf,
                file_name=f"figure.{file_extension}",
                mime=mime
            )
    else:
        st.write("Please input data to generate the plot.")
