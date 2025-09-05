import streamlit as st
import re
from collections import Counter

def parse_input(input_text):
    elements = re.split(r'[,;\s\t]+', input_text)
    orig_list = list(filter(None, elements))
    elements_list = list(set(orig_list))
    elements_list.sort()
     # Find duplicates and their counts
    counter = Counter(orig_list)
    duplicates = {item: count for item, count in counter.items() if count > 1}
    
    return orig_list, elements_list, duplicates  # Return original list, unique items, and duplicates


st.markdown("Remove duplicates")

with st.form("Input_items"):
    group_elements = st.text_area(
                        "Input list",
                        help="Enter items separated by space, comma, semicolon, or tab"
                    )
    submitted = st.form_submit_button("Submit")
orig_list, union, duplicates = parse_input(group_elements)
st.write(f"Number of original items:{len(orig_list)}")
st.write(f"Number of unique items:{len(union)}")
if duplicates:
    st.write("Duplicated items:")
    for item, count in duplicates.items():
        st.write(f"- '{item}' appears {count} times")
st.write("Tab separated list:")
st.write("\t".join(union))
st.write("Comma separated list:")
st.write(f"""'{"','".join(union)}'""")
