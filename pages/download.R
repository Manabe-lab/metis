import streamlit as st
import requests
from bs4 import BeautifulSoup

def get_all_links(url):
    try:
        response = requests.get(url)
        soup = BeautifulSoup(response.text, 'html.parser')
        links = [a['href'] for a in soup.find_all('a', href=True)]
        return links
    except Exception as e:
        st.write(f"An error occurred: {e}")
        return []

url = st.text_input('Enter a URL', '')
if url:
    links = get_all_links(url)
    st.write(links)
