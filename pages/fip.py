import streamlit as st
import requests
from bs4 import BeautifulSoup
from ftplib import FTP
import os
from pathlib import Path
import re
import subprocess
import tarfile
import sys
import requests
import time


def list_dirs(dir_path):
    return [name for name in os.listdir(dir_path) if os.path.isdir(os.path.join(dir_path, name))]


if 'key_count' not in st.session_state: # new directory を入れるst.text_inputのkeyをupdateしないと値がクリアーされず、何度も追加される。
   st.session_state['key_count'] =  1

if 'SCALAorcellxgene' not in st.session_state:
    st.session_state['SCALAorcellxgene'] = "SCALA"

st.markdown("### Download a file for SCALA and cellxgene\n")
SCALAorcellxgene = st.radio(
    "Download a file for SCALA or cellxgene",
    ["SCALA", "cellxgene"], label_visibility = 'collapsed',
    index=0)


if 'dir_path' not in st.session_state:
    st.session_state['dir_path'] = '/home/lab/sftp-data/SCALA-data/SCALA-download/Personal_folders' # specify your
    st.session_state['home_dir'] = '/home/lab/sftp-data/SCALA-data/SCALA-download/Personal_folders'


if SCALAorcellxgene == "SCALA" and st.session_state['SCALAorcellxgene'] == "cellxgene":
    st.session_state['home_dir']  = '/home/lab/sftp-data/SCALA-data/SCALA-download/Personal_folders'
    st.session_state['dir_path'] = st.session_state['home_dir']
    st.session_state['SCALAorcellxgene'] = "SCALA"
elif SCALAorcellxgene == "cellxgene" and st.session_state['SCALAorcellxgene'] == "SCALA":
    st.session_state['home_dir']  = '/home/lab/sftp-data/SCALA-data/SCALA-download/cellxgene-data/Uploaded'
    st.session_state['dir_path'] = st.session_state['home_dir']
    st.session_state['SCALAorcellxgene'] = "cellxgene"


st.markdown("### Select the folder to store the file\n")
st.markdown("##### Now selected directory is " + st.session_state.dir_path.replace(st.session_state['home_dir'] ,''))
dirs = list_dirs(st.session_state['dir_path'])
dirs.sort()
dirs.insert(0, '')
dirs.insert(1, '..')
st.write("##### Select a directory:")
choice = st.selectbox('Select a directory', dirs, key =st.session_state['key_count'],  label_visibility = 'collapsed' )
#if choice == '..':
#    st.session_state['dir_path'] = os.path.dirname(st.session_state['dir_path'])
#elif choice != '..':
if (choice !=  '') and (choice != '..'):
    st.session_state['dir_path'] = os.path.join(st.session_state['dir_path'], choice)
    st.session_state['key_count'] += 2
elif (choice == '..') and (len(st.session_state['dir_path']) > 62):
    st.session_state['dir_path'] = os.path.dirname(st.session_state['dir_path'])
    st.session_state['key_count'] += 2
st.markdown("###### To make a new directory in " + st.session_state.dir_path.replace(st.session_state['home_dir'] ,'') + ':')
new_dir_name = st.text_input('Type a name for the new directory', value= "", key = st.session_state['key_count'] + 1) #keyを指定しないと入力した値がずっと残る。
if st.button("Make a new directory"):
    if new_dir_name:
        st.session_state['key_count'] += 2
        new_dir_path = os.path.join(st.session_state['dir_path'], new_dir_name)
        if not os.path.exists(new_dir_path):
            os.makedirs(new_dir_path)
            st.write(f'New directory `{new_dir_name}` has been created.')
            st.session_state["dir_path"] = new_dir_path
st.write("")
target_dir = st.session_state["dir_path"].replace(st.session_state['home_dir'] ,"")
if st.button(f'Select this folder: {target_dir}', type = 'primary') or st.session_state["dir_path"] != st.session_state['home_dir'] :
    if st.session_state["dir_path"] == st.session_state['home_dir'] :
        sys.exit(1)
    st.markdown("#### Downloaded file will be in " +  st.session_state['SCALAorcellxgene'] + ": " + st.session_state['dir_path'].replace(st.session_state['home_dir'] ,""))

    st.markdown('---\n')

    uploaded_file = st.file_uploader("Choose a file to upload", type=None)


    # Button to trigger file upload and save
    if st.button("Upload") and uploaded_file is not None:
       with open(uploaded_file.name, 'wb') as f:
            f.write(uploaded_file)

