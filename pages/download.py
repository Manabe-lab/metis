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

def get_all_links(url):
    try:
        response = requests.get(url)
        soup = BeautifulSoup(response.text, 'html.parser')
        links = [a['href'] for a in soup.find_all('a', href=True)]
        rawtar = ''
        if gse + '_RAW.tar' in response.text:
            rawtar = gse + '_RAW.tar'
        return links, rawtar
    except Exception as e:
        st.write(f"An error occurred: {e}")
        return []

def list_ftp_files(path):
    # Connect to the FTP server
    with FTP('ftp.ncbi.nlm.nih.gov') as ftp:
         # List the files
        ftp.login(user="anonymous", passwd="")
        files = ftp.nlst(path)
        return files

def download_ftp_file(path, save_dir, targz_files, tar_files, gz_files):
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login(user="anonymous", passwd="")
    down_file = path.split('/')[-1]
    abs_path = save_dir + '/' + down_file
    with open( abs_path, 'wb') as fp:
        ftp.retrbinary('RETR ' + path, fp.write)
    # gz fileのとき
    if down_file.endswith('.tar.gz'):
        targz_files.append(abs_path)
    elif down_file.endswith('.gz'): # tar.gzではないgz
        if ("barcodes.tsv.gz" not in down_file) and ("features.tsv.gz" not in down_file) and ("matrix.mtx.gz" not in down_file):
            gz_files.append(abs_path)
    if down_file.endswith('.tar'):
        tar_files.append(abs_path)
#        subprocess.call('pigz -d ' + abs_path)
    return targz_files, tar_files, gz_files

def expand_gz(path, ):
    for i in path:
        os.system('pigz -d -f ' + i)

def expand_rawtar(path, save_dir, targz_files, tar_files, gz_files):
    st.write('expand_rawtar')
    st.write(path)
    for i in path:
        # まず内容をみる
        st.write("expanding " + i)
        with tarfile.open(name=i, mode="r") as mytar:
            members = mytar.getmembers()
      #      st.write("tar file contents:")
      #      st.write([x.name for x in members])
        for x in members:
            down_file = x.name
            if down_file.endswith('.tar.gz'):
                targz_files.append(save_dir + "/" + down_file)
            elif down_file.endswith('.gz'): # tar.gzではないgz
                if ("barcodes.tsv.gz" not in down_file) and ("features.tsv.gz" not in down_file) and ("matrix.mtx.gz" not in down_file):
                    gz_files.append(save_dir + "/" + down_file)
            if down_file.endswith('.tar'):
                tar_files.append(save_dir + "/" + down_file)
        os.system('tar -xvf ' + i + ' -C ' + save_dir)
        os.remove(i)
    return targz_files, tar_files, gz_files

def expand_tar(path, save_dir): #gz_filesはglobal
    for i in path:
        # まず内容をみる
        with tarfile.open(name=i, mode="r") as mytar:
            members = mytar.getmembers()
            st.write("tar file contents:")
            st.write([x.name for x in members])
        for x in members:
            tar_content = x.name
            if tar_content.endswith('.gz'):
                tar_content_plain = ''.join(tar_content)
                if ("barcodes.tsv.gz" not in tar_content_plain) and ("features.tsv.gz" not in tar_content_plain) and ("matrix.mtx.gz" not in tar_content_plain):
                    gz_files.append(  save_dir + '/' + tar_content)
        os.system('tar -xvf ' + i + ' -C ' + save_dir)
        os.remove(i)

def expand_targz(path, save_dir):
    for i in path:
        with tarfile.open(name=i, mode="r:gz") as mytar: # gzしてあるのも読めるはず
            members = mytar.getmembers()
            st.write("tar file contents:")
            st.write([x.name for x in members])
        for x in members:
            tar_content = x.name
            if tar_content.endswith('.gz'):
                if ("barcodes.tsv.gz" not in tar_content) and ("features.tsv.gz" not in tar_content) and ("matrix.mtx.gz" not in tar_content):
                    gz_files.append(  save_dir + '/' + tar_content)
        os.system('tar -xvzf ' + i + ' -C ' + save_dir)
        os.remove(i)

def list_dirs(dir_path):
    return [name for name in os.listdir(dir_path) if os.path.isdir(os.path.join(dir_path, name))]


@st.cache_data
def requests_get(url):
    response = requests.get(url, verify=False)
    return response

@st.cache_data
def get_retry(url, retry_times=3, errs=[500, 502, 503]):
    for t in range(retry_times + 1):
        r = requests.get(url, verify=False)
        if t < retry_times:
            if r.status_code in errs:
                time.sleep(2)
                continue
        return r


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

    GSEorURL = st.radio(
    "Download a file in GSE/GSM or URL-specified file",
    ["GSE/GSM", "URL"],
    index=0)


    if GSEorURL == 'URL':
        st.markdown("### Target file URL\n")
        url = st.text_input('Enter target file URL (link)', '')
        st.write("You can copy the link by right-clicking on the link on the Web page.")
        if url:
            file_name = os.path.basename(url)
            st.write("Will download: " + file_name )
            if st.button('Start downloading'):
                st.write("Downloading...")
                response = get_retry(url) # requests_get(url)
                saveFilePath = os.path.join(st.session_state['dir_path'], file_name)
                with st.spinner('Writing to disk...'):
                    try:
                        with open(saveFilePath, 'wb') as saveFile:
                            saveFile.write(response.content)
                    except:
                        st.write("Error in downloading the file")
                    else:
                        st.write("Done downloading!")



    if GSEorURL == 'GSE/GSM':
        st.markdown("### GSE/GSM file download\n")
        gse = st.text_input('Enter GSE/GSM (GSE123456)', '')

        if gse:
            url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + gse
            links, rawtar = get_all_links(url)
            ftp_links = []
            if rawtar != "":
                ftp_links.append(rawtar)
            for i in links:
                if "ftp://" in i:
                    ftp_links.append(i)

            st.write("Files:")

            s = ''
            for i in ftp_links:
                s += "- " + i + "\n"
            st.markdown(s)

            ftp_choice = st.selectbox(
            'Choose file or directory', ftp_links)

            # Streamlit application
            st.markdown("### Download file\n")

            if ftp_choice != rawtar:
                try:
                    ftp_path = ftp_choice.replace("ftp://ftp.ncbi.nlm.nih.gov/","") # ftpのpathを除き
                    ftp_path = re.sub(r'\/[^\/]+$','',ftp_path)#最後のpathも除く
                    ftp_file_list = list_ftp_files(ftp_path)
                    st.write("Files to download")
                    s = ''
                    for i in ftp_file_list:
                        s += "- " + i + "\n"
                    st.markdown(s)
                except Exception as e:
                    st.error(f'Error getting file list: {e}')
            else:
                st.write("Download " + rawtar)

            if st.button('Start downloading'):
                if  ftp_choice != rawtar:
                    try:
                        targz_files = []
                        gz_files = []
                        tar_files = []
                        for i in ftp_file_list:
                            targz_files, tar_files, gz_files = download_ftp_file(i, st.session_state['dir_path'], targz_files, tar_files, gz_files)
                        if len(targz_files) > 0:
                            st.write('targz file:')
                            st.write(", ".join([x.split('/')[-1] for x in targz_files]))
                            expand_targz(targz_files, st.session_state['dir_path'])
                        if len(tar_files) > 0:
                            st.write('tar file:')
                            st.write(", ".join([x.split('/')[-1] for x in tar_files]))
                            expand_tar(tar_files, st.session_state['dir_path'])
                        if len(gz_files) > 0:
                            st.write('gz file:')
                            st.write(", ".join([x.split('/')[-1] for x in gz_files]))
                            expand_gz(gz_files)
                        st.success('File downloaded successfully!')
                    except Exception as e:
                        st.error(f'Error downloading file: {e}')
                else:
                    addr = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=" + gse + "&format=file"
                    st.write("Download link:")
                    st.write(addr)
                    try:
                        urlData = requests.get(addr).content
                        with open( st.session_state['dir_path'] + rawtar ,mode='wb') as f: # wb でバイト型を書き込める
                            f.write(urlData)
                        rawtar_files = [st.session_state['dir_path'] + rawtar]
                        st.write('tar_files')
                        st.write(rawtar_files)
                        targz_files = []
                        gz_files = []
                        tar_files = []
                        targz_files, tar_files, gz_files = expand_rawtar(rawtar_files, st.session_state['dir_path'], targz_files, tar_files, gz_files)
                        if len(targz_files) > 0:
                            st.write('targz file:')
                            st.write(", ".join([x.split('/')[-1] for x in targz_files]))
                            expand_targz(targz_files, st.session_state['dir_path'])
                        if len(tar_files) > 0:
                            st.write('tar file:')
                            st.write(", ".join([x.split('/')[-1] for x in tar_files]))
                            expand_tar(tar_files, st.session_state['dir_path'])
                        if len(gz_files) > 0:
                            st.write('gz file:')
                            st.write(", ".join([x.split('/')[-1] for x in gz_files]))
                            expand_gz(gz_files)
                        st.success('File downloaded successfully!')
                    except Exception as e:
                        st.error(f'Error downloading file: {e}')