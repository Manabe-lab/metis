import streamlit as st
from streamlit_file_browser import st_file_browser
import os

st.header('File browser')
event = st_file_browser( key='A',
	path = '/home/lab/sftp-data/SCALA-data/SCALA-download/Personal_folders',
	show_new_folder = True,
	show_choose_file = False,
	show_download_file = True,
	show_upload_file=True,
	use_static_file_server = False,
	)
if event is not None:
	if event['type'] == 'SELECT_FILE':
		st.write( event)



