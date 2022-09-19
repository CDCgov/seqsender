#!/usr/bin/env python3

import ftplib
import argparse
from argparse import RawTextHelpFormatter
import sys
from datetime import datetime
import requests
import json
import os
import pandas as pd
import yaml
from Bio import SeqIO
import xml.etree.ElementTree as ET
import shutil
import re

work_dir = os.getcwd()

def main():
	parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,
		description='Updating the upload log file for submissions')

	parser.add_argument("--log",
		help="Upload log file",
		required=True)

	args = parser.parse_args()
	
	upload_log_file = args.log

	upload_log = pd.read_csv(upload_log_file, header = 0, dtype = str, engine = "python", encoding="windows-1254", index_col=False)
	upload_log = upload_log.fillna("")

	upload_log['directory'] = [re.sub(work_dir, "./", name) for name in upload_log['directory']]
	upload_log['directory'] = [re.sub("//*", "/", name) for name in upload_log['directory']]

	upload_log['config'] = [re.sub(os.path.dirname(os.path.abspath(__file__)), "./", name) for name in upload_log['config']]
	upload_log['config'] = [re.sub("//*", "/", name) for name in upload_log['config']]

	upload_log.to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "upload_log.csv"), header = True, index = False)

if __name__ == "__main__":
    main()
