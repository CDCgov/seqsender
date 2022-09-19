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
import pathlib   

# Get working and home directory
work_dir = os.getcwd()
home_dir = os.path.expanduser('~')

def main():

	parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,
		description='Updating the upload log file for submissions')

	parser.add_argument("--log",
		help="Upload log file",
		required=True)

	parser.add_argument("--config",
		help="Config file",
		required=True)

	parser.add_argument("--unique_name",
		help="Config file",
		required=True)

	args = parser.parse_args()
	
	upload_log_file = args.log; config_file = args.config; unique_name = args.unique_name;

	if os.path.exists(upload_log_file):

		# Read in the log file
		upload_log = pd.read_csv(upload_log_file, header = 0, dtype = str, engine = "python", encoding="windows-1254", index_col=False)
		upload_log = upload_log.fillna("")

		# Read in the config file
		with open(config_file, "r") as f:
			config_dict = yaml.safe_load(f)

		# Get the submission directory from config file
		submission_dir = config_dict["general"]["submission_directory"]

		# Change config file to user's working directory
		upload_log.at[upload_log.index[upload_log['name'] == unique_name], 'directory'] = os.path.abspath(re.sub("^~", home_dir, os.path.join(submission_dir, unique_name)))
		upload_log.at[upload_log.index[upload_log['name'] == unique_name], 'config'] = os.path.abspath(re.sub("^~", home_dir, config_file))

		# print(upload_log)

		# Copy upload log back to users
		upload_log.to_csv(os.path.join(work_dir, "upload_log.csv"), header = True, index = False)

	else:

		print("\n" + "Error: upload log file does not exist at:" + upload_log_file + "\n")
		sys.exit(1)
		
if __name__ == "__main__":
    main()
