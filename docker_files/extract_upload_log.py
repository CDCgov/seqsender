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
		description='Extracting the upload log file for submissions')

	parser.add_argument("--log",
		help="Upload log file",
		required=True)

	args = parser.parse_args()
	
	upload_log_file = args.log

	if os.path.exists(upload_log_file):

		# Read in the upload log file
		upload_log = pd.read_csv(upload_log_file, header = 0, dtype = str, engine = "python", encoding="windows-1254", index_col=False)
		upload_log = upload_log.fillna("")

		# Extract output directory
		upload_log['directory'] = [re.sub("^./", "", name) for name in upload_log['directory']]
		upload_log['directory'] = [re.sub("//*", "/", os.path.join(work_dir, name)) for name in upload_log['directory']]

		# Extract config directory
		upload_log['config'] = [re.sub("^./", "", name) for name in upload_log['config']]
		upload_log['config'] = [re.sub("//*", "/",  os.path.join(work_dir, name)) for name in upload_log['config']]

		# Copy the config files to pipeline directory
		for i in range(len(upload_log)):

			config_file = upload_log.loc[i, "config"]	

			if os.path.exists(config_file):
				
				with open(config_file, "r") as f:
					config_dict = yaml.safe_load(f)

				# Get the submission directory from config file
				submission_dir = config_dict["general"]["submission_directory"]

				# If submission directory is empty, use default output directory
				if submission_dir == "":
					submission_directory = os.path.join(work_dir, "output")
				else:
					submission_directory = re.sub("^./", "", submission_dir)
					submission_directory = re.sub("//*", "/", os.path.join(work_dir, submission_directory))

				# Directory to save config file
				outfile = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config_files", os.path.basename(config_file))

				# Save config file
				os.system("bash %s/save_config.sh %s %s %s" % (os.path.dirname(os.path.abspath(__file__)), config_file, submission_directory, outfile))

				# Update config file to docker path
				upload_log.loc[i, 'config'] = outfile
			else:
				print("Error: The config file path provided in upload_log.csv for " + upload_log.loc[i, 'name'] + " does not exist at: \n" + "./" + upload_log.iloc[i]['config'], file=sys.stderr)
				sys.exit(1)

		upload_log.to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "upload_log.csv"), header = True, index = False)

if __name__ == "__main__":
    main()
