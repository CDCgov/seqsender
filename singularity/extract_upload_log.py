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

		# Config file in pipeline directory
		config_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config_files")

		# Check if config_files directory exists in pipeline directory
		if os.path.exists(config_dir) == False:
			os.makedirs(config_dir)

		# Copy the config files to pipeline directory
		for i in range(len(upload_log)):

			config_file = upload_log.loc[i, "config"]	

			if os.path.exists(config_file):		

				outfile = os.path.join(config_dir, os.path.basename(config_file))

				# Save config file back to user 
				try:
					shutil.copy(config_file, outfile)
				except shutil.SameFileError:
					pass

				# Update config file to docker path
				upload_log.loc[i, 'config'] = outfile

		upload_log.to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "upload_log.csv"), header = True, index = False)

	else:

		print("\n" + "Error: upload log file does not exist at:" + upload_log_file + "\n", file=sys.stderr)

if __name__ == "__main__":
    main()
