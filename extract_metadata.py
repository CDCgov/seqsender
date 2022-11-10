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
		description='Making sure metadata file exists for submissions')

	parser.add_argument("--metadata",
		help="Metadata file",
		required=True)

	parser.add_argument("--config_file",
		help="Config file",
		required=True)

	args = parser.parse_args()
	
	metadata_file = args.metadata; config_file = args.config_file;

	# print(metadata_file); print(config_file)

	if os.path.exists(metadata_file) and os.path.exists(config_file):

		# Read in the config file
		with open(config_file, "r") as f:
			config_dict = yaml.safe_load(f)
		
		# Read in the metadata file
		metadata = pd.read_csv(metadata_file, header = 0, dtype = str, sep = config_dict["general"]["metadata_file_sep"], engine = "python", encoding="windows-1254", index_col=False)
		metadata = metadata.fillna("")

		sra_file_path_1 = [os.path.abspath(re.sub("^~", home_dir, file)) for file in metadata['sra_file_path_1']]
		sra_file_path_2 = [os.path.abspath(re.sub("^~", home_dir, file)) for file in metadata['sra_file_path_2']]

		# print(sra_file_path_1); print(sra_file_path_2);

		metadata['sra_file_path_1'] = sra_file_path_1
		metadata['sra_file_path_2'] = sra_file_path_2

		metadata.to_csv(metadata_file, header = True, index = False, sep = "\t")

	else:

		if not os.path.exists(metadata_file):
			print("\n" + "Error: metadata file does not exist at:" + metadata_file + "\n")
			sys.exit(1)

		if not os.path.exists(config_file):
			print("\n" + "Error: config file does not exist at:" + config_file + "\n")
			sys.exit(1)

if __name__ == "__main__":
    main()
