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

# print(work_dir)
# print(home_dir)

def main():

	parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,
		description='Making sure config file exists for submissions.')

	parser.add_argument("--config",
		help="Config file",
		required=True)

	args = parser.parse_args()

	config_file = args.config	

	if os.path.exists(config_file):
		
		# Read in the config file
		with open(config_file, "r") as f:
			config_dict = yaml.safe_load(f)

		# Get the submission directory from config file
		submission_dir = config_dict["general"]["submission_directory"]

		# If submission directory is empty, use 'output' as default directory
		if submission_dir is None or submission_dir == "":
			submission_dir = os.path.join(work_dir, "output")
		else:
			submission_dir = re.sub("^~", home_dir, submission_dir)
			submission_dir = os.path.abspath(submission_dir)
		
		# print(submission_dir)

		# Save config file to pipeline directory
		config_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config_files")

		# Check if config_files directory exists in pipeline directory
		if os.path.exists(config_dir) == False:
			os.makedirs(config_dir)		

		# Create a outfile for the pipeline
		outfile = os.path.join(config_dir, os.path.basename(config_file))

		# Save config file to config_files directory
		os.system("bash %s/save_config.sh %s %s %s" % (os.path.dirname(os.path.abspath(__file__)), config_file, submission_dir, outfile))

	else:

		print("\n" + "Error: config file does not exist at:" + config_file + "\n")
		sys.exit(1)

if __name__ == "__main__":
	main()
