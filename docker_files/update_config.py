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
		description='Making sure config file exists for submissions.')

	parser.add_argument("--config",
		help="Config file",
		required=True)

	parser.add_argument("--unique_name",
		help="Unique identififer",
		required=False)

	args = parser.parse_args()
	
	config_file = args.config; unique_name = args.unique_name;

	if os.path.exists(config_file):

		# Read in the config file
		with open(config_file, "r") as f:
			config_dict = yaml.safe_load(f)

		# Get the submission directory from config file
		submission_dir = config_dict["general"]["submission_directory"]

		# If submission directory is empty, use default output directory
		submission_directory = re.sub(work_dir, "./", submission_dir)
		submission_directory = re.sub("//*", "/", submission_directory)

		# Create a config directory
		config_dir = os.path.join(work_dir, "config_files")

		# Check if config_files directory exists
		if os.path.exists(config_dir) == False:
			os.makedirs(config_dir)	

		# Directory to save config file
		outfile = os.path.join(config_dir, os.path.basename(config_file))

		# Save config file
		os.system("bash %s/save_config.sh %s %s %s" % (os.path.dirname(os.path.abspath(__file__)), config_file, submission_directory, outfile))

		# Notify users where config file is stored on user's directory
		print("\n" + "Using config file stored at: " + "./config_files/" + os.path.basename(config_file) + "\n")

		# Check the sra path in the output directory
		sra_file_path = os.path.join(submission_dir, unique_name, "biosample_sra", "sra_file_path.txt")

		if os.path.exists(sra_file_path):
			os.system('cat %s | sed "s,%s,.,g" > %s' % (sra_file_path, work_dir, sra_file_path))


if __name__ == "__main__":
	main()
