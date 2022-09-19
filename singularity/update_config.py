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

	args = parser.parse_args()
	
	config_file = args.config

	if os.path.exists(config_file):

		# Create a config directory
		config_dir = os.path.join(work_dir, "config_files")

		# Check if config_files directory exists
		if os.path.exists(config_dir) == False:
			os.makedirs(config_dir)	

		# Directory to save config file
		outfile = os.path.join(config_dir, os.path.basename(config_file))

		# Save config file back to user 
		try:
			shutil.copy(config_file, outfile)
		except shutil.SameFileError:
			pass

		# Notify users where config file is stored on user's directory
		print("\n" + "Using config file stored at: " + outfile + "\n")

	else:

		print("\n" + "Error: config file does not exist at:" + config_file + "\n", file=sys.stderr)

if __name__ == "__main__":
	main()
