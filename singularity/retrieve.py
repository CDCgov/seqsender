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
		description='Return a sample file used by seqsender to submit submissions such as config file, metadata file, and required columns file.')

	parser.add_argument("--sample_file",
		help="Name of a file to retrieve: \noptions: {config/metadata/required_columns}",
		required=True)

	parser.add_argument("--outdir",
		help="Output directory to return the file. Default is 'template' directory under current working directory.",
		required=False,
		default="")

	args = parser.parse_args()
	
	file = args.sample_file
	output = args.outdir

	if output is None or output == "":
		output = os.path.join(work_dir, "template")
	else:
		output = re.sub("^~", home_dir, output)
		output = os.path.abspath(output)
		
	# Check if output directory exists	
	if os.path.exists(output) == False:
		os.makedirs(output)	

	if file == "config":

		source = os.path.join(os.path.dirname(os.path.abspath(__file__)), "template", "default_config_template.yaml")
		destination = os.path.join(output, "default_config.yaml")
		shutil.copy(source, destination)

	elif file == "metadata":

		source = os.path.join(os.path.dirname(os.path.abspath(__file__)), "template", "metadata_template.tsv")
		destination = os.path.join(output, "test_metadata.tsv")
		shutil.copy(source, destination)

	elif file == "required_columns":

		source = os.path.join(os.path.dirname(os.path.abspath(__file__)), "template", "required_columns_template.yaml")
		destination = os.path.join(output, "required_columns.yaml")
		shutil.copy(source, destination)

	else:

		print("\n" + "Error: There is no file in the submission pipeline with that name." + "\n")
		sys.exit(1)

if __name__ == "__main__":
    main()
