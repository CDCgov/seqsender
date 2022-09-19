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
		description='Return a sample file used by seqsender to submit submissions such as config file, metadata file, and required columns file.')

	parser.add_argument("--sample_file",
		help="Name of a file to retrieve: config/metadata/required_columns",
		required=True)

	parser.add_argument("--outdir",
		help="Output directory to return the file",
		required=False,
		default="")

	args = parser.parse_args()
	
	file = args.sample_file
	output = args.outdir

	if output is None or output == "":
		fpath = work_dir
	else:
		fpath = re.sub('^./', "", output)
		fpath = re.sub('//*', "/", os.path.join(work_dir, fpath))
	
	if os.path.exists(fpath) == False:
		os.makedirs(fpath)

	if file == "config":

		source = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config_files", "default_config.yaml")
		destination = os.path.join(fpath, "default_config.yaml")
		shutil.copy(source, destination)

	elif file == "metadata":

		source = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_input", "test_metadata.tsv")
		destination = os.path.join(fpath, "test_metadata.tsv")
		shutil.copy(source, destination)

	elif file == "required_columns":

		source = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config_files", "required_columns.yaml")
		destination = os.path.join(fpath, "required_columns.yaml")
		shutil.copy(source, destination)

	else:

		print("There is no file in the submission pipeline with this name.")

if __name__ == "__main__":
    main()
