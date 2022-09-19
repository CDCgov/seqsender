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
		description='Making sure metadata file exists for submissions')

	parser.add_argument("--metadata",
		help="Metadata file",
		required=True)

	args = parser.parse_args()
	
	metadata_file = args.metadata

	if os.path.exists(metadata_file):
		
		# Read in the metadata file
		metadata = pd.read_csv(metadata_file, header = 0, dtype = str, sep = "\t", engine = "python", encoding="windows-1254", index_col=False)
		metadata = metadata.fillna("")

		sra_file_path_1 = [re.sub("^./", "", name) for name in metadata['sra_file_path_1']]
		sra_file_path_1 = [re.sub("//*", "/", os.path.join(work_dir, name)) for name in sra_file_path_1]
		sra_file_path_2 = [re.sub("^./", "", name) for name in metadata['sra_file_path_2']]
		sra_file_path_2 = [re.sub("//*", "/", os.path.join(work_dir, name)) for name in sra_file_path_2]

		for i in range(len(sra_file_path_1)):		
			if os.path.exists(sra_file_path_1[i]) == False:
				print("Error: sra_file_path_1 provided in metadata file for " + metadata.iloc[i]['fasta_name'] + " does not exist at: \n" + "./" + metadata.iloc[i]['sra_file_path_1'], file=sys.stderr)
				sys.exit(1)

			if os.path.exists(sra_file_path_2[i]) == False:
				print("Error: sra_file_path_2 provided in metadata file for " + metadata.iloc[i]['fasta_name'] + " does not exist at: \n" + "./" + metadata.iloc[i]['sra_file_path_2'], file=sys.stderr)
				sys.exit(1)

		metadata['sra_file_path_1'] = sra_file_path_1
		metadata['sra_file_path_2'] = sra_file_path_2

		# Save metadata to pipeline directory
		metadata_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

		if os.path.exists(metadata_dir) == False:
			os.makedirs(metadata_dir) 

		metadata.to_csv(os.path.join(metadata_dir, "metadata.tsv"), header = True, index = False, sep = "\t")

if __name__ == "__main__":
    main()
