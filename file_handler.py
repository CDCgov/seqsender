#!/usr/bin/env python3

###########################	Description	##################################
# Functions to handle validating file locations and loading files
################################################################################

import sys
import os
from typing import Dict
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import List

# Validate file exists or error out
def validate_file(file_type: str, file_path: str):
	if not os.path.isfile(file_path):
		print(f"Error: Input {file_type.replace('_',' ')} does not exist at: {file_path}", file=sys.stderr)
		sys.exit(1)

# Take dictionary of all files to validate
def validate_files(file_dict: Dict[str, str], organism):
	for file_type, file_path in file_dict.items():
		if file_path:
			validate_file(file_type, file_path)

def validate_directory(name: str, path: str):
	if not os.path.exists(path):
		print(f"There is no {name} at: {path}", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def validate_gisaid_installer(submission_dir: str, organism: str):
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	if not os.path.isfile(gisaid_cli_path_option_one) and not os.path.isfile(gisaid_cli_path_option_two):
		print(f"Error: There is not a GISAID CLI for {organism} located at: {gisaid_cli_path_option_one}", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at {gisaid_cli_path_option_one}", file=sys.stderr)
		sys.exit(1)

# Create directory and don't error out if it already exists
def create_directory(path: str):
	os.makedirs(path, exist_ok = True)

# Load yaml file or error out
def load_yaml(yaml_type: str, yaml_path: str):
	with open(yaml_path, "r") as file:
		try:
			config_dict = yaml.load(file, Loader = yaml.FullLoader)
		except:
			print(f"Error: {yaml_type} is incorrect. File must be a valid yaml format.", file=sys.stderr)
			sys.exit(1)
	return config_dict

# Load csv into pandas df and clean then return
def load_csv(file_path: str) -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load fasta file into pandas df and clean then return
def load_fasta_file(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "r") as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in record:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Save submission xml
def save_xml(submission_xml: bytes, submission_files_dir: str) -> None:
	# Save string as submission.xml
	with open(os.path.join(submission_files_dir, "submission.xml"), "wb") as f:
		f.write(submission_xml)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_files_dir, "submission.xml")):
		time.sleep(10)

# Create fasta file based on database
def create_fasta(organism: str, database: List[str], metadata: pd.DataFrame, submission_files_dir: str) -> None:
	records = []
	for index, row in metadata.iterrows():
		records.append(SeqRecord(row["fasta_sequence_orig"], id = row["sequence_name"], description = ""))
	with open(os.path.join(submission_files_dir, "sequence.fsa"), "w+") as f:
		SeqIO.write(records, f, "fasta")
