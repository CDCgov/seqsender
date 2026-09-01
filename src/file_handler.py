#!/usr/bin/env python3

###########################	Description	##################################
# Functions to handle validating file locations and loading files
################################################################################

import sys
import os
from typing import Dict, Any
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import Optional
import shutil
import yaml
import time

from src.settings import SAMPLE_NAME_DATABASE_PREFIX, PROG_DIR

def copy_file(source: str, destination: str):
	shutil.copy(source, destination)

# Validate file exists or error out
def validate_file(file_type: str, file_path: str):
	if not os.path.isfile(file_path):
		print(f"Error: Input {file_type.replace('_',' ')} does not exist at: {file_path}", file=sys.stderr)
		sys.exit(1)

def validate_directory(name: str, path: str):
	if not os.path.exists(path):
		print(f"There is no {name} at: {path}", file=sys.stderr)
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

def save_yaml(config_dict: dict[str, Any], yaml_path: str):
	with open(yaml_path, "w") as file:
		try:
			yaml.safe_dump(config_dict, file, sort_keys = False)
		except:
			print(f"Error: Config file was unable to be written to {yaml_path}. Please check your file permissions.", file=sys.stderr)
			sys.exit(1)

# Is a entire pandas row made of just whitespace, empty strings, or None
def is_row_empty(row: pd.Series) -> bool:
	return all(cell is None or (isinstance(cell, str) and cell.strip() == "") for cell in row)

# Load csv into pandas df and clean then return
def load_csv(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
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
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Save submission xml
def save_xml(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save pandas df to csv
def save_csv(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def create_fasta(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
