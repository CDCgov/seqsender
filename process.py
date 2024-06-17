
# Python Libraries
import pathlib
import pandas as pd
import sys
import os
import time
import re
import xmltodict
import xml.etree.ElementTree as ET
import requests
import yaml
from distutils.util import strtobool
from datetime import datetime
from pandera import pandera, DataFrameSchema, Column, Check, Index, MultiIndex
import ftplib
import json
import importlib
from cerberus import Validator
from typing import List, Set, Dict, Any, Union

# Local imports
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import file_handler
from config.seqsender_schema import schema as seqsender_schema
from config.seqsender_upload_log_schema import schema as seqsender_upload_log_schema

# Get program directory
PROG_DIR: str = os.path.dirname(os.path.abspath(__file__))

# Check the config file
def get_config(config_file: str, database: List[str]) -> Dict[str, Any]:
	# Determine required database
	submission_portals = []
	if "BIOSAMPLE" in database or "SRA" in database or "GENBANK" in database:
		submission_portals.append("ncbi")
	if "GISAID" in database:
		submission_portals.append("gisaid")
	# Check if list empty
	if not submission_portals:
		print("Error: Submission portals list cannot be empty.", file=sys.stderr)
		sys.exit(1)
	submission_schema = "_".join(submission_portals)
	# Read in user config file
	config_dict = file_handler.load_yaml(yaml_type = "Config file", yaml_path = config_file)
	# Check if yaml forms dictionary
	print(config_dict)
	if type(config_dict) is dict:
		schema = eval(open(os.path.join(PROG_DIR, "config", "config_file", (submission_schema + "_schema.py")), 'r').read())
		database_specific_config_updates(schema, database)
		validator = Validator(schema)
		# Validate based on schema
		if validator.validate(config_dict, schema) is False:
			print("Error: Config file is not properly setup. Please correct config file based on issue below:", file=sys.stderr)
			print(validator.errors, file=sys.stderr)
			sys.exit(1)
		else:
			return config_dict["Submission"]
	else:
		print("Error: Config file is incorrect. File must be a valid yaml format.", file=sys.stderr)
		sys.exit(1)

def database_specific_config_updates(schema: Dict[str, Any], database: List[str]) -> Dict[str, Any]:
	# Update seqsender base schema to include needed checks
	if "BIOSAMPLE" in database:
		schema["Submission"]["schema"]["NCBI"]["schema"]["BioSample_Package"]["required"] = True
		schema["Submission"]["schema"]["NCBI"]["schema"]["BioSample_Package"]["nullable"] = False
	return schema

# Read in metadata file
def get_metadata(database: List[str], organism: str, metadata_file: str, config_dict: Dict[str, Any]) -> pd.DataFrame:
	# Read in metadata file
	metadata = file_handler.load_csv(metadata_file)
	# Update seqsender base schema to include needed checks
	if "BioSample" in database or "SRA" in database:
		seqsender_schema.update_columns({"bioproject":{"checks":Check.str_matches(r"^(?!\s*$).+"),"nullable":False,"required":True}})
	biosample_schema = sra_schema = genbank_schema = genbank_cmt_schema = genbank_src_schema = gisaid_schema = None
	# Import schemas
	schemas_dict = dict()
	if "BIOSAMPLE" in database:
		schemas_dict["BioSample"] = importlib.import_module("config.biosample." + config_dict["NCBI"]["BioSample_Package"].strip().replace(".", "_")).schema
	if "SRA" in database:
		schemas_dict["SRA"] = importlib.import_module("config.sra_schema").schema
	if "GENBANK" in database:
		schemas_dict["GenBank"] = importlib.import_module("config.genbank.genbank_schema").schema
		if "cmt-" in metadata.columns:
			schemas_dict["GenBank comment"] = importlib.import_module("config.genbank.genbank_cmt_schema").schema
		if "src-" in metadata.columns:
			schemas_dict["GenBank source"] = importlib.import_module("config.genbank.genbank_src_schema").schema
	if "GISAID" in database:
		schemas_dict["GISAID"] = importlib.import_module("config.gisaid.gisaid_" + organism + "_schema").schema
	# Validate metadata on schema's
	error_msg_list: List[Union[str, pandera.errors.SchemaErrors]] = []
	# Validate required columns for seqsender
	try:
		seqsender_schema.validate(metadata, lazy = True)
	except pandera.errors.SchemaErrors as schema_error:
		error_msg_list.append("Error: Seqsender columns are incorrect:")
		error_msg_list.append(schema_error)
	# Validate required columns for databases
	for schema_name, schema in schemas_dict.items():
		try:
			schema.validate(metadata, lazy = True)
		except pandera.errors.SchemaErrors as schema_error:
			error_msg_list.append(f"Error: {schema_name} columns are incorrect:")
			error_msg_list.append(schema_error)
	if error_msg_list:
		for error_msg in error_msg_list:
			print(error_msg, file=sys.stderr)
		sys.exit(1)
	return metadata

# Check user credentials information
def check_credentials(config_dict: Dict[str, Any], database: str) -> None:
	# Check username
	if "Username" not in config_dict.keys():
		print("Error: there is no Submission > " + database + " > Username information in config file.", file=sys.stderr)
		sys.exit(1)
	elif ("Username" in config_dict.keys()) and ((config_dict["Username"] is not None) and (config_dict["Username"] != "")):
		pass
	else:
		print("Error: Submission > " + database + " > Username in the config file cannot be empty.", file=sys.stderr)
		sys.exit(1)
	# Check password
	if "Password" not in config_dict.keys():
		print("Error: there is no Submission > " + database + " > Password information in config file.", file=sys.stderr)
		sys.exit(1)
	elif ("Password" in config_dict.keys()) and ((config_dict["Password"] is not None) and (config_dict["Password"] != "")):
		pass
	else:
		print("Error: Submission > " + database + " > Password in the config file cannot be empty.", file=sys.stderr)
		sys.exit(1)
	# Check client-id if database is GISAID
	if database != "GISAID":
		return
	elif "Client-Id" not in config_dict.keys():
		print("Error: there is no Submission > " + database + " > Client-Id information in config file.", file=sys.stderr)
		sys.exit(1)
	elif ("Client-Id" in config_dict.keys() and ((config_dict["Client-Id"] is not None) and (config_dict["Client-Id"] != ""))):
		pass
	else:
		print("Error: Submission > " + database + " > Client-Id in the config file cannot be empty.", file=sys.stderr)
		sys.exit(1)

# Update dataframeschema based on today's date
# Uses regex to allow for all valid date formats: YYYY, YYYY-MM, YYYY-MM-DD
def datetime_schema_regex() -> str:
	# January, March, April, May, June, July/August, September, October, November, December
	specific_month_ending_regex = ["[0][1][-][3][0-1])$","[0][3][-][3][0-1])$","[0][4][-][3][0])$","[0][5][-][3][0-1])$","[0][6][-][3][0])$","[0][7-8][-][3][0-1])$","[0][9][-][3][0])$","[1][0][-][3][0-1])$","[1][1][-][3][0])$","[1][2][-][3][0-1])$"]
	# Prefilled to cover 1900's
	datetime_regex = "([1][9]\d{2})$|([1][9]\d{2}[-][0]\d)$|([1][9]\d{2}[-][1][0-2])$|([1][9]\d{2}[-][0]\d[-][0-2]\d)$|([1][9]\d{2}[-][1][0-2][-][0-2]\d)$|([1][9]\d{2}[-][0][1][-][3][0-1])$|([1][9]\d{2}[-][0][3][-][3][0-1])$|([1][9]\d{2}[-][0][4][-][3][0])$|([1][9]\d{2}[-][0][5][-][3][0-1])$|([1][9]\d{2}[-][0][6][-][3][0])$|([1][9]\d{2}[-][0][7-8][-][3][0-1])$|([1][9]\d{2}[-][0][9][-][3][0])$|([1][9]\d{2}[-][1][0][-][3][0-1])$|([1][9]\d{2}[-][1][1][-][3][0])$|([1][9]\d{2}[-][1][2][-][3][0-1])$"
	# Previous years for 2000's
	today = datetime.now()
	year = today.strftime("%Y")
	# Covers just year dates
	datetime_regex = datetime_regex + "|([2][0][0-" + year[2] +"][0-" + year[3] + "])$"
	# Covers previous decades and all their months in correct format
	prev_decade = "([2][0][0-" + str(int(year[2]) - 1) + "]\d"
	all_month_regex = ["[0][1-9])$","[1][0-2])$","[0][1-9][-][0][1-9])$","[0][1-9][-][1-2]\d)$","[1][0-2][-][0][1-9])$","[1][0-2][-][1-2]\d)$"]
	datetime_regex = datetime_regex + "|" + "|".join([(prev_decade + "[-]" + x) for x in all_month_regex])
	datetime_regex = datetime_regex + "|" + "|".join([(prev_decade + "[-]" + x) for x in specific_month_ending_regex])
	# Covers previous individual years and all their months in correct format
	if int(year[3]) > 0:
		prev_year = "([2][0][" + str(int(year[2])) + "][0-" + str(int(year[3]) - 1) + "]"
		datetime_regex = datetime_regex + "|" + "|".join([(prev_year + "[-]" + x) for x in all_month_regex])
		datetime_regex = datetime_regex + "|" + "|".join([(prev_year + "[-]" + x) for x in specific_month_ending_regex])
	# Covers all months and days in current year
	curr_year = "([2][0][" + year[2] + "][" + year[3] + "]"
	month = today.strftime("%m")
	month_list = []
	# Cover all previous months
	if int(month) >= 10:
		month_list.append("[0]\d)$")
		month_list.append("[0]\d[-][0][1-9])$")
		month_list.append("[0]\d[-][1-2]\d)$")
		if int(month[1]) == 0:
			month_list.append("[1][0])$")
		elif int(month[1]) == 1:
			month_list.append("[1][0-1])$")
			month_list.append("[1][0][-][0][1-9])$")
			month_list.append("[1][0][-][1-2]\d)$")
		else:
			month_list.append("[1][0-2])$")
			month_list.append("[1][0-1][-][0][1-9])$")
			month_list.append("[1][0-1][-][1-2]\d)$")
	elif int(month) > 1:
		month_list.append("[0][1-" + str(int(month)) + "])$")
		month_list.append("[0][" + str(int(month) - 1) + "][-][0][1-9])$")
		month_list.append("[0][" + str(int(month) - 1) + "][-][1-2]\d)$")
	else:
		month_list.append("[0][1])$")
	for month_val in specific_month_ending_regex:
		# Getting previous month end date
		if int(month) > int(month_val[4]):
			month_list.append(month_val)
	# Covers current month and day
	day = today.strftime("%d")
	if int(day) >= 10:
		month_list.append("[" + str(int(month[0])) + "][" + str(int(month[1])) + "][-][0][1-9])$")
		if int(day) >= 20:
			month_list.append("[" + str(int(month[0])) + "][" + str(int(month[1])) + "][-][1-" + str(int(day[0]) - 1) + "]\d)$")
		month_list.append("[" + str(int(month[0])) + "][" + str(int(month[1])) + "][-][" + str(int(day[0])) + "][0-" + str(int(day[1])) + "])$")
	else:
		month_list.append("[" + str(int(month[0])) + "][" + str(int(month[1])) + "][-][0][1-" + str(int(day[1])) + "])$")
	datetime_regex = datetime_regex + "|" + "|".join([(curr_year + "[-]" + x) for x in month_list])
	return datetime_regex

# Check sample names in metadata file are listed in fasta file
def process_fasta_samples(metadata: pd.DataFrame, fasta_file: str) -> pd.DataFrame:
	fasta_df = file_handler.load_fasta_file(fasta_file)
	# Check duplicates in fasta_df
	duplicated_df = fasta_df[fasta_df.duplicated(subset = ["fasta_name_orig"], keep = False)]
	if not duplicated_df.empty:
		print("Error: Sequences in fasta file must be unique at: " + fasta_file + "\nDuplicate Sequences\n" + duplicated_df["fasta_sequence_orig"].to_string(index=False), file=sys.stderr)
		sys.exit(1)
	# Validate duplicates don't appear on merge
	try:
		merged_df = metadata.merge(fasta_df, how = "outer", left_on = "sequence_name", right_on = "fasta_name_orig", validate = "1:1")
	except:
		print("Error: Unable to merge fasta file to metadata file. Please validate there are not duplicate sequences in both files.", file=sys.stderr)
		sys.exit(1)
	# Check if fasta file has sequences not in metadata
	if merged_df["sequence_name"].isnull().any():
		print("Error: Sequences in fasta file do not have an associated sequence in metadata file. Please update sequences below:\n" + merged_df[merged_df["sequence_name"].isnull()]["fasta_name_orig"].to_string(), file=sys.stderr)
		sys.exit(1)
	# Check if metadata has sequences not in fasta file
	if merged_df["fasta_name_orig"].isnull().any():
		print("Error: Sequences in metadata file do not have an associated sequence in fasta file. Please update sequences below:\n" + merged_df[merged_df["fasta_name_orig"].isnull()]["sequence_name"].to_string(), file=sys.stderr)
		sys.exit(1)
	return merged_df
