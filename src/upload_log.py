#!/usr/bin/env python3

###########################	Description	##################################
# Functions for handling update log
################################################################################

import pandas as pd
from typing import Optional, Any, Union
from pandera import pandera, DataFrameSchema, Column, Check, Index, MultiIndex
import os
import sys
from datetime import datetime
from cryptography.fernet import InvalidToken

import src.genbank_handler as genbank_handler
import src.biosample_sra_handler as biosample_sra_handler
import src.ncbi_handler as ncbi_handler
import src.tools as tools
import src.file_handler as file_handler
from config.seqsender.upload_log_schema import schema as upload_schema
from config.seqsender.submission_status_report.biosample_submission_status_report_schema import schema as status_report_bs_schema
from config.seqsender.submission_status_report.sra_submission_status_report_schema import schema as status_report_sra_schema
from config.seqsender.submission_status_report.genbank_submission_status_report_schema import schema as status_report_gb_schema

from src.settings import SAMPLE_NAME_DATABASE_PREFIX, BIOSAMPLE_SUBMISSION_STATUS_COLUMNS, SRA_SUBMISSION_STATUS_COLUMNS, GENBANK_SUBMISSION_STATUS_COLUMNS, SUBMISSION_LOG_COLUMNS

# create new submission_status.csv based on databases submitting to
def create_submission_status_csv(database: list[str], metadata: pd.DataFrame, submission_dir: str) -> None:
	submission_status_file = os.path.join(submission_dir, "submission_status_report.csv")
	database_columns: list[str]  = []
	sample_name_columns: list[str] = []
	ordered_database_columns: list[str]  = []
	if "BIOSAMPLE" in database:
		database_columns += BIOSAMPLE_SUBMISSION_STATUS_COLUMNS
		sample_name_columns += [f"{SAMPLE_NAME_DATABASE_PREFIX['BIOSAMPLE']}sample_name"]
		ordered_database_columns += [f"{SAMPLE_NAME_DATABASE_PREFIX['BIOSAMPLE']}sample_name"] + BIOSAMPLE_SUBMISSION_STATUS_COLUMNS
	if "SRA" in database:
		database_columns += SRA_SUBMISSION_STATUS_COLUMNS
		sample_name_columns += [f"{SAMPLE_NAME_DATABASE_PREFIX['SRA']}sample_name"]
		ordered_database_columns += [f"{SAMPLE_NAME_DATABASE_PREFIX['SRA']}sample_name"] + SRA_SUBMISSION_STATUS_COLUMNS
	if "GENBANK" in database:
		database_columns += GENBANK_SUBMISSION_STATUS_COLUMNS
		sample_name_columns += [f"{SAMPLE_NAME_DATABASE_PREFIX['GENBANK']}sample_name"]
		ordered_database_columns += [f"{SAMPLE_NAME_DATABASE_PREFIX['GENBANK']}sample_name"] + GENBANK_SUBMISSION_STATUS_COLUMNS
	sample_name_df = metadata[sample_name_columns].copy()
	sample_name_df = sample_name_df.assign(**{col: "" for col in database_columns})
	sample_name_df = sample_name_df.reindex(columns = ordered_database_columns)
	file_handler.save_csv(df=sample_name_df, file_path=submission_status_file)

# Validate data in submission_status.csv file is correctly formatted
def validate_submission_status_df(metadata: pd.DataFrame, database: list[str]) -> None:
	error_msg_list: list[pandera.errors.SchemaErrors] = []
	if "BIOSAMPLE" in database:
		try:
			status_report_bs_schema.validate(metadata, lazy = True)
		except pandera.errors.SchemaErrors as schema_error:
			error_msg_list.append(schema_error)
	if "SRA" in database:
		try:
			status_report_sra_schema.validate(metadata, lazy = True)
		except pandera.errors.SchemaErrors as schema_error:
			error_msg_list.append(schema_error)
	if "GENBANK" in database:
		try:
			status_report_gb_schema.validate(metadata, lazy = True)
		except pandera.errors.SchemaErrors as schema_error:
			error_msg_list.append(schema_error)
	if error_msg_list:
		print(metadata.head())
		tools.pretty_print_pandera_errors(file="submission_status_report.csv", error_msgs=error_msg_list)
		sys.exit(1)

# Update values of existing submission_status.csv file
def update_submission_status_csv(submission_dir: str, update_database: str, update_df: pd.DataFrame) -> None:
	# Check if there are updates to be made
	if update_df.empty:
		print(f"Error: Unable to update 'submission_status.csv' for '{update_database}' at '{submission_dir}'. The log file may be empty.", file=sys.stderr)
	# Pop off directory if inside database directory
	if os.path.split(submission_dir)[-1] in ["BIOSAMPLE", "SRA", "GENBANK"]:
		submission_dir = os.path.dirname(submission_dir)
	submission_status_file = os.path.join(submission_dir, "submission_status_report.csv")
	file_handler.validate_file(file_type="submission status report", file_path=submission_status_file)
	df = file_handler.load_csv(file_path=submission_status_file)
	# Identify database schemas based on <database_prefix>_sample_name
	all_databases = [database_name for database_name, database_prefix in SAMPLE_NAME_DATABASE_PREFIX.items() if (database_prefix + "sample_name") in df.columns]
	validate_submission_status_df(metadata=df, database=all_databases)
	if f"{SAMPLE_NAME_DATABASE_PREFIX[update_database]}sample_name" in update_df:
		sample_column = (SAMPLE_NAME_DATABASE_PREFIX[update_database] + "sample_name")
	elif f"{SAMPLE_NAME_DATABASE_PREFIX[update_database]}segment_name" in update_df:
		sample_column = (SAMPLE_NAME_DATABASE_PREFIX[update_database] + "segment_name")
	df = df.set_index(sample_column, drop = False)
	update_df = update_df.set_index(sample_column)
	df.update(update_df)
	df = df.reset_index(drop = True)
	validate_submission_status_df(metadata=df, database=all_databases)
	file_handler.save_csv(df=df, file_path=submission_status_file)

# Create new row in submission_log.csv for database submission
def create_submission_log(database: str, organism: str, submission_name: str, submission_dir: str, database_dir: str, config_file: str, submission_status: str, submission_id: str, submission_type: str) -> None:
	# if log exists load existing one
	if os.path.isfile(os.path.join(submission_dir, "submission_log.csv")):
		df = load_submission_log(submission_dir=submission_dir)
	else:
		# If log doesn't exist create it
		df = pd.DataFrame(columns = SUBMISSION_LOG_COLUMNS)
	# Create new field
	new_entry = {'Submission_Name': submission_name,
				 'Organism': organism,
				 'Database': database,
				 'Submission_Type': submission_type,
				 'Submission_Date': datetime.now().strftime("%Y-%m-%d"),
				 'Submission_ID': submission_id,
				 'Submission_Status': submission_status,
				 'Submission_Directory': database_dir,
				 'Config_File': config_file,
				 'Update_Date': datetime.now().strftime("%Y-%m-%d")
				}
	df.loc[len(df)] = new_entry # type: ignore
	# Remove duplicates and keep latest update
	df = df.drop_duplicates(subset = ["Submission_Name", "Organism", "Database", "Submission_Type", "Config_File"], keep = "last", ignore_index = True)
	file_handler.save_csv(df=df, file_path=submission_dir, file_name="submission_log.csv")

# Update values of existing submission_log.csv file
def update_submission_log(database: str, organism: str, submission_name: str, submission_log_dir: str, submission_dir: str, submission_status: str, submission_id: str, submission_type: str) -> None:
	df = load_submission_log(submission_dir=submission_log_dir)
	# Select correct fields
	df_partial = df.loc[(df["Organism"] == organism) & (df["Database"] == database) & (df["Submission_Directory"] == submission_dir) & (df["Submission_Name"] == submission_name) & (df["Submission_Type"] == submission_type)]
	# Update existing field
	if df_partial.shape[0] > 0:
		df.loc[df_partial.index.values, "Submission_ID"] = submission_id
		df.loc[df_partial.index.values, "Submission_Status"] = submission_status
		df.loc[df_partial.index.values, "Update_Date"] = datetime.now().strftime("%Y-%m-%d")
	else:
		print(f"Error: '{submission_name}' '{database}' is not present in the submission log at '{submission_log_dir}'.", file=sys.stderr)
		sys.exit(1)
	file_handler.save_csv(df=df, file_path=submission_log_dir, file_name="submission_log.csv")

# Validate submission_dir and log exists and then validate syntax is correct
def load_submission_log(submission_dir: str) -> pd.DataFrame:
	submission_log_file = os.path.join(submission_dir, "submission_log.csv")
	file_handler.validate_file(file_type = "submission_log", file_path = submission_log_file)
	df = file_handler.load_csv(file_path = submission_log_file)
	# Drop no longer used columns if present
	df = df.drop(columns=["Table2asn", "GFF_File"], errors="ignore")
	# Remove duplicates and keep latest update
	df = df.drop_duplicates(subset = ["Submission_Name", "Organism", "Database", "Submission_Type", "Config_File"], keep = "last", ignore_index = True)
	# Force uppercase on columns if they are required to be uppercase
	df[["Organism", "Database", "Submission_Type", "Submission_Status", "Submission_ID"]] = df[["Organism", "Database", "Submission_Type", "Submission_Status", "Submission_ID"]].astype(str).apply(lambda col: col.str.upper())
	try:
		upload_schema.validate(df, lazy = True)
	except pandera.errors.SchemaErrors as schema_error:
		print("Error: Upload log columns are incorrect. Cannot process submissions.", file=sys.stderr)
		tools.pretty_print_pandera_errors(file=submission_log_file, error_msgs=[schema_error])
		sys.exit(1)
	return df

# Validate files listed inside submission_log.csv
def validate_fields_exist(df: pd.DataFrame):
	submission_dir = df["Submission_Directory"].iloc[0]
	submission_status_file = os.path.join(os.path.dirname(submission_dir), "submission_status_report.csv")
	config_file = df["Config_File"].iloc[0]
	file_handler.validate_directory(name = "directory", path = submission_dir)
	file_handler.validate_file(file_type = "submission_status_report", file_path = submission_status_file)
	file_handler.validate_file(file_type = "config file", file_path = config_file)

# Process submission status of existing biosample/sra database submission
def process_biosample_sra(submission_name: str, database: str, organism: str, submission_log_dir: str, submission_dir: str, curr_status: str, config_dict: dict[str, Any], submission_type: str) -> tuple[bool, str]:
	if curr_status == "PROCESSED":
		return True, curr_status
	report_file = ncbi_handler.get_ncbi_report(database=database, submission_name=submission_name, submission_dir=submission_dir, config_dict=config_dict, submission_type=submission_type)
	if report_file:
		new_submission_status, submission_id = biosample_sra_handler.process_biosample_sra_report(report_file=report_file, database=database, submission_dir=submission_dir)
	else:
		new_submission_status = curr_status
		submission_id = "PENDING"
	update_submission_log(database=database, organism=organism, submission_name=submission_name, submission_log_dir=submission_log_dir, submission_dir=submission_dir, submission_status=new_submission_status, submission_id=submission_id, submission_type=submission_type)
	if new_submission_status == "PROCESSED":
		return True, new_submission_status
	return False, new_submission_status

# Upload log submit GenBank database submission after previous required database submission
def upload_log_submit_genbank(genbank_type: str, submission_name: str, organism: str, submission_log_dir: str, submission_dir: str, config_dict: dict[str, Any], submission_type:str) -> str:
	if genbank_type == "GENBANK-TBL2ASN":
		submission_id = genbank_handler.create_table2asn(submission_name=submission_name, submission_dir=submission_dir)
		if submission_id == "VALIDATED":
			submission_status = ncbi_handler.email_table2asn(submission_name=submission_name, submission_dir=submission_dir, config_dict=config_dict, submission_type=submission_type)
		else:
			submission_status = "PENDING"
	elif genbank_type == "GENBANK-FTP":
		genbank_handler.create_zip(submission_name=submission_name, submission_dir=submission_dir)
		ncbi_handler.submit_ncbi(database="GENBANK", submission_name=submission_name, submission_dir=submission_dir, config_dict=config_dict, submission_type=submission_type)
		submission_id, submission_status = "PENDING", "SUBMITTED"
	else:
		print(f"Error: {genbank_type} is not a valid GenBank submission option.", file=sys.stderr)
		sys.exit(1)
	update_submission_log(database=genbank_type, organism=organism, submission_name=submission_name, submission_log_dir=submission_log_dir, submission_dir=submission_dir, submission_status=submission_status, submission_id=submission_id, submission_type=submission_type)
	return submission_status

# Process submission status of existing genbank database submission
def process_genbank(genbank_type: str, submission_name: str, submission_log_dir: str, submission_dir: str, curr_status: str, organism: str, config_dict: dict[str, Any], submission_type: str, linking_databases: dict[str, bool]) -> tuple[bool, str]:
	if curr_status in ["PROCESSED", "EMAILED"]:
		return True, curr_status
	elif curr_status == "WAITING" and submission_ready(submission_requirements=linking_databases, config_dict=config_dict, database="GENBANK"):
		if config_dict["Link_Sample_Between_NCBI_Databases"]:
			genbank_handler.update_genbank_files(linking_databases=linking_databases, organism=organism, submission_dir=submission_dir, config_dict=config_dict)
		new_submission_status = upload_log_submit_genbank(genbank_type=genbank_type, submission_name=submission_name, organism=organism, submission_log_dir=submission_log_dir, submission_dir=submission_dir, config_dict=config_dict, submission_type=submission_type)
		return False, new_submission_status
	elif curr_status == "WAITING":
		return False, curr_status
	elif curr_status == "PENDING" and genbank_type == "GENBANK-TBL2ASN":
		submission_id = genbank_handler.create_table2asn(submission_name=submission_name, submission_dir=submission_dir)
		if submission_id == "VALIDATED":
			new_submission_status = ncbi_handler.email_table2asn(submission_name=submission_name, submission_dir=submission_dir, config_dict=config_dict, submission_type=submission_type)
		else:
			new_submission_status = "PENDING"
		update_submission_log(database=genbank_type, organism=organism, submission_name=submission_name, submission_log_dir=submission_log_dir, submission_dir=submission_dir, submission_status=new_submission_status, submission_id=submission_id, submission_type=submission_type)
		return False, new_submission_status
	else:
		report_file = ncbi_handler.get_ncbi_report("GENBANK", submission_name, submission_dir, config_dict, submission_type)
		if report_file:
			new_submission_status, submission_id = genbank_handler.process_genbank_report(report_file=report_file, submission_dir=submission_dir)
		else:
			new_submission_status = curr_status
			submission_id = "PENDING"
		update_submission_log(database=genbank_type, organism=organism, submission_name=submission_name, submission_log_dir=submission_log_dir, submission_dir=submission_dir, submission_status=new_submission_status, submission_id=submission_id, submission_type=submission_type)
		if new_submission_status == "PROCESSED":
			return True, new_submission_status
		return False, new_submission_status

# Determine if upload log can submission to database
def submission_ready(submission_requirements: dict[str,bool], config_dict: dict[str, Any], database: str) -> bool:
	if submission_requirements["BIOSAMPLE"] and submission_requirements["SRA"]:
		return True
	else:
		return False

# Create SeqSender dict of the status for each DB under one submission name
def create_submission_requirements_dict(group_df: pd.DataFrame) -> dict[str, bool]:
	submission_requirements = dict()
	database_list = group_df["Database"].tolist()
	for database in ["BIOSAMPLE", "SRA"]:
		if database in database_list:
			if group_df.loc[group_df["Database"] == database, "Submission_Status"].iloc[0] == "PROCESSED":
				submission_requirements[database] = True
			else:
				submission_requirements[database] = False
		else:
			submission_requirements[database] = True
	if "GENBANK-FTP" in database_list:
		if group_df.loc[group_df["Database"] == "GENBANK-FTP", "Submission_Status"].iloc[0] == "PROCESSED":
			submission_requirements["GENBANK"] = True
		else:
			submission_requirements["GENBANK"] = False
	elif "GENBANK-TBL2ASN" in database_list:
		if group_df.loc[group_df["Database"] == "GENBANK-TBL2ASN", "Submission_Status"].iloc[0] == "PROCESSED":
			submission_requirements["GENBANK"] = True
		else:
			submission_requirements["GENBANK"] = False
	else:
		submission_requirements["GENBANK"] = True
	return submission_requirements

# Update all databases listed under one submission_name
def update_grouped_submission(group_df: pd.DataFrame, submission_log_dir: str, decrypt_key: str):
	validate_fields_exist(df=group_df)
	submission_requirements = create_submission_requirements_dict(group_df=group_df)
	# Reset index
	group_df = group_df.reset_index()
	submission_name = group_df.at[0, "Submission_Name"]
	submission_type = group_df.at[0, "Submission_Type"]
	submission_organism = group_df.at[0, "Organism"]
	submission_dir = group_df.at[0, "Submission_Directory"]
	databases = group_df["Database"].tolist()
	config_dict = tools.get_config(config_file=group_df.at[0, "Config_File"], databases=databases, decrypt_key=decrypt_key)
	if "BIOSAMPLE" in databases:
		biosample_status = group_df.loc[group_df["Database"] == "BIOSAMPLE", "Submission_Status"].iloc[0]
		submission_dir = group_df.loc[group_df["Database"] == "BIOSAMPLE", "Submission_Directory"].iloc[0]
		submission_requirements["BIOSAMPLE"], biosample_status = process_biosample_sra(submission_name=submission_name, organism=submission_organism, database="BIOSAMPLE", curr_status=biosample_status, submission_log_dir=submission_log_dir, submission_dir=submission_dir, config_dict=config_dict["NCBI"], submission_type=submission_type)
		print(f"\tBioSample: {biosample_status}")
	if "SRA" in databases:
		sra_status = group_df.loc[group_df["Database"] == "SRA", "Submission_Status"].iloc[0]
		submission_dir = group_df.loc[group_df["Database"] == "SRA", "Submission_Directory"].iloc[0]
		submission_requirements["SRA"], sra_status = process_biosample_sra(submission_name=submission_name, organism=submission_organism, database="SRA", curr_status=sra_status, submission_log_dir=submission_log_dir, submission_dir=submission_dir, config_dict=config_dict["NCBI"], submission_type=submission_type)
		print(f"\tSRA: {sra_status}")
	# Same requirements for GENBANK-FTP and GENBANK-TBL2ASN
	if any("GENBANK" in database for database in databases):
		if "GENBANK-FTP" in databases:
			genbank_type = "GENBANK-FTP"
		elif "GENBANK-TBL2ASN" in databases:
			genbank_type = "GENBANK-TBL2ASN"
		else:
			print(f"Error: Incorrect database option for GenBank in 'submission_log.csv' databases '{databases}' for '{submission_name}'.")
			sys.exit(1)
		genbank_status = group_df.loc[group_df["Database"] == genbank_type, "Submission_Status"].iloc[0]
		submission_dir = group_df.loc[group_df["Database"] == genbank_type, "Submission_Directory"].iloc[0]
		submission_requirements["GENBANK"], genbank_status = process_genbank(genbank_type=genbank_type, submission_name=submission_name, submission_log_dir=submission_log_dir, submission_dir=submission_dir, curr_status=genbank_status, organism=submission_organism, config_dict=config_dict["NCBI"], submission_type=submission_type, linking_databases=submission_requirements)
		print(f"\tGenBank: {genbank_status}")

# Update submission log, if given submission_name only update that specific submission
def update_submission_status(submission_dir: str, submission_name: Optional[str], decrypt_key: set[str]) -> None:
	df = load_submission_log(submission_dir)
	grouped_submissions = df.groupby(["Submission_Name", "Organism", "Submission_Type", "Config_File"])
	print("Checking Submissions:")
	for name, group in grouped_submissions:
		if not group["Submission_Status"].isin(["PROCESSED", "EMAILED"]).all() and submission_name is None or name[0] == submission_name:
			print(f"Submission: {name[0]}")
			for key in decrypt_key:
				try:
					update_grouped_submission(group_df=group, submission_log_dir=submission_dir, decrypt_key=key)
				except InvalidToken:
					print(f"Unable to access config file with key for {name}. Continuing to next key/submission.", file=sys.stderr)
				except Exception as e:
					print(f"Error: Unable to process {name} because:\n{e}", file=sys.stderr)
					break # If error is not related to decrypt_key then skip to next submission
	print("\nUpdating submissions complete.")
