#!/usr/bin/env python3

###########################	Description	##################################
# Functions for handling update log
################################################################################

import file_handler
import pandas as pd
from config.seqsender_upload_log_schema import schema as upload_schema
from typing import List, Optional
from pandera import pandera, DataFrameSchema, Column, Check, Index, MultiIndex
import os
import sys

# create the submission df for biosample
def create_submission_status_csv(database: List[str], sequence_names: pd.DataFrame, submission_status_file: str) -> None:
	status_submission_df = sequence_names
	if "BIOSAMPLE" in database:
		status_submission_df["biosample_status"] = ""
		status_submission_df["biosample_accession"] = ""
		status_submission_df["biosample_message"] = ""
	if "SRA" in database:
		status_submission_df["sra_status"] = ""
		status_submission_df["sra_accession"] = ""
		status_submission_df["sra_message"] = ""
	if "GENBANK" in database:
		status_submission_df["genbank_status"] = ""
		status_submission_df["genbank_accession"] = ""
		status_submission_df["genbank_message"] = ""
	if "GISAID" in database:
		status_submission_df["gisaid_accession_epi_isl_id"] = ""
		status_submission_df["gisaid_accession_epi_id"] = ""
		status_submission_df["gisaid_message"] = ""
	# Save df
	status_submission_df.to_csv(submission_status_file, header = True, index = False)

# Create submission log csv
def create_submission_log(database: str, submission_position: int, organism: str, submission_name: str, submission_dir: str, config_file: str, submission_status: str, submission_id: str, submission_type: str) -> None:
	# If file doesn't exist create it
	if os.path.isfile(os.path.join(submission_dir, "submission_log.csv")) == True:
		df = pd.read_csv(os.path.join(submission_dir, "submission_log.csv"), header = 0, dtype = str, engine = "python", encoding="utf-8", index_col=False)
	else:
		df = pd.DataFrame(columns = ["Submission_Name", "Organism", "Database", "Submission_Position", "Submission_Type", "Submission_Date", "Submission_Status", "Submission_Directory", "Config_File", "Table2asn", "GFF_File", "Update_Date"])
	# Fill in the log field if it exists, otherwise create new
	df_partial = df.loc[(df["Organism"] == organism) & (df["Database"] == database) & (df["Submission_Directory"] == submission_dir) & (df["Submission_Name"] == submission_name) & (df["Submission_Type"] == submission_type)]
	# Update existing field
	if df_partial.shape[0] > 0:
		df.loc[df_partial.index.values, "Submission_Position"] = submission_position
		df.loc[df_partial.index.values, "Submission_Status"] = submission_id + ";" + submission_status
		df.loc[df_partial.index.values, 'Update_Date'] = datetime.now().strftime("%Y-%m-%d")
	else:
		# Create new field
		status = submission_id + ";" + submission_status
		new_entry = {'Submission_Name': submission_name,
					 'Organism': organism,
					 'Database': database,
					 'Submission_Type': submission_type,
					 'Submission_Date': datetime.now().strftime("%Y-%m-%d"),
					 'Submission_Status': status,
					 'Submission_Directory': submission_dir,
					 'Config_File': config_file,
					 'Update_Date': datetime.now().strftime("%Y-%m-%d")
					}
		df.loc[len(df)] = new_entry # type: ignore[call-overload]
	df.to_csv(os.path.join(submission_dir, "submission_log.csv"), header = True, index = False)

# Validate submission_dir and log exists
def load_upload_log(submission_dir: str) -> pd.DataFrame:
	file_handler.validate_directory(name = "directory", path = submission_dir)
	submission_log_file = os.path.join(submission_dir, "submission_log.csv")
	file_handler.validate_file("submission_log", submission_log_file)
	df = file_handler.load_csv(submission_log_file)
	try:
		upload_schema.validate(df, lazy = True)
	except pandera.errors.SchemaErrors as schema_error:
		print("Error: Upload log columns are incorrect. Cannot process submissions.", file=sys.stderr)
		print(schema_error, file=sys.stderr)
		sys.exit(1)
	return df

def update_individual_submission(submission_dir, submission_name):

# Update submission log
def update_submission_status(submission_dir: str, submission_name: Optional[str]) -> None:
	df = load_upload_log(submission_dir)
	grouped_submissions = df.groupby(["Submission_Name", "Organism", "Submission_Type", "Submission_Directory", "Config_File"])
	for name, group in grouped_submissions:
		print("Checking Submission:")			
		print(f"Group_name: {name}")
		print(group)
	sys.exit(1)
	df_partial = df.loc[(df["Organism"] == organism) & (df["Submission_Name"] == submission_name) & (df["Submission_Directory"] == submission_dir) & (df["Submission_Type"] == submission_type)]
	if df_partial.shape[0] == 0:
		print("Error: Submission name: " + submission_name + " for "+organism+" "+submission_type+"-data is not found in the submission log file.", file=sys.stderr)
		print("Error: Either a submission has not been made or an entry has been moved.", file=sys.stderr)
		sys.exit(1)
	# Order get a list of submitting databases
	df_partial = df_partial.sort_values(by=["Submission_Position"])
	database = df_partial["Database"].tolist()
	# Output message
	print("\n"+"Checking submission status for:"+"\n", file=sys.stdout)
	print("Submission name: " + submission_name, file=sys.stdout)
	print("Submission organism: " + organism, file=sys.stdout)
	print("Submission type: " + submission_type, file=sys.stdout)
	# Check the status of each database in its order of submission
	for database_name in database:
		print("\n" + "Submission database: " + database_name, file=sys.stdout)
		df = pd.read_csv(submission_log_file, header = 0, dtype = str, engine = "python", encoding="utf-8", index_col=False).sort_values('Submission_Position', ascending=True)
		df_processing = df[(df["Organism"] == organism) & (df["Database"] == database_name) & (df["Submission_Directory"] == submission_dir) & (df["Submission_Name"] == submission_name) & (df["Submission_Type"] == submission_type)]
		df_processing = df_processing.reset_index(drop=True)
		submission_dir = df_processing["Submission_Directory"][0]
		submission_position = df_processing["Submission_Position"][0]
		submission_id, submission_status = df_processing["Submission_Status"][0].strip().split(";")
		config_file = df_processing["Config_File"][0]
		table2asn = df_processing["Table2asn"][0]
		gff_file = df_processing["GFF_File"][0]
		# Check if submission files exist in parent directory
		submission_files_dir = os.path.join(submission_dir, submission_name, "submission_files", database_name)
		if os.path.exists(submission_files_dir) == False:
			print("Error: Submission files for "+submission_name+" does not exist at "+submission_files_dir, file=sys.stderr)
			sys.exit(1)
		# Check if submission report status csv exists
		submission_status_file = os.path.join(submission_dir, submission_name, "submission_report_status.csv")
		if os.path.isfile(submission_status_file) == False:
			print("Error: Submission status report for "+submission_name+" does not exist at "+submission_status_file, file=sys.stderr)
			sys.exit(1)
		# Check if config file exists
		if os.path.isfile(config_file) == False:
			print("Error: Config file for "+submission_name+" does not exist at "+config_file, file=sys.stderr)
			sys.exit(1)
		else:
			config_dict = get_config(config_file=config_file, database=database)
		# IF GISAID in a list of submitting databases, check if CLI is downloaded and store in the correct directory
		gisaid_cli = None
		if "GISAID" in database_name:
			gisaid_cli = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
		# Check if the gisaid_cli exists
		if (gisaid_cli is not None) and os.path.isfile(gisaid_cli) == False:
			print("There is no GISAID CLI package for " + organism + " located at "+ gisaid_cli, file=sys.stderr)
			print("Please download the CLI package from GISAID platform", file=sys.stderr)
			print("Then place a copy of the CLI binary at "+ gisaid_cli, file=sys.stderr)
			sys.exit(1)
		# Check the status of the submission
		if "processed-ok" in submission_status:
			print("Submission status: " + submission_status, file=sys.stdout)
		else:
			# Pull download submission report and update its status
			if database_name in ["BIOSAMPLE", "SRA", "GENBANK"]:
				# If report exists, processing the report and output status of the submission
				if database_name in ["BIOSAMPLE", "SRA"]:
					report_file = report.get_ncbi_process_report(database=database_name, submission_name=submission_name, submission_files_dir=submission_files_dir, config_dict=config_dict['NCBI'], submission_type=submission_type)
					if report_file is not None and os.path.isfile(report_file):
						submission_status, submission_id = report.process_biosample_sra_report(report_file=report_file, submission_status_file=submission_status_file)
				elif database_name == "GENBANK":
					# Update submission
					if "---" in submission_status:
						# Check if biosample, sra, and gisaid are in a list of submitting databases
						if int(submission_position) == 1:
							other_submitting_db = [x for x in database if x in ["BIOSAMPLE", "SRA"]]
						elif int(submission_position) == 2:
							other_submitting_db = [x for x in database if x in ["BIOSAMPLE", "SRA", "GISAID"]]
						# Update biosample, sra, gisaid accession on genbank submission
						if len(other_submitting_db) > 0:
							all_status = []
							for db in other_submitting_db:
								db_df = df.loc[df["Database"] == db]
								db_df = db_df.reset_index(drop=True)
								db_status = db_df["Submission_Status"][0]
								# If the status of biosample or sra is processed-ok, then go ahead and submit to Genbank
								if "processed-ok" in db_status:
									all_status += [1]
									report.update_genbank_files(database=database, organism=organism, submission_files_dir=submission_files_dir, submission_status_file=submission_status_file)
								else:
									all_status += [0]
							# Submit via Table2asn
							if all(all_status):
								if table2asn == True:
									if os.path.isfile(str(gff_file)) == False:
										gff_file = None
									submission_id, submission_status = create.create_genbank_table2asn(submission_name=submission_name, submission_files_dir=submission_files_dir, gff_file=gff_file)
									if submission_status == "validated":
										submission_status = submit.sendmail(database=database_name, submission_name=submission_name, submission_dir=submission_dir, config_dict=config_dict['NCBI'], test=test)
								else:
									# Submit via FTP
									create.create_genbank_zip(submission_name=submission_name, submission_files_dir=submission_files_dir)
									submit.submit_ncbi(database=database_name, submission_name=submission_name, submission_dir=submission_dir, config_dict=config_dict["NCBI"], submission_type=submission_type)
									submission_status = "submitted"
									submission_id = "pending"
					else:
						# If report exists, processing the report and output status of the submission
						report_file = report.get_ncbi_process_report(database=database_name, submission_name=submission_name, submission_files_dir=submission_files_dir, config_dict=config_dict['NCBI'], submission_type=submission_type)
						if report_file is not None and os.path.isfile(report_file):
							submission_status, submission_id = report.process_genbank_report(report_file=report_file, submission_status_file=submission_status_file, submission_files_dir=submission_files_dir)
			elif database_name == "GISAID":
				if ("---" in submission_status) and (int(submission_position) == 1):
					assert isinstance(gisaid_cli, str)
					submission_status = submit.submit_gisaid(organism=organism, database=database_name, submission_dir=submission_dir, submission_name=submission_name, config_dict=config_dict["GISAID"], gisaid_cli=gisaid_cli, submission_status_file=submission_status_file, submission_type=submission_type)
					submission_id = ""
				elif ("---" in submission_status) and ("GENBANK" in database) and (int(submission_position) == 2):
					db_status = df[(df["Database"] == "GENBANK"), "Submission_Status"][0]
					if "processed-ok" in db_status:
						report.update_gisaid_files(organism=organism, submission_files_dir=submission_files_dir, submission_status_file=submission_status_file)
						assert isinstance(gisaid_cli, str)
						submission_status = submit.submit_gisaid(organism=organism, database=database_name, submission_dir=submission_dir, submission_name=submission_name, config_dict=config_dict["GISAID"], gisaid_cli=gisaid_cli, submission_status_file=submission_status_file, submission_type=submission_type)
						submission_id = ""
			# Update status in the submission log
			create.create_submission_log(database=database_name, submission_position=submission_position, organism=organism, submission_name=submission_name, submission_dir=submission_dir, config_file=config_file, submission_status=submission_status, submission_id=submission_id, submission_type=submission_type)
			# Print out the submission status
			print("Submission status: " + submission_status, file=sys.stdout)
