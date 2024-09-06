#!/usr/bin/env python3

# Python Libraries
import pathlib
import os
import sys
from datetime import datetime
import subprocess
import argparse
import pandas as pd
from distutils.util import strtobool
from typing import List, Dict, Set, Optional, Tuple, Any, TypedDict

# Local imports
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import setup
import file_handler
import argument_handler
import ncbi_handler
import genbank_handler
import biosample_sra_handler
import gisaid_handler
import upload_log
import tools

from settings import VERSION

import sys

# Define current time
STARTTIME = datetime.now()

# Get execution time
def get_execution_time() -> None:
	print(f"\nTotal runtime (HRS:MIN:SECS): {str(datetime.now() - STARTTIME)}")

# Setup needed requirements for running
def prep(database: List[str], organism: str, submission_dir: str, submission_name: str, config_file: str, metadata_file: str, fasta_file: Optional[str], gff_file: Optional[str], table2asn: bool, skip_validation: bool = False) -> Tuple[str, Dict[str, Any], pd.DataFrame]:
	# Create the appropriate files
	File_Dict = TypedDict("File_Dict", {"config_file": str, "metadata_file": str, "fasta_file": Optional[str], "gff_file": Optional[str]})
	file_dict: File_Dict = {
		"config_file" : config_file,
		"metadata_file" : metadata_file,
		"fasta_file" : fasta_file,
		"gff_file" : gff_file
	}
	# Check if submission directory exists
	file_handler.validate_directory(name = "submission directory", path = submission_dir)
	# Validate files
	for file_type, file_path in file_dict.items():
		if file_type in ["fasta_file", "gff_file"] and file_path is None:
			# If not provided
			continue
		elif isinstance(file_path, str) and (os.path.sep in file_path or (os.path.altsep and os.path.altsep in file_path)):
			# If full file path given
			file_handler.validate_file(file_type=file_type, file_path=file_path)
		else:
			# If only file name given and stored at submission_dir
			assert isinstance(file_path, str)
			updated_path = os.path.join(submission_dir, submission_name, file_path)
			file_handler.validate_file(file_type=file_type, file_path=updated_path)
			file_dict[file_type] = updated_path # type: ignore
	# load config file
	config_dict = tools.get_config(config_file=file_dict["config_file"], database=database)
	# warn if user is submitting biosample & sra together with Link_Sample_Between_NCBI_Databases off
	if not config_dict["NCBI"]["Link_Sample_Between_NCBI_Databases"] and 'SRA' in database and 'BIOSAMPLE' in database:
		print("Warning: You are submitting to BioSample and SRA together, but Link_Sample_Between_NCBI_Databases is off. If you did not intend to do this, SRA submission is likely to fail.\nIf this was unintentional, cancel the submission and set 'Link_Sample_Between_NCBI_Databases: on' in your config yaml.")
	# load metadata file
	metadata = tools.get_metadata(database=database, organism=organism, metadata_file=file_dict["metadata_file"], config_dict=config_dict, skip_validation=skip_validation)
	# Load fasta file into metadata
	if file_dict["fasta_file"]:
		metadata = tools.process_fasta_samples(metadata=metadata, fasta_file=file_dict["fasta_file"])
	# Prepping submission files for each given database
	for database_name in database:
		print(f"Creating submission files for {database_name}", file=sys.stdout)
		database_dir = os.path.join(submission_dir, submission_name, "submission_files", database_name)
		file_handler.create_directory(database_dir)
		if database_name in ["BIOSAMPLE", "SRA"]:
			biosample_sra_handler.create_biosample_sra_submission(organism=organism, database=database_name, submission_name=submission_name, submission_dir=database_dir, config_dict=config_dict["NCBI"], metadata=metadata)
		elif database_name == "GENBANK":
			genbank_handler.create_genbank_submission(organism=organism, submission_name=submission_name, submission_dir=database_dir, config_dict=config_dict["NCBI"], metadata=metadata, gff_file=file_dict["gff_file"], table2asn=table2asn)
		elif database_name == "GISAID":
			gisaid_handler.create_gisaid_files(organism=organism, database=database_name, submission_name=submission_name, submission_dir=database_dir, config_dict=config_dict["GISAID"], metadata=metadata)
		else:
			print(f"Error: Database {database_name} is not a valid database selection.", file=sys.stderr)
			sys.exit(1)
		print(f"Files are stored at: {database_dir}", file=sys.stdout)
	return (file_dict["config_file"], config_dict, metadata)

# Setup needed requirements for running
def submit(database: List[str], organism: str, submission_dir: str, submission_name: str, config_file: str, metadata_file: str, fasta_file: Optional[str], gff_file: Optional[str], table2asn: bool = False, test: bool = False, skip_validation: bool = False) -> None:
	# IF database is GISAID, check if CLI is in the correct directory
	if "GISAID" in database:
		file_handler.validate_gisaid_installer(submission_dir, organism)
	config_file_path, config_dict, metadata = prep(database=database, organism=organism, submission_dir=submission_dir, submission_name=submission_name, config_file=config_file, metadata_file=metadata_file, fasta_file=fasta_file, gff_file=gff_file, table2asn=table2asn, skip_validation=skip_validation)
	print("", file=sys.stdout)
	upload_log.create_submission_status_csv(database=database, metadata=metadata, submission_dir=os.path.join(submission_dir, submission_name, "submission_files"))
	submission_type = tools.get_submission_type(test=test)
	for database_name in database:
		submission_status, submission_id = "WAITING", "PENDING"
		database_dir = os.path.join(submission_dir, submission_name, "submission_files", database_name)
		file_handler.validate_directory(name="submission directory", path=database_dir)
		if database_name in ["BIOSAMPLE", "SRA"]:
			ncbi_handler.submit_ncbi(submission_name=submission_name, submission_dir=database_dir, database=database_name, config_dict=config_dict["NCBI"], submission_type=submission_type)
			submission_status = "SUBMITTED"
		elif "GENBANK" in database_name:
			sub_pos = tools.get_submission_position(config_dict=config_dict, database="GENBANK")
			link_ncbi = config_dict["NCBI"]["Link_Sample_Between_NCBI_Databases"]
			if "BIOSAMPLE" in database or "SRA" in database:
				ncbi_other_databases = True
			else:
				ncbi_other_databases = False
			if table2asn:
				database_name = "GENBANK-TBL2ASN"
			else:
				database_name = "GENBANK-FTP"
			if (sub_pos is None or sub_pos == 1) and (not ncbi_other_databases or not link_ncbi):
				if table2asn:
					submission_status = ncbi_handler.email_table2asn(submission_name=submission_name, submission_dir=database_dir, config_dict=config_dict["NCBI"], submission_type=submission_type)
				else:
					ncbi_handler.submit_ncbi(submission_name=submission_name, submission_dir=database_dir, database="GENBANK", config_dict=config_dict["NCBI"], submission_type=submission_type)
					submission_status = "SUBMITTED"
		elif "GISAID" in database_name:
			sub_pos = tools.get_submission_position(config_dict=config_dict, database="GISAID")
			if sub_pos is None or sub_pos == 1:
				submission_status = gisaid_handler.submit_gisaid(organism=organism, submission_dir=database_dir, submission_name=submission_name, config_dict=config_dict["GISAID"], submission_type=submission_type)
		else:
			print(f"Error: Database selection {database_name} is not valid.", file=sys.stderr)
			sys.exit(1)
		upload_log.create_submission_log(database=database_name, organism=organism, submission_name=submission_name, submission_dir=submission_dir, database_dir=database_dir, config_file=config_file_path, submission_status=submission_status, submission_id=submission_id, submission_type=submission_type)

def main():
	"""The main routine"""
	print("")
	parser = argument_handler.args_parser()
	args = parser.parse_args()

	# Parse the command argument
	command = args.command

	# Determine databases selected
	if command in ["prep", "submit", "test_data"]:
		database = []
		if args.biosample:
			database += [args.biosample]
		if args.sra:
			database += [args.sra]
		if args.genbank:
			database += [args.genbank]
		if args.gisaid:
			print("here")
			database += [args.gisaid]
		if len(database) == 0:
			print("ERROR: Missing a required database selection. See USAGE below.", file=sys.stderr)
			parser.print_help()
			sys.exit(0)

	if command in ["prep", "submit", "submission_status", "test_data"]:
		submission_dir = os.path.abspath(args.submission_dir)

	# Execute the command
	if command == "prep":
		prep(organism=args.organism, database=database, submission_name=args.submission_name, submission_dir=submission_dir, config_file=args.config_file, metadata_file=args.metadata_file, fasta_file=args.fasta_file, gff_file=args.gff_file, table2asn=args.table2asn, skip_validation=args.skip_validation)
	elif command == "submit":
		submit(organism=args.organism, database=database, submission_name=args.submission_name, submission_dir=submission_dir, config_file=args.config_file, metadata_file=args.metadata_file, fasta_file=args.fasta_file, gff_file=args.gff_file, table2asn=args.table2asn, test=args.test, skip_validation=args.skip_validation)
	elif command == "submission_status":
		upload_log.update_submission_status(submission_dir=submission_dir, submission_name=args.submission_name)
	elif command == "test_data":
		setup.create_test_data(organism=args.organism, database=database, submission_dir=submission_dir)
	elif command == "version":
		print("Version: {VERSION}", file=sys.stdout)
	elif command == "update_biosample":
		print("Updating BioSample requirements.", file=sys.stdout)
		setup.download_biosample_xml_list()
	else:
		# If no command display help
		parser.print_help()
		sys.exit(0)

	# Print out the execution time
	get_execution_time()

if __name__ == "__main__":
	main()
