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
from typing import List, Dict, Set, Optional, Tuple, Any

# Local imports
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import process
import setup
import file_handler
import argument_handler
import ncbi_handler
import genbank_handler
import biosample_sra_handler
import gisaid_handler
import upload_log


# Get program directory
PROG_DIR: str = os.path.dirname(os.path.abspath(__file__))

# Define version of seqsender
VERSION: str = "1.1.0 (Beta)"

# Define current time
STARTTIME = datetime.now()

# Define organsim choices
ORGANISM_CHOICES: List[str] = ["FLU", "COV", "POX", "ARBO", "OTHER"]

# Define database choices
DATABASE_CHOICES: List[str] = ["BIOSAMPLE", "SRA", "GENBANK", "GISAID"]

# Get execution time
def get_execution_time() -> None:
	print(f"\nTotal runtime (HRS:MIN:SECS): {str(datetime.now() - STARTTIME)}")

# Setup needed requirements for running
def prep(database: List[str], organism: str, submission_dir: str, submission_name: str, config_file: str, metadata_file: str, fasta_file: Optional[str], gff_file: Optional[str]) -> Tuple[Dict[str, Any], pd.DataFrame]:
	# Create the appropriate files
	file_dict = {
		"config_file" : os.path.join(submission_dir, submission_name, config_file),
		"metadata_file" : os.path.join(submission_dir, submission_name, metadata_file),
		"fasta_file" : os.path.join(submission_dir, submission_name, str(fasta_file)) if fasta_file is not None else None,
		"gff_file" : os.path.join(submission_dir, submission_name, str(gff_file)) if gff_file is not None else None
	}
	file_handler.validate_files(file_dict, organism)
	submission_status_file = os.path.join(submission_dir, submission_name, "submission_report_status.csv")
	# Check if submission directory exists
	file_handler.validate_directory(name = "submission directory", path = submission_dir)
	# Determine whether this is a test or production submission
	# if test is True:
	# 	submission_type = "Test"
	# else:
	# 	submission_type = "Production"
	# load config file
	config_dict = process.get_config(config_file=config_file, database=database)
	# load metadata file
	metadata = process.get_metadata(database=database, organism=organism, metadata_file=metadata_file, config_dict=config_dict)
	# Load fasta file into metadata
	if fasta_file:
		metadata = process.process_fasta_samples(metadata = metadata, fasta_file = fasta_file)
	# Prepping submission files for each given database
	for database_name in database:
		print(f"Creating submission files for {database_name}", file=sys.stdout)
		submission_files_dir = os.path.join(submission_dir, submission_name, "submission_files", database)
		file_handler.create_directory(submission_files_dir)
		if database_name in ["BIOSAMPLE", "SRA"]:
			biosample_sra_handler.create_biosample_sra_submission(organism=organism, database=database_name, submission_name=submission_name, submission_dir=submission_dir, config_dict=config_dict["NCBI"], metadata=metadata)
		elif database_name == "GENBANK":
			genbank_handler.create_genbank_submission(organism=organism, submission_name=submission_name, submission_dir=submission_dir, config_dict=config_dict["NCBI"], metadata=metadata, gff_file=gff_file)
		elif database_name == "GISAID":
			gisaid_handler.create_gisaid_submission(organism=organism, database=database_name, submission_name=submission_name, submission_dir=submission_dir, config_dict=config_dict["GISAID"], metadata=metadata)
		else:
			print(f"Error: Database {database_name} is not a valid database selection.", file=sys.stderr)
			sys.exit(1)
		print(f"Files are stored at: {submission_files_dir}", file=sys.stdout)
	return (config_dict, metadata)

# Setup needed requirements for running
def start(database: List[str], organism: str, submission_dir: str, submission_name: str, config_file: str, metadata_file: str, fasta_file: Optional[str], gff_file: Optional[str], test: bool = False) -> None:
	# IF database is GISAID, check if CLI is in the correct directory
	if "GISAID" in database:
		file_handler.validate_gisaid_installer(submission_dir, organism)
	config_dict, metadata = prep(database, organism, submission_dir, submission_name, config_file, metadata_file, fasta_file, gff_file)
	create.create_submission_status_csv(database=database, sequence_names=sequence_names, submission_status_file=submission_status_file)
	for database_name in database:
		if database_name in ["BIOSAMPLE", "SRA", "GENBANK"]:
			submit.submit_ncbi(submission_name=submission_name, submission_dir=submission_dir, database=database_name, config_dict=config_dict["NCBI"], submission_type=submission_type)
		elif "GISAID" in database_name:
			submit.submit_gisaid(organism=organism, database=database_name, submission_dir=submission_dir, submission_name=submission_name, config_dict=config_dict["GISAID"], gisaid_cli=gisaid_cli,  submission_status_file=submission_status_file, submission_type=submission_type)
		else:
			print(f"Error: Database selection {database_name} is not valid.", file=sys.stderr)
			sys.exit(1)
		create.create_submission_log(database=database_name, submission_position=submission_position, organism=organism, submission_name=submission_name, submission_dir=submission_dir, config_file=config_file, submission_status=submission_status, submission_id=submission_id, submission_type=submission_type)

def main():
	"""The main routine"""
	parser, submit_prep_subparser  = argument_handler.args_parser(ORGANISM_CHOICES)
	args = parser.parse_args()

	# Parse the command argument
	command = args.command

	# Determine databases
	if command in ["prep", "submit", "test_data"]:
		database = []
		if args.biosample != None:
			database += [args.biosample]
		if args.sra:
			database += [args.sra]
		if args.genbank:
			database += [args.genbank]
		if args.gisaid:
			database += [args.gisaid]

	# Clean submission dir
	submission_dir = os.path.abspath(args.submission_dir)

	# Execute the command
	if command == "prep":
		prep(organism=args.organism, database=database, submission_name=args.submission_name, submission_dir=submission_dir, config_file=args.config_file, metadata_file=args.metadata_file, fasta_file=args.fasta_file, gff_file=args.gff_file)
	elif command == "submit":
		# If database is not given, display help
		if len(database) == 0:
			print("\n"+"ERROR: Missing a database selection. See USAGE below."+"\n", file=sys.stdout)
			submit_prep_subparser.print_help()
			sys.exit(0)
		start(organism=args.organism, database=database, submission_name=args.submission_name, submission_dir=args.submission_dir, config_file=args.config_file, metadata_file=args.metadata_file, fasta_file=args.fasta_file, gff_file=args.gff_file, test=args.test)
	elif command == "check_submission_status":
		upload_log.update_submission_status(submission_dir=args.submission_dir, submission_name=args.submission_name)
	elif command == "test_data":
		# If database is not given, display help
		if len(database) == 0:
			print("\n"+"ERROR: Missing a database selection. See USAGE below."+"\n", file=sys.stdout)
			submit_prep_subparser.print_help()
			sys.exit(0)
		setup.create_test_data(organism=args.organism, database=database, submission_dir=args.submission_dir)
	elif command == "version":
		print("\n"+"Version: " + VERSION, file=sys.stdout)
		sys.exit(0)
	elif command == "update_biosample":
		print("Updating BioSample requirements.", file=sys.stdout)
		setup.download_biosample_xml_list()
		sys.exit(0)
	else:
		# If no command display help
		parser.print_help()
		sys.exit(0)

	# Print out the execution time
	get_execution_time()

if __name__ == "__main__":
	main()
