#!/usr/bin/env python3

# Python Libraries
import pathlib
import os
import sys
from datetime import datetime
import subprocess
import argparse
from distutils.util import strtobool

# Local imports
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import create
import process
import report
import submit
import setup

# Get program directory
PROG_DIR = os.path.dirname(os.path.abspath(__file__))

# Define version of seqsender
VERSION = "1.1.0 (Beta)"

# Define current time
STARTTIME = datetime.now()

# Define organsim choices
ORGANISM_CHOICES = ["FLU", "COV"]

# Define database choices
DATABASE_CHOICES = ["BIOSAMPLE", "SRA", "GENBANK", "GISAID"]

# Get execution time
def get_execution_time():
	print(f"\nTotal runtime (HRS:MIN:SECS): {str(datetime.now() - STARTTIME)}")

# Setup needed requirements for running
def start(command, database, organism, submission_dir, submission_name, config_file, metadata_file, fasta_file=None, table2asn=False, gff_file=None, test=False):
	# Create the appropriate files
	submission_dir = os.path.abspath(submission_dir)
	config_file = os.path.join(submission_dir, submission_name, config_file)
	metadata_file = os.path.join(submission_dir, submission_name, metadata_file)
	fasta_file = os.path.join(submission_dir, submission_name, str(fasta_file)) if fasta_file is not None else None
	gff_file = os.path.join(submission_dir, submission_name, str(gff_file)) if gff_file is not None else None
	submission_status_file = os.path.join(submission_dir, submission_name, "submission_report_status.csv")
	# Check if submission directory exists
	if os.path.exists(submission_dir) == False:
		print("There is no submission directory at " + submission_dir, file=sys.stderr)
		sys.exit(1)
	# Check if config file exists
	if os.path.isfile(config_file) == False:
		print("There is no config file at " + config_file, file=sys.stderr)
		sys.exit(1)
	# Check if metadata file exists
	if os.path.isfile(metadata_file) == False:
		print("There is no metadata at " + metadata_file, file=sys.stderr)
		sys.exit(1)
	# If fasta file is provided, check if file exists
	if (fasta_file is not None) and (os.path.isfile(fasta_file) == False):
		print("There is no fasta file at " + fasta_file, file=sys.stderr)
		sys.exit(1)
	# If table2asn is true, if gff file is given, check if file exists
	if (table2asn == True) and (gff_file is not None) and (os.path.isfile(gff_file) == False):
		print("Error: gff file does not exist at: " + gff_file, file=sys.stderr)
		sys.exit(1)
	# IF database is GISAID, check if CLI is downloaded and store in the correct directory
	gisaid_cli = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI") if "GISAID" in database else None
	# Check if the gisaid_cli exists
	if (gisaid_cli is not None) and os.path.isfile(gisaid_cli) == False:
		print("There is no GISAID CLI package for " + organism + " located at "+ gisaid_cli, file=sys.stderr)
		print("Please download the GISAID " + organism + " CLI package from the GISAID platform", file=sys.stderr)
		print("Extract the zip file and place a copy of the CLI binary at "+ gisaid_cli, file=sys.stderr)
		sys.exit(1)	
	# Determine whether this is a test or production submission
	if test is True:
		submission_type = "Test"
	else:
		submission_type = "Production"
	# Check config file
	config_dict = process.get_config(config_file=config_file, database=database)
	# Check metadata file
	metadata = process.get_metadata(database=database, organism=organism, metadata_file=metadata_file)
	# Create identifier for each database to store submitting samples in submission status worksheet
	identifier_columns = dict()
	# Prepping submission files for each given database
	for database_name in database:
		if database_name in ["BIOSAMPLE", "SRA", "GENBANK"]:
			identifier_columns.update({"ncbi-spuid": "ncbi-sample_name"})
			create.create_ncbi_submission(organism=organism, database=database_name, submission_name=submission_name, submission_dir=submission_dir, config_dict=config_dict["NCBI"], metadata=metadata, fasta_file=fasta_file, table2asn=table2asn, gff_file=gff_file)
			if "GENBANK" in database_name:
				identifier_columns.update({"gb-seq_id": "ncbi-sequence_name"})
		elif "GISAID" in database_name:
			if "FLU" in organism:
				identifier_columns.update({"gs-Isolate_Name": "gs-sample_name"})
				identifier_columns.update({"gs-seq_id": "gs-sequence_name"})
			elif "COV" in organism:
				identifier_columns.update({"gs-virus_name": "gs-sample_name"})
			create.create_gisaid_submission(organism=organism, database=database_name, submission_name=submission_name, submission_dir=submission_dir, config_dict=config_dict["GISAID"], metadata=metadata, fasta_file=fasta_file)
		else:
			print("Error: Database " + database_name + " is not a valid database selection.", file=sys.stderr)
			sys.exit(1)
	# Extract the samples information from metadata
	sequence_names = metadata[identifier_columns.keys()].copy()
	sequence_names = sequence_names.rename(columns=identifier_columns)
	create.create_submission_status_csv(database=database, sequence_names=sequence_names, submission_status_file=submission_status_file)
	# Submit to databases
	if command == "submit":
		for database_name in database:
			# BioSample and SRA can be submitted together but to add accessions to GenBank they must be fully processed
			if database_name in ["BIOSAMPLE", "SRA"]:	
				if ("GISAID" in database) and (int(config_dict["GISAID"]["Submission_Position"]) == 1):
					submission_position = 2
				else:
					submission_position = 1
				submission_status = "submitted"
				submission_id = "pending"
				submit.submit_ncbi(submission_name=submission_name, submission_dir=submission_dir, database=database_name, config_dict=config_dict["NCBI"], submission_type=submission_type)
			elif "GENBANK" in database_name:
				if ("GISAID" in database) and (int(config_dict["GISAID"]["Submission_Position"]) == 1):
					submission_position = 2
				else:
					submission_position = 1
				# Don't submit to GENBANK yet if BIOSAMPLE OR SRA is given or if GISAID is given and the submission_position is set to primary
				if any(x in ["BIOSAMPLE", "SRA"] for x in database) or (submission_position == 2):
					submission_status = "---"
					submission_id = "---"
				elif table2asn == True:
					# GenBank Table2asn submission
					submit.sendmail(database=database_name, submission_name=submission_name, submission_dir=submission_dir, config_dict=config_dict["NCBI"], test=test)
					submission_status = "processed-ok"
					submission_id = "Table2asn"
				else:
					# GenBank FTP submission
					submit.submit_ncbi(database=database_name, submission_name=submission_name, submission_dir=submission_dir, config_dict=config_dict["NCBI"], submission_type=submission_type)
					submission_status = "submitted"
					submission_id = "pending"
			elif "GISAID" in database_name:
				# Don't submit to GISAID yet if any of NCBI databases is given and the submission position is set to secondary
				if  any(x in ["BIOSAMPLE", "SRA", "GENBANK"] for x in database) and int(config_dict["GISAID"]["Submission_Position"]) == 2:
					submission_position = 2
					submission_status = "---"
					submission_id = ""
				else:
					submission_position = 1
					submission_status = submit.submit_gisaid(organism=organism, database=database_name, submission_dir=submission_dir, submission_name=submission_name, config_dict=config_dict["GISAID"], gisaid_cli=gisaid_cli,  submission_status_file=submission_status_file, submission_type=submission_type)
					submission_id = ""
			create.create_submission_log(database=database_name, submission_position=submission_position, organism=organism, submission_name=submission_name, submission_dir=submission_dir, config_file=config_file, submission_status=submission_status, submission_id=submission_id, table2asn=table2asn, gff_file=gff_file, submission_type=submission_type)

# Define program arguments and commands
def args_parser():
	"""
	Argument parser setup and build.
	"""
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
									description='Automate the process of batch uploading consensus sequences and metadata to databases of your choices')
	database_parser = argparse.ArgumentParser(add_help=False)
	organism_parser = argparse.ArgumentParser(add_help=False)
	submission_name_parser = argparse.ArgumentParser(add_help=False)
	submission_dir_parser = argparse.ArgumentParser(add_help=False)
	config_file_parser = argparse.ArgumentParser(add_help=False)
	file_parser = argparse.ArgumentParser(add_help=False)
	table2asn_parser = argparse.ArgumentParser(add_help=False)
	gff_parser = argparse.ArgumentParser(add_help=False)
	test_parser = argparse.ArgumentParser(add_help=False)

	database_parser.add_argument("--biosample", "-b",
		action="store_const",
		const="BIOSAMPLE",
		default="",
		help="Submit to BioSample database.")
	database_parser.add_argument("--sra", "-s",
		action="store_const",
		const="SRA",
		default="",
		help="Submit to SRA database.")
	database_parser.add_argument("--genbank", "-n",
		action="store_const",
		const="GENBANK",
		default="",
		help="Submit to Genbank database.")
	database_parser.add_argument("--gisaid", "-g",
		action="store_const",
		const="GISAID",
		default="",
		help="Submit to GISAID database.")
	organism_parser.add_argument("--organism",
		help="Type of organism data",
		choices=ORGANISM_CHOICES,
		type=str.upper,
		default=ORGANISM_CHOICES[0],
		required=True)
	submission_name_parser.add_argument("--submission_name",
		help='Name of the submission',
		required=True)	
	submission_dir_parser.add_argument("--submission_dir",
		help='Directory to where all required files (such as metadata, fasta, etc.) are stored',
		required=True)	
	config_file_parser.add_argument("--config_file",
		help="Config file stored in submission directory",
		required=True)
	file_parser.add_argument("--metadata_file",
		help="Metadata file stored in submission directory",
		required=True)
	test_parser.add_argument("--test",
		help="Whether to perform a test submission.",
		required=False,
		action="store_const",
		default=False,
		const=True)

	# Parse the database argument
	database_args = database_parser.parse_known_args()[0]

	# If genbank and/or gisaid in the database list, must provide fasta file
	if any(x in database_args for x in ["genbank", "gisaid"]):
		file_parser.add_argument("--fasta_file",
			help="Fasta file stored in submission directory",
			required=True)
	else:
		file_parser.add_argument("--fasta_file",
			help="Fasta file stored in submission directory",
			required=False)	
			
	# If genbank in the database list, determine whether to prepare table2asn submission		
	if any(x in database_args for x in ["genbank"]):
		table2asn_parser.add_argument("--table2asn",
			help="Whether to prepare a Table2asn submission.",
			required=False,
			action="store_const",
			default=False,
			const=True)	
		# Optional: add annotation to table2asn submission
		gff_parser.add_argument("--gff_file",
			help="An annotation file to add to a Table2asn submission",
			required=False)

	# Create the submodule commands
	subparser_modules = parser.add_subparsers(dest='command')

	# prep command
	prep_module = subparser_modules.add_parser(
		'prep',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description='Prepare submission files for future uploads',
		parents=[database_parser, organism_parser, submission_name_parser, submission_dir_parser, config_file_parser, file_parser, table2asn_parser, gff_parser, test_parser]
	)

	# submit command
	submit_module = subparser_modules.add_parser(
		'submit',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description='Create submission files and then batch uploading them to databases of choices.',
		parents=[database_parser, organism_parser, submission_name_parser, submission_dir_parser, config_file_parser, file_parser, table2asn_parser, gff_parser, test_parser]
	)

	# check_submission_status command
	update_module = subparser_modules.add_parser(
		'check_submission_status',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description='Check existing process of a submission', 
		parents=[submission_dir_parser, submission_name_parser, organism_parser, test_parser]
	)

	# template command
	template_module = subparser_modules.add_parser(
		'template',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description='Return a set of files (e.g., config file, metadata file, fasta files, etc.) that are needed to make a submission',
		parents=[database_parser, organism_parser, submission_dir_parser, submission_name_parser]
	)

	# version command
	version_module = subparser_modules.add_parser(
		'version',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description='Show version and exit'
	)

	return(parser, prep_module)

def main():
	"""The main routine"""
	parser, submit_prep_subparser  = args_parser()
	args = parser.parse_args()

	# Parse the command argument
	command = args.command

	# Determine whether required files that needed in the command
	database = []
	if "biosample" in args:
		database += [args.biosample]  
	if "sra" in args:
		database += [args.sra]  
	if "genbank" in args:
		database += [args.genbank]  
	if "gisaid" in args:
		database += [args.gisaid] 	
	if "organism" in args:
		organism = args.organism
	if "submission_name" in args:
		submission_name = args.submission_name
	if "submission_dir" in args:
		submission_dir = args.submission_dir
	if "config_file" in args:
		config_file = args.config_file
	if "metadata_file" in args:
		metadata_file = args.metadata_file

	# Determine if optional arguments are given. If not, set a DEFAULT value to each argument
	# fasta
	if "fasta_file" in args:
		fasta_file = args.fasta_file
	else:
		fasta_file = None
	# table2asn
	if "table2asn" in args:
		table2asn = args.table2asn
	else:
		table2asn = False
	# gff_file
	if "gff_file" in args:
		gff_file = args.gff_file
	else:
		gff_file = None
	# test submission
	if "test" in args:
		test = args.test
	else:
		test = False

	# Get database list	
	database = [x for x in database if x]
	
	# Execute the command
	if command in ["prep", "submit"]:
		# If database is not given, display help
		if len(database) == 0:
			print("\n"+"ERROR: Missing a database selection. See USAGE below."+"\n", file=sys.stdout)
			submit_prep_subparser.print_help()
			sys.exit(0)
		start(command=command, organism=organism, database=database, submission_name=submission_name, submission_dir=submission_dir, config_file=config_file, metadata_file=metadata_file, fasta_file=fasta_file, table2asn=table2asn, gff_file=gff_file, test=test)
	elif command == "check_submission_status":
		process.update_submission_status(submission_dir=submission_dir, submission_name=submission_name, organism=organism, test=test)
	elif command == "template":
		# If database is not given, display help
		if len(database) == 0:
			print("\n"+"ERROR: Missing a database selection. See USAGE below."+"\n", file=sys.stdout)
			submit_prep_subparser.print_help()
			sys.exit(0)
		setup.create_zip_template(organism=organism, database=database, submission_dir=submission_dir, submission_name=submission_name)
	elif command == "version":
		print("\n"+"Version: " + VERSION, file=sys.stdout)
		sys.exit(0)
	else:
		# If no command display help
		parser.print_help()
		sys.exit(0)

	# Print out the execution time	
	get_execution_time()

if __name__ == "__main__":
	main()
