#!/usr/bin/env python3

###########################    Description    ##################################
# Parsers for reading seqsender input
################################################################################

import argparse
from typing import List

def args_parser(organism_choices: List[str]):
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
		choices=organism_choices,
		default="",
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

	# Parse arguments that change inputs
	database_args = database_parser.parse_known_args()[0]

	# If genbank and/or gisaid in the database list, must provide fasta file
	if "GENBANK" in database_args or "GISAID" in database_args:
		file_parser.add_argument("--fasta_file",
			help="Fasta file stored in submission directory",
			required=True)
	else:
		file_parser.add_argument("--fasta_file",
			help="Fasta file stored in submission directory",
			default=None,
			required=False)

	# Optional: add annotation to table2asn submission
	gff_parser.add_argument("--gff_file",
		help="An annotation file to add to a Table2asn submission",
		default=None,
		required=False)

	# Create the submodule commands
	subparser_modules = parser.add_subparsers(dest='command')

	# prep command
	prep_module = subparser_modules.add_parser(
		'prep',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description='Prepare submission files for future uploads',
		parents=[database_parser, organism_parser, submission_name_parser, submission_dir_parser, config_file_parser, file_parser, gff_parser, test_parser]
	)

	# submit command
	submit_module = subparser_modules.add_parser(
		'submit',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description='Create submission files and then batch uploading them to databases of choices.',
		parents=[database_parser, organism_parser, submission_name_parser, submission_dir_parser, config_file_parser, file_parser, gff_parser, test_parser]
	)

	# check_submission_status command
	update_module = subparser_modules.add_parser(
		'submission_status',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description='Check existing process of either all submissions or a single submission if given a name.',
		parents=[submission_dir_parser, submission_name_parser]
	)

	# Generate test data command
	test_output_module = subparser_modules.add_parser(
		'test_data',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description='Return a set of files (e.g., config file, metadata file, fasta files, etc.) that are needed to make a submission for a specific organism/database.',
		parents=[database_parser, organism_parser, submission_dir_parser]
	)

	# biosample xml download command
	biosample_xml_module = subparser_modules.add_parser(
		'update_biosample',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description='Updates the Pandera schema validation for all BioSample packages based off the BioSample Package XML.'
	)

	# version command
	version_module = subparser_modules.add_parser(
		'version',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		description='Show version and exit'
	)

	return(parser, prep_module)
