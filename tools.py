
# Python Libraries
import os
import sys
import importlib
import pathlib
import pandas as pd
from settings import PROG_DIR
from typing import List, Dict, Any, Optional, Union, Set
import file_handler
import json
import pandera
from pandera import Check
from datetime import datetime, timedelta
from cerberus import Validator
import re

from config.seqsender.seqsender_schema import schema as seqsender_schema
from settings import SCHEMA_EXCLUSIONS, BIOSAMPLE_REGEX, SRA_REGEX, GISAID_REGEX, GENBANK_REGEX, GENBANK_REGEX_CMT, GENBANK_REGEX_SRC, GENBANK_DEPRECATED_COLUMNS

# Check the config file
def get_config(config_file: str, databases: List[str]) -> Dict[str, Any]:
	# Determine required database
	submission_portals = set()
	for database in databases:
		if "BIOSAMPLE" in database or "SRA" in database or "GENBANK" in database:
			submission_portals.add("ncbi")
		if "GISAID" in database:
			submission_portals.add("gisaid")
	# Check if list empty
	if not submission_portals:
		print("Error: Submission portals list cannot be empty.", file=sys.stderr)
		sys.exit(1)
	submission_schema_file = get_submission_schema_config_name(submission_portals=submission_portals)
	# Read in user config file
	config_dict = file_handler.load_yaml(yaml_type = "Config file", yaml_path = config_file)
	# Check if yaml forms dictionary
	if type(config_dict) is dict:
		schema = eval(open(os.path.join(PROG_DIR, "config", "seqsender", "config_file", submission_schema_file), 'r').read())
		database_specific_config_schema_updates(schema, databases)
		validator = Validator(schema)
		# Validate based on schema
		if validator.validate(config_dict, schema) is False:
			print("Error: Config file is not properly setup. Please correct config file based on issue below:", file=sys.stderr)
			print(json.dumps(validator.errors, indent = 4), file=sys.stderr)
			sys.exit(1)
		else:
			if "GENBANK" in databases and "GISAID" in databases:
				validate_submission_position(config_dict=config_dict)
			config_dict = parse_hold_date(config_dict=config_dict)
			return config_dict["Submission"]
	else:
		print("Error: Config file is incorrect. File must be a valid yaml format.", file=sys.stderr)
		sys.exit(1)

def get_submission_schema_config_name(submission_portals: Set[str]) -> str:
	submission_schema_file_name = ""
	if "ncbi" in submission_portals:
		submission_schema_file_name += "ncbi_"
	if "gisaid" in submission_portals:
		submission_schema_file_name += "gisaid_"
	submission_schema_file_name += "schema.py"
	return submission_schema_file_name

def validate_submission_position(config_dict: Dict[str, Any]):
	genbank_position = get_submission_position(config_dict=config_dict, database="GENBANK")
	gisaid_position = get_submission_position(config_dict=config_dict, database="GISAID")
	if (gisaid_position is None and genbank_position is not None) or (gisaid_position is not None and genbank_position is None) or (isinstance(gisaid_position, int) and isinstance(genbank_position, int) and gisaid_position == genbank_position):
		print(f"Error: Config file is incorrect. Submission position for GISAID '{gisaid_position}' and GenBank '{genbank_position}' must both be either left empty, or set to '1' and '2' based on submission preference.", file=sys.stderr)
		sys.exit(1)

def get_submission_type(test: bool) -> str:
	if test:
		return "TEST"
	else:
		return "PRODUCTION"

def get_submission_position(config_dict: Dict[str, Any], database: str) -> Optional[int]:
	if database in ["BIOSAMPLE", "SRA", "GENBANK"]:
		parent_database = "NCBI"
	elif database == "GISAID":
		parent_database = "GISAID"
	else:
		print(f"Error: database {database} is not a valid selection.", file=sys.stderr)
		sys.exit(1)
	if "Submission" in config_dict:
		config_dict = config_dict["Submission"]
	if parent_database in config_dict:
		config_dict = config_dict[parent_database]
	if "Submission_Position" in config_dict and isinstance(config_dict["Submission_Position"], int):
		return config_dict["Submission_Position"]
	else:
		return None

def database_specific_config_schema_updates(schema: Dict[str, Any], database: List[str]) -> Dict[str, Any]:
	# Update seqsender base schema to include needed checks
	if "BIOSAMPLE" in database:
		schema["Submission"]["schema"]["NCBI"]["schema"]["BioSample_Package"]["required"] = True
		schema["Submission"]["schema"]["NCBI"]["schema"]["BioSample_Package"]["nullable"] = False
	if "GENBANK" in database:
		schema["Submission"]["schema"]["NCBI"]["schema"]["Publication_Title"]["required"] = True
		schema["Submission"]["schema"]["NCBI"]["schema"]["Publication_Title"]["nullable"] = False
		schema["Submission"]["schema"]["NCBI"]["schema"]["Publication_Status"]["required"] = True
		schema["Submission"]["schema"]["NCBI"]["schema"]["Publication_Status"]["nullable"] = False
	return schema

# Parse Config file specified release date for NCBI field
def parse_hold_date(config_dict: Dict[str, Any]):
	if "NCBI" in config_dict["Submission"] and "Specified_Release_Date" in config_dict["Submission"]["NCBI"] and config_dict["Submission"]["NCBI"]["Specified_Release_Date"] and config_dict["Submission"]["NCBI"]["Specified_Release_Date"].strip() != "":
		release_date_string = config_dict["Submission"]["NCBI"]["Specified_Release_Date"].strip().lower()
		try:
			if re.search(r"\d+\s*(days|weeks|years)", release_date_string):
				today = pd.Timestamp.now().date()
				numeric_value = int(release_date_string.strip().split(" ")[0])
				if "days" in release_date_string:
					time_delta = pd.DateOffset(days = numeric_value)
				elif "weeks" in release_date_string:
					time_delta = pd.DateOffset(weeks = numeric_value)
				elif "months" in release_date_string:
					time_delta = pd.DateOffset(months = numeric_value)
				else:
					print(f"Error: Config file field 'Specified_Release_Date', with value '{config_dict['Submission']['NCBI']['Specified_Release_Date']}' is not a valid time delta. To be a valid timedelta format, it must be formatted as '<num> days', '<num> weeks', or '<num> months'. The prefix <num> is a numeric value which will be used to create a specified release date based on when the SeqSender 'submit' command is ran with the time delta specified added.", file=sys.stderr)
					sys.exit(1)
				calculated_date = today + time_delta
				config_dict["Submission"]["NCBI"]["Specified_Release_Date"] = calculated_date.strftime("%Y-%m-%d")
				return config_dict
			elif re.search(r"\d{4}-\d{2}-\d{2}", release_date_string):
				parsed_date = datetime.strptime(release_date_string, "%Y-%m-%d")
				if parsed_date.date() > datetime.today().date():
					config_dict["Submission"]["NCBI"]["Specified_Release_Date"] = parsed_date.strftime("%Y-%m-%d")
					return config_dict
				else:
					print(f"Error: Config file field 'Specified_Release_Date', with value '{config_dict['Submission']['NCBI']['Specified_Release_Date']}' is not a valid date. To be a valid date format, it must be formatted as 'YYYY-MM-DD', with zero padding and it must be later than today {datetime.now().strftime('%Y-%m-%d')}.", file=sys.stderr)
					sys.exit(1)
		except ValueError:
			pass
		print(f"Error: Unable to parse config file field 'Specified_Release_Date', with value '{config_dict['Submission']['NCBI']['Specified_Release_Date']}'. For field to be valid it must be 'None', formatted as 'YYYY-MM-DD', with zero padding and it must be later than today {datetime.now().strftime('%Y-%m-%d')}.", file=sys.stderr)
		sys.exit(1)
	return config_dict

# Error out if deprecated submission column names detected
def warn_deprecated_columns(database: List[str], metadata: pd.DataFrame) -> None:
	if "GENBANK" in database:
		deprecated_columns = [col for col in GENBANK_DEPRECATED_COLUMNS if col in metadata.columns]
		if deprecated_columns:
			print(f"Error: GenBank columns '{deprecated_columns}' are deprecated and are no longer supported by GenBank. Please remove them before submission.", file=sys.stderr)
			sys.exit(1)

# Read in metadata file
def get_metadata(database: List[str], organism: str, metadata_file: str, config_dict: Dict[str, Any], skip_validation: bool = False) -> pd.DataFrame:
	# Read in metadata file
	metadata = file_handler.load_csv(metadata_file)
	warn_deprecated_columns(database = database, metadata = metadata)
	# Update seqsender base schema to include needed checks
	if "BIOSAMPLE" in database or "SRA" in database:
		seqsender_schema.update_columns({"bioproject":{"checks":Check.str_matches(r"^(?!\s*$).+"),"nullable":False,"required":True}})
	biosample_schema = sra_schema = genbank_schema = genbank_cmt_schema = genbank_src_schema = gisaid_schema = None
	# Import schemas
	schemas_dict = dict()
	if "BIOSAMPLE" in database:
		schemas_dict["BioSample"] = (BIOSAMPLE_REGEX, importlib.import_module("config.biosample." + config_dict["NCBI"]["BioSample_Package"].strip().replace(".", "_")).schema)
	if "SRA" in database:
		schemas_dict["SRA"] = (SRA_REGEX, importlib.import_module("config.sra.sra_schema").schema)
	if "GENBANK" in database:
		schemas_dict["GenBank"] = ((GENBANK_REGEX + "|^sequence_name$"), importlib.import_module("config.genbank.genbank_schema").schema)
		# if [col_name for col_name in metadata if col_name.startswith("cmt-")]:
		# 	schemas_dict["GenBank comment"] = (GENBANK_REGEX_CMT, importlib.import_module("config.genbank.genbank_cmt_schema").schema)
		if [col_name for col_name in metadata if isinstance(col_name, str) and col_name.startswith("src-")]:
			if organism == "FLU":
				schemas_dict["GenBank source"] = (GENBANK_REGEX_SRC, importlib.import_module("config.genbank.genbank_flu_src_schema").schema)
			else:
				schemas_dict["GenBank source"] = (GENBANK_REGEX_SRC, importlib.import_module("config.genbank.genbank_src_schema").schema)
	if "GISAID" in database:
		schemas_dict["GISAID"] = ((GISAID_REGEX + "|^sequence_name$"), importlib.import_module("config.gisaid.gisaid_" + organism + "_schema").schema)
	if skip_validation == False:
		# Validate metadata on schema's
		error_msg_list: List[pandera.errors.SchemaErrors] = []
		# Validate required columns for seqsender
		try:
			seqsender_schema.validate(metadata, lazy = True)
		except pandera.errors.SchemaErrors as schema_error:
			error_msg_list.append(schema_error)
		# Validate required columns for databases
		for schema_name, tuple in schemas_dict.items():
			try:
				regex, schema = tuple
				database_specific_metadata = metadata.filter(regex=regex).copy().drop_duplicates()
				schema.validate(database_specific_metadata, lazy = True)
			except pandera.errors.SchemaErrors as schema_error:
				error_msg_list.append(schema_error)
		if error_msg_list:
			pretty_print_pandera_errors(file=metadata_file, error_msgs=error_msg_list)
			sys.exit(1)
	return metadata

def pretty_print_pandera_errors(file: str, error_msgs: List[pandera.errors.SchemaErrors]):
	print(f"Error: file {file} has the following error('s):\n(Note: Index position is calculated excluding column headers and the first row index value starting at '1'.)\n")
	for specific_schema_error in error_msgs:
		duplicate_errors = dict()
		for error in specific_schema_error.failure_cases.itertuples():
			# Missing column
			if error.check == "column_in_dataframe":
				print(f"Error: Missing required column '{error.failure_case}', ensure the file has not been modified and retry.", file=sys.stderr)
			# Column requires specific values capitalization matters
			elif re.search("isin\(\['.*'(, '.*')+\]\)", error.check):
				print(f"Error: Column '{error.column}' has an incorrect value at index '{(error.index + 1)}'. This field can only contain the values '{(error.check.replace('isin(', '')[:-1])}', you provided '{error.failure_case}'.", file=sys.stderr)
			# Column cannot have null values or empty strings
			elif error.check == "str_matches('^(?!\\s*$).+')":
				print(f"Error: Column '{error.column}' has an empty field at index '{(error.index + 1)}' that is required. This field cannot be left blank.", file=sys.stderr)
			# Column requires specific values capitalization does not matter
			elif re.search("str_matches\(\'\(\?i\)\(\\\\W\|\^\)\(.*\|.*\)\(\\\\W\|\$\)\'\)", error.check):
				accepted_values = error.check.replace("str_matches('(?i)(\\W|^)(", "").replace(")(\W|$)')", "").split("|")
				print(f"Error: Column '{error.column}' at index '{(error.index + 1)}' has the value '{error.failure_case}'. This field must be one of the accepted values: {accepted_values}.", file=sys.stderr)
			# Column submission group, among a group of columns at least one must contain a non null value
			elif re.search("\(lambda df: ~\(df\[\".*\"\].isnull\(\)( & df\[\".*\"\].isnull\(\))+\), ignore_na = False\)", error.check):
				column_group = error.check.replace("(lambda df: ~(df[\"", "").replace("\"].isnull() & df[\"", "+").replace("\"].isnull()), ignore_na = False)", "").split("+")
				print(f"Error: In column group, every sample must have at least one non-null value in at least one of the following columns: '{column_group}'.", file=sys.stderr)
			# Column has minimum or maximum character length requirements
			elif re.search("str_length\(.*\)", error.check):
				min_value, max_value = error.check.replace("str_length(", "").replace(")", "").split(", ")
				print(f"Error: Column '{error.column}' has a character limit of minimum '{min_value}', maximum '{max_value}'. The value '{error.failure_case}' at index '{(error.index + 1)}' does not meet these requirements.", file=sys.stderr)
			# Unique submission_log.csv column that takes capitalization controlled "PENDING" and "SUBMITTED"
			# but also takes in a raw value from the NCBI report file generated and accounts for possible whitespace
			elif error.check == "str_matches('^(PENDING|SUBMITTED|\WSUB\d*\W)$')":
				accepted_values = error.check.replace("str_matches('^(", "").replace("\W", "").replace(")$')", "").replace("\d*", "<numeric_values>").split("|")
				print(f"Error: Column '{error.column}' at index '{(error.index + 1)}' has a value '{error.failure_case}'. This field must be one of the accepted values: '{accepted_values}'.", file=sys.stderr)
			# Dates in column are incorrectly formatted
			elif error.check == "invalid_date_format":
				print(f"Error: Column '{error.column}' at index '{(error.index + 1)}' has a value '{error.failure_case}'. This field must be a valid date format based on ISO 8601: '[\"YYYY-MM-DD\", \"YYYY-MM\", or \"YYYY\"]'.", file=sys.stderr)
			# Check columns that must have identical value (i.e. NCBI GUI submission portal title)
			elif error.schema_context.lower() == "column" and error.column in ["bs-title", "bs-comment", "sra-title", "sra-comment", "gb-title", "gb-comment"]:
				print(f"Error: Column '{error.column}' must have the same value for every row as it is only used once and applies to the entire submission. This field is an internal NCBI field for the NCBI submission portal website (https://submit.ncbi.nlm.nih.gov/subs/) to aid you in identifying your submissions.", file=sys.stderr)
			# Check ordered columns
			elif error.check == "column_ordered":
				print(f"Error: Column '{error.failure_case}' is incorrectly ordered for file '{file}'.", file=sys.stderr)
			# Check sra file names
			elif error.check == "no_regex_column_match('sra-file_[2-9]\d*')":
				print("Error: Column 'sra-file_#' is required, where # is the numeric value of the file for the SRA sample. (i.e. sra-file_1)", file=sys.stderr)
			# Collect all duplicate values and print them at the end to group index positions together
			elif error.check == "field_uniqueness":
				column_key = f"Column '{error.column}' with value '{error.failure_case}'"
				if column_key not in duplicate_errors:
					duplicate_errors[column_key] = [(error.index + 1)]
				else:
					duplicate_errors[column_key].append((error.index + 1))
				continue
			# If unable to parse pandera error message
			else:
				print(f"Error: Unable to pretty print pandera error message for data. Error message is:", file=sys.stderr)
				print(error)
				print(f"{error.schema_context}: {error.column}", file=sys.stderr)
				if error.index:
					print(f"Index: {(error.index + 1)}", file=sys.stderr)
				print(f"Value: {error.failure_case}", file=sys.stderr)
				print(f"Validator: {error.check}", file=sys.stderr)
				print("If you would like to contribute to SeqSender. Make a issue on github reporting this error case to have a descriptive version of this error added.", file=sys.stderr)
			print("", file=sys.stderr)
		if duplicate_errors:
			for column_value_msg, indices in duplicate_errors.items():
				print(f"Error: {column_value_msg} is duplicated at indices: '{indices}'.", file=sys.stderr)
				print("", file=sys.stderr)

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

# Collect all schema files from config subdirectory
# Convert file name to importable string
def get_all_schema_files():
	# Load all schema files
	files = os.walk(os.path.join(PROG_DIR, "config"))
	schemas = []
	for root, dir_list, file_list in files:
		# Get seqsender subdirectory structure
		subdirectory = root.split("/seqsender/config")[-1].replace("/", ".")
		root_end = "config."
		if subdirectory and subdirectory != ".":
			root_end = "config" + subdirectory + "."
		# Collect schemas
		schemas += [(root_end + file.replace(".py", "")) for file in file_list if file.endswith(".py")]
	return schemas

# Import pandera schema as module
def load_schema(schema_name: str):
	return importlib.import_module(schema_name).schema

# Process schema file to retrieve field information
def process_schema(schema):
	schema_contents = []
	for column_name, context in schema.columns.items():
		required_field = context.required
		description_field = context.description
		# Convert required boolean to string
		if required_field:
			required_field = "Required"
		else:
			required_field = "Optional"
		# Update required field if required group of columns
		if description_field and "At least one required: Group" in description_field:
			required_field = "At least one field required. Group: " + description_field.split("Group: \"")[-1].split("\".")[0]
		# Update SRA wildcard field for raw files
		if column_name == "sra-file_[2-9]\d*":
			column_name = "sra-file_#"
		schema_contents.append({"column_name": column_name, "required_column": required_field, "description": description_field})
	return pd.DataFrame(schema_contents)

# Update all shiny templates based on all schema files in config folder
def update_all_schema_templates():
	schemas_list = get_all_schema_files()
	for schema_name in schemas_list:
		# Skip schemas in exclusion list
		if schema_name in SCHEMA_EXCLUSIONS:
			continue
		try:
			schema = load_schema(schema_name)
		except Exception as e:
			print("Warning: Unable to load schema \"" + schema_name + "\".", file=sys.stderr)
			print(e, file=sys.stderr)
			continue
		template = dict()
		try:
			metadata_template = process_schema(schema)
		except:
			print("Warning: Unable to process schema into metadata template.", file=sys.stderr)
			print(e, file=sys.stderr)
			continue
		metadata_template.to_csv(os.path.join(PROG_DIR, "shiny", "templates", schema_name.replace("_", ".") + "_template.csv"), header = True, index = False)
