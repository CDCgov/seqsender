
# Python Libraries
import os
import sys
import importlib
import pathlib
import pandas as pd
from settings import PROG_DIR
from typing import List, Dict, Any, Optional, Union
import file_handler
import json
import pandera
from pandera import Check
from datetime import datetime, timedelta
from cerberus import Validator
import re

from config.seqsender.seqsender_schema import schema as seqsender_schema
from seqsender import logger
from settings import PROG_DIR, SCHEMA_EXCLUSIONS, BIOSAMPLE_REGEX, SRA_REGEX, GISAID_REGEX, GENBANK_REGEX, GENBANK_REGEX_CMT, GENBANK_REGEX_SRC, GENBANK_DEPRECATED_COLUMNS

# Check the config file
def get_config(config_file: str, databases: List[str]) -> Dict[str, Any]:
	logger.info("Loading config file...")
	# Determine required database
	submission_portals = set()
	for database in databases:
		if "BIOSAMPLE" in database or "SRA" in database or "GENBANK" in database:
			submission_portals.add("ncbi")
		if "GISAID" in database:
			submission_portals.add("gisaid")
	# Check if list empty
	if not submission_portals:
		logger.error("Submission portals list cannot be empty.")
		sys.exit(1)
	submission_schema = "_".join(submission_portals)
	# Read in user config file
	config_dict = file_handler.load_yaml(yaml_type = "Config file", yaml_path = config_file)
	# Check if yaml forms dictionary
	logger.debug("Validating config file...")
	if type(config_dict) is dict:
		schema = eval(open(os.path.join(PROG_DIR, "config", "seqsender", "config_file", (submission_schema + "_schema.py")), 'r').read())
		database_specific_config_schema_updates(schema, databases)
		validator = Validator(schema)
		# Validate based on schema
		if validator.validate(config_dict, schema) is False:
			logger.error(f"Config file is not properly setup. Please correct config file based on issue below:\n{json.dumps(validator.errors, indent = 4)}")
			sys.exit(1)
		else:
			if "GENBANK" in databases and "GISAID" in databases:
				validate_submission_position(config_dict=config_dict)
			config_dict = parse_hold_date(config_dict=config_dict)
			logger.success("Config file successfully loaded.")
			return config_dict["Submission"]
	else:
		logger.error("Config file is incorrect. File must be a valid yaml format.")
		sys.exit(1)

def validate_submission_position(config_dict: Dict[str, Any]):
	genbank_position = get_submission_position(config_dict=config_dict, database="GENBANK")
	gisaid_position = get_submission_position(config_dict=config_dict, database="GISAID")
	if (gisaid_position is None and genbank_position is not None) or (gisaid_position is not None and genbank_position is None) or (isinstance(gisaid_position, int) and isinstance(genbank_position, int) and gisaid_position == genbank_position):
		logger.error(f"Config file is incorrect. Submission position for GISAID '{gisaid_position}' and GenBank '{genbank_position}' must both be either left empty, or set to '1' and '2' based on submission preference.")
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
		logger.error(f"Database {database} is not a valid selection.")
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
					logger.error(f"Config file field 'Specified_Release_Date', with value '{config_dict['Submission']['NCBI']['Specified_Release_Date']}' is not a valid time delta. To be a valid timedelta format, it must be formatted as '<num> days', '<num> weeks', or '<num> months'. The prefix <num> is a numeric value which will be used to create a specified release date based on when the SeqSender 'submit' command is ran with the time delta specified added.")
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
					logger.error(f"Config file field 'Specified_Release_Date', with value '{config_dict['Submission']['NCBI']['Specified_Release_Date']}' is not a valid date. To be a valid date format, it must be formatted as 'YYYY-MM-DD', with zero padding and it must be later than today {datetime.now().strftime('%Y-%m-%d')}.")
					sys.exit(1)
		except ValueError:
			pass
		logger.error(f"Unable to parse config file field 'Specified_Release_Date', with value '{config_dict['Submission']['NCBI']['Specified_Release_Date']}'. For field to be valid it must be 'None', formatted as 'YYYY-MM-DD', with zero padding and it must be later than today {datetime.now().strftime('%Y-%m-%d')}.")
		sys.exit(1)
	return config_dict

# Error out if deprecated submission column names detected
def check_deprecated_columns(database: List[str], metadata: pd.DataFrame) -> None:
	logger.debug("Checking for deprecated columns...")
	if "GENBANK" in database:
		deprecated_columns = [col for col in GENBANK_DEPRECATED_COLUMNS if col in metadata.columns]
		if deprecated_columns:
			logger.error(f"GenBank columns '{deprecated_columns}' are deprecated and are no longer supported by GenBank. Please remove them before submission.")
			sys.exit(1)
	logger.debug("No deprecated columns found.")

def import_schema(database: str, schema_name: str, schema_regex: str, custom_validation_dir: str):
	default_schema_path = schema_path = "config." + database.strip().lower() + "." + schema_name
	if custom_validation_dir:
		custom_schema_path = os.path.join(custom_validation_dir, schema_name + ".py")
		logger.debug(f"Checking for custom '{database}' schema file: {custom_schema_path}")
		if os.path.isfile(custom_schema_path):
			logger.debug(f"Found custom '{database}' schema file.")
			schema_path = "config.custom.ID.genbank_cmt_schema"
		else:
			logger.debug(f"Custom database schema for '{database}' not found, reverting to default schema.")
	try:
		logger.debug(f"Attempting to import schema module '{schema_name}' at path: {schema_path}")
		schema_import = (schema_regex, importlib.import_module(schema_path).schema)
		logger.debug("Schema successfully imported.")
	except ModuleNotFoundError:
		logger.error(f"Database {schema_name} schema not found.")
		sys.exit(1)
	except ImportError as e:
		logger.error(f"Issue when importing database schema '{schema_name}':\n{e}")
		sys.exit(1)
	except Exception as e:
		logger.critical(f"An unexpected error occured when importing database schema '{schema_name}':\n{e}")
		sys.exit(1)
		logger.debug("Schema module imported.")
	return schema_import

def load_custom_validation(custom_validation):
	cleanup_sys_path = False
	custom_validation_dir = os.path.join(PROG_DIR, "config", "custom", custom_validation)
	logger.debug(f"Checking custom validation location: '{custom_validation_dir}'")
	if os.path.isdir(custom_validation_dir) == False:
		logger.warn(f"Custom validation '{custom_validation}' enabled but directory '{custom_validation_dir}' not found")
	if custom_validation_dir not in sys.path:
		logger.debug(f"Adding custom validation directory to sys.path")
		sys.path.append(custom_validation_dir)
		cleanup_sys_path = True
	logger.success(f"Custom validation added: '{custom_validation}'")
	return custom_validation_dir, cleanup_sys_path

def load_schemas(metadata: pd.DataFrame, organism: str, database: List[str], config_dict: Dict[str, Any], custom_validation_dir: str) -> Dict[str, Any]:
	# Import schemas
	schemas = dict()
	if "BIOSAMPLE" in database:
		biosample_schema = config_dict["NCBI"]["BioSample_Package"].strip().replace(".", "_")
		schemas["BioSample"] = import_schema(database="BioSample", schema_name=biosample_schema, schema_regex=BIOSAMPLE_REGEX, custom_validation_dir=custom_validation_dir)
	if "SRA" in database:
		schemas["SRA"] = import_schema(database="SRA", schema_name="sra_schema", schema_regex=SRA_REGEX, custom_validation_dir=custom_validation_dir)
	if "GENBANK" in database:
		schemas["GenBank"] = import_schema(database="GenBank", schema_name="genbank_schema", schema_regex=(GENBANK_REGEX + "|^sequence_name$"), custom_validation_dir=custom_validation_dir)
		if [col_name for col_name in metadata if isinstance(col_name, str) and col_name.startswith("cmt-")]:
			schemas["GenBank comment"] = import_schema(database="GenBank", schema_name="genbank_cmt_schema", schema_regex=GENBANK_REGEX_CMT, custom_validation_dir=custom_validation_dir)
		if [col_name for col_name in metadata if isinstance(col_name, str) and col_name.startswith("src-")]:
			src_schema = "genbank_src_schema"
			if organism == "FLU":
				src_schema = "genbank_flu_src_schema"
			schemas["GenBank source"] = import_schema(database="GenBank", schema_name=src_schema, schema_regex=GENBANK_REGEX_SRC, custom_validation_dir=custom_validation_dir)
	if "GISAID" in database:
		schemas["GISAID"] = import_schema(database="GISAID", schema_name=("gisaid_" + organism + "_schema"), schema_regex=(GISAID_REGEX + "|^sequence_name$"), custom_validation_dir=custom_validation_dir)
	return schemas

def validate_metadata(metadata: pd.DataFrame, database: List[str], schemas: Dict[str, Any]) -> List[pandera.errors.SchemaErrors]:
	# Update seqsender base schema to include needed checks
	if "BIOSAMPLE" in database or "SRA" in database:
		seqsender_schema.update_columns({"bioproject":{"checks":Check.str_matches(r"^(?!\s*$).+"),"nullable":False,"required":True}})
	biosample_schema = sra_schema = genbank_schema = genbank_cmt_schema = genbank_src_schema = gisaid_schema = None
	# Validate metadata on schema's
	error_msg_list: List[pandera.errors.SchemaErrors] = []
	# Validate required columns for seqsender
	logger.debug("Validating metadata against schemas...")
	try:
		logger.debug("Validating metadata: SeqSender schema")
		seqsender_schema.validate(metadata, lazy = True)
		logger.debug("Metadata validated.")
	except pandera.errors.SchemaErrors as schema_error:
		error_msg_list.append(schema_error)
		logger.error("Metadata failed SeqSender schema.")
	# Validate required columns for databases
	for schema_name, tuple in schemas.items():
		try:
			logger.debug(f"Validating metadata: {schema_name} schema")
			regex, schema = tuple
			database_specific_metadata = metadata.filter(regex=regex).copy().drop_duplicates()
			schema.validate(database_specific_metadata, lazy = True)
			logger.debug("Metadata validated.")
		except pandera.errors.SchemaErrors as schema_error:
			error_msg_list.append(schema_error)
			logger.error(f"Metadata failed {schema_name} schema.")
	return error_msg_list

def expand_metadata_columns(metadata: pd.DataFrame) -> pd.DataFrame:
	expand_columns = [col for col in metadata.columns if "|" in col]
	for multi_name in expand_columns:
		for col_name in multi_name.split("|"):
			col_name = col_name.strip()
			metadata[col_name] = metadata[multi_name]
	metadata = metadata.drop(columns=expand_columns)
	return metadata

# Read in metadata file
def get_metadata(database: List[str], organism: str, metadata_file: str, config_dict: Dict[str, Any], skip_validation: bool = False, custom_validation: str = "") -> pd.DataFrame:
	# Read in metadata file
	logger.info("Loading metadata file...")
	metadata = file_handler.load_csv(metadata_file)
	metadata = expand_metadata_columns(metadata=metadata)
	if skip_validation == False:
		check_deprecated_columns(database = database, metadata = metadata)
		custom_validation_dir, cleanup_sys_path = "", False
		if custom_validation:
			custom_validation_dir, cleanup_sys_path = load_custom_validation(custom_validation)
		if custom_validation and cleanup_sys_path:
			logger.debug("Removing custom validation directory from sys.path")
			sys.path.remove(custom_validation_dir)
		schemas = load_schemas(metadata=metadata, organism=organism, database=database, config_dict=config_dict, custom_validation_dir=custom_validation_dir)
		validation_errors = validate_metadata(metadata=metadata, database=database, schemas=schemas)
		if validation_errors:
			pretty_print_pandera_errors(file=metadata_file, error_msgs=validation_errors)
			sys.exit(1)
		logger.debug("Metadata validated against schemas.")
	else:
		logger.warning("Metadata validation skipped.")
	logger.success("Metadata successfully loaded.")
	return metadata

def pretty_print_pandera_errors(file: str, error_msgs: List[pandera.errors.SchemaErrors]):
	logger.error(f"File {file} has the following error('s):\n(Note: Index position is calculated excluding column headers and the first row index value starting at '1'.)\n")
	for specific_schema_error in error_msgs:
		duplicate_errors = dict()
		for error in specific_schema_error.failure_cases.itertuples():
			# Missing column
			if error.check == "column_in_dataframe":
				logger.error(f"Missing required column '{error.failure_case}', ensure the file has not been modified and retry.")
			# Column requires specific values capitalization matters
			elif re.search("isin\(\['.*'(, '.*')+\]\)", error.check):
				logger.error(f"Column '{error.column}' has an incorrect value at index '{(error.index + 1)}'. This field can only contain the values '{(error.check.replace('isin(', '')[:-1])}', you provided '{error.failure_case}'.")
			# Column cannot have null values or empty strings
			elif error.check == "str_matches('^(?!\\s*$).+')":
				logger.error(f"Column '{error.column}' has an empty field at index '{(error.index + 1)}' that is required. This field cannot be left blank.")
			# Column requires specific values capitalization does not matter
			elif re.search("str_matches\(\'\(\?i\)\(\\\\W\|\^\)\(.*\|.*\)\(\\\\W\|\$\)\'\)", error.check):
				accepted_values = error.check.replace("str_matches('(?i)(\\W|^)(", "").replace(")(\W|$)')", "").split("|")
				logger.error(f"Column '{error.column}' at index '{(error.index + 1)}' has the value '{error.failure_case}'. This field must be one of the accepted values: {accepted_values}.")
			# Column submission group, among a group of columns at least one must contain a non null value
			elif re.search("\(lambda df: ~\(df\[\".*\"\].isnull\(\)( & df\[\".*\"\].isnull\(\))+\), ignore_na = False\)", error.check):
				column_group = error.check.replace("(lambda df: ~(df[\"", "").replace("\"].isnull() & df[\"", "+").replace("\"].isnull()), ignore_na = False)", "").split("+")
				logger.error(f"In column group, every sample must have at least one non-null value in at least one of the following columns: '{column_group}'.")
			# Column has minimum or maximum character length requirements
			elif re.search("str_length\(.*\)", error.check):
				min_value, max_value = error.check.replace("str_length(", "").replace(")", "").split(", ")
				logger.error(f"Column '{error.column}' has a character limit of minimum '{min_value}', maximum '{max_value}'. The value '{error.failure_case}' at index '{(error.index + 1)}' does not meet these requirements.")
			# Unique submission_log.csv column that takes capitalization controlled "PENDING" and "SUBMITTED"
			# but also takes in a raw value from the NCBI report file generated and accounts for possible whitespace
			elif error.check == "str_matches('^(PENDING|SUBMITTED|\WSUB\d*\W)$')":
				accepted_values = error.check.replace("str_matches('^(", "").replace("\W", "").replace(")$')", "").replace("\d*", "<numeric_values>").split("|")
				logger.error(f"Column '{error.column}' at index '{(error.index + 1)}' has a value '{error.failure_case}'. This field must be one of the accepted values: '{accepted_values}'.")
			# Dates in column are incorrectly formatted
			elif error.check == "invalid_date_format":
				logger.error(f"Column '{error.column}' at index '{(error.index + 1)}' has a value '{error.failure_case}'. This field must be a valid date format based on ISO 8601: '[\"YYYY-MM-DD\", \"YYYY-MM\", or \"YYYY\"]'.")
			# Check columns that must have identical value (i.e. NCBI GUI submission portal title)
			elif error.schema_context.lower() == "column" and error.column in ["bs-title", "bs-comment", "sra-title", "sra-comment", "gb-title", "gb-comment"]:
				logger.error(f"Column '{error.column}' must have the same value for every row as it is only used once and applies to the entire submission. This field is an internal NCBI field for the NCBI submission portal website (https://submit.ncbi.nlm.nih.gov/subs/) to aid you in identifying your submissions.")
			# Check ordered columns
			elif error.check == "column_ordered":
				logger.error(f"Column '{error.failure_case}' is incorrectly ordered for file '{file}'.")
			# Check sra file names
			elif error.check == "no_regex_column_match('sra-file_[2-9]\d*')":
				logger.error("Column 'sra-file_#' is required, where # is the numeric value of the file for the SRA sample. (i.e. sra-file_1)")
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
				logger.error(f"Unable to pretty print pandera error message for data. Error message is:\n{error}")
				logger.error(f"{error.schema_context}: {error.column}")
				if error.index:
					logger.error(f"Index: {(error.index + 1)}")
				logger.error(f"Value: {error.failure_case}")
				logger.error(f"Validator: {error.check}")
				logger.error("If you would like to contribute to SeqSender. Make a issue on github reporting this error case to have a descriptive version of this error added.\n")
		if duplicate_errors:
			for column_value_msg, indices in duplicate_errors.items():
				logger.error(f"{column_value_msg} is duplicated at indices: '{indices}'.\n")

# Check user credentials information
def check_credentials(config_dict: Dict[str, Any], database: str) -> None:
	# Check username
	logger.debug(f"Checking {database} credentials...")
	if "Username" not in config_dict.keys():
		logger.error("There is no Submission > " + database + " > Username information in config file.")
		sys.exit(1)
	elif ("Username" in config_dict.keys()) and ((config_dict["Username"] is not None) and (config_dict["Username"] != "")):
		pass
	else:
		logger.error("Submission > " + database + " > Username in the config file cannot be empty.")
		sys.exit(1)
	# Check password
	if "Password" not in config_dict.keys():
		logger.error("There is no Submission > " + database + " > Password information in config file.")
		sys.exit(1)
	elif ("Password" in config_dict.keys()) and ((config_dict["Password"] is not None) and (config_dict["Password"] != "")):
		pass
	else:
		logger.error("Submission > " + database + " > Password in the config file cannot be empty.")
		sys.exit(1)
	# Check client-id if database is GISAID
	if database != "GISAID":
		return
	elif "Client-Id" not in config_dict.keys():
		logger.error("There is no Submission > " + database + " > Client-Id information in config file.")
		sys.exit(1)
	elif ("Client-Id" in config_dict.keys() and ((config_dict["Client-Id"] is not None) and (config_dict["Client-Id"] != ""))):
		pass
	else:
		logger.error("Submission > " + database + " > Client-Id in the config file cannot be empty.")
		sys.exit(1)
	logger.debug("Credentials found.")

# Check sample names in metadata file are listed in fasta file
def process_fasta_samples(metadata: pd.DataFrame, fasta_file: str) -> pd.DataFrame:
	logger.info("Loading fasta file...")
	fasta_df = file_handler.load_fasta_file(fasta_file)
	# Check duplicates in fasta_df
	logger.debug("Checking for duplicates in fasta...")
	duplicated_df = fasta_df[fasta_df.duplicated(subset = ["fasta_name_orig"], keep = False)]
	if not duplicated_df.empty:
		logger.error(f"Sequences in fasta file must be unique at: {fasta_file}")
		logger.error(f"Duplicate Sequences\n{duplicated_df['fasta_sequence_orig'].to_string(index=False)}")
		sys.exit(1)
	logger.debug("No duplicates in fasta.")
	# Validate duplicates don't appear on merge
	logger.debug("Merging fasta with metadata...")
	try:
		merged_df = metadata.merge(fasta_df, how = "outer", left_on = "sequence_name", right_on = "fasta_name_orig", validate = "1:1")
	except:
		logger.error("Unable to merge fasta file to metadata file. Please validate there are not duplicate sequences in both files.")
		sys.exit(1)
	# Check if fasta file has sequences not in metadata
	if merged_df["sequence_name"].isnull().any():
		logger.error("Sequences in fasta file do not have an associated sequence in metadata file. Please update sequences below:\n" + merged_df[merged_df["sequence_name"].isnull()]["fasta_name_orig"].to_string())
		sys.exit(1)
	# Check if metadata has sequences not in fasta file
	if merged_df["fasta_name_orig"].isnull().any():
		logger.error("Sequences in metadata file do not have an associated sequence in fasta file. Please update sequences below:\n" + merged_df[merged_df["fasta_name_orig"].isnull()]["sequence_name"].to_string())
		sys.exit(1)
	logger.debug("Fasta successfully merged with metadata.")
	logger.success("Fasta successfully loaded.")
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
			logger.warning("Unable to load schema \"" + schema_name + "\".")
			logger.warning(e)
			continue
		template = dict()
		try:
			metadata_template = process_schema(schema)
		except Exception as e:
			logger.warning("Unable to process schema into metadata template.")
			logger.warning(e)
			continue
		metadata_template.to_csv(os.path.join(PROG_DIR, "shiny", "templates", schema_name.replace("_", ".") + "_template.csv"), header = True, index = False)
