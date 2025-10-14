
# Python Libraries
import pathlib
import sys
from zipfile import ZipFile
import ftplib
import io
import os
import json
import subprocess
import socket
import pandas as pd
import shutil
import platform
import urllib
import gzip
import stat
import requests
import time
import xmltodict
import xml.etree.ElementTree as ET
from io import BytesIO
from typing import List, Dict, Any

# Local imports
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import tools
from settings import NCBI_FTP_HOST

# Get program directory
PROG_DIR: str = os.path.dirname(os.path.abspath(__file__))
# BioSample atribute html prefix
BIOSAMPLE_HTML_PREFIX: str = "https://www.ncbi.nlm.nih.gov/biosample/docs/packages"
# BioSample atribute html suffix
BIOSAMPLE_HTML_SUFFIX: str = "/?format=xml"
# Schema file header
SCHEMA_HEADER: str = r"""from pandera import DataFrameSchema, Column, Check, Index, MultiIndex

schema = DataFrameSchema(
	columns={
		\"bs-sample_name\": Column(
			dtype=\"object\",
			checks=[
				Check.str_matches(r"^(?!\s*$).+"),
			],
			nullable=False,
			unique=True,
			coerce=False,
			required=True,
			description=\"Identifier name used for BioSample. Max length is 50 characters.\",
			title=\"sample_name\",
		),
		\"bs-sample_title\": Column(
			dtype=\"object\",
			checks=[
				Check.str_matches(r"^(?!\s*$).+"),
			],
			nullable=False,
			unique=False,
			coerce=False,
			required=True,
			description=\"Descriptive title for sample.\",
			title=\"sample title\",
		),
		\"bs-sample_description\": Column(
			dtype=\"object\",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description=\"Optional description for sample.\",
			title=\"sample description\",
		),"""
SCHEMA_COLUMNS_FOOTER: str = """
		\"bs-title\": Column(
			dtype=\"object\",
			checks=[
				Check(lambda s: s.nunique() == 1),
			],
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description=\"Optional internal field for how the BioSample submission should be named when viewed from the NCBI submission portal. If not provided, when performing submissions <--submission_name> with the suffix \\\"-BS\\\" will be used instead.\",
			title=\"biosample submission portal name\",
		),
		\"bs-comment\": Column(
			dtype=\"object\",
			checks=[
				Check(lambda s: s.nunique() == 1),
			],
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description=\"Optional internal field explaining the purpose of the submission for when interacting and resolving submission issues with NCBI.\",
			title=\"biosample submission portal description\",
		)"""

TEST_CONNECTIONS = {"HTTP": {"website":"http://www.google.com", "database": "GENERAL", "error_msg": "Possible internet connectivity issues; unable to connect to 'http://www.google.com'."},
"HTTPS": {"website": "https://www.google.com", "database": "GENERAL", "error_msg": "Possible internet connectivity issues; unable to connect to 'https://www.google.com'."},
"NCBI": {"website": "https://www.ncbi.nlm.nih.gov", "database": "NCBI", "error_msg": "Unable to connect to 'https://www.ncbi.nlm.nih.gov'; ensure NCBI services are running and you are able to connect to them before proceeding."},
"NCBI API": {"website": "https://submit.ncbi.nlm.nih.gov", "database": "NCBI", "error_msg": "Unable to connect to 'https://submit.ncbi.nlm.nih.gov'; ensure NCBI services are running and you are able to connect to them before proceeding."},
"GISAID": {"website": "https://www.epicov.org/epi3/start", "database": "GISAID", "error_msg": "Unable to connect to 'https://www.epicov.org/epi3'; ensure GISAID services are running and you are able to connect to them before proceeding."},
"GISAID": {"website": "https://gisaid.org/", "database": "GISAID", "error_msg": "Unable to connect to 'https://www.epicov.org/epi3'; ensure GISAID services are running and you are able to connect to them before proceeding."}
}

# Create example data for testing
def create_test_data(organism: str, database: List[str], submission_dir: str) -> None:
	if organism not in ["FLU", "COV"]:
		print("SeqSender currently only has test data available for the organisms \"FLU\" and \"COV\" currently, more test sets will be added with later versions. ", file=sys.stdout)
		sys.exit(0)
	# Create output directory
	submission_dir = os.path.abspath(submission_dir)
	out_dir = os.path.join(submission_dir, organism + "_TEST_DATA")
	os.makedirs(out_dir, exist_ok = True)
	# Create sra directory
	out_sra_dir = os.path.join(out_dir, "raw_reads")
	# Create a list of files to output
	out_metadata_file = os.path.join(out_dir, "metadata.csv")
	out_config_file = os.path.join(out_dir, "config.yaml")
	out_sequence_file = os.path.join(out_dir, "sequence.fasta")
	out_fastq_1_r1_file = os.path.join(out_sra_dir, "fastq_1_R1.fastq.gz")
	out_fastq_1_r2_file = os.path.join(out_sra_dir, "fastq_1_R2.fastq.gz")
	out_fastq_2_r1_file = os.path.join(out_sra_dir, "fastq_2_R1.fastq.gz")
	out_fastq_2_r2_file = os.path.join(out_sra_dir, "fastq_2_R2.fastq.gz")
	# Create a list of test files to output
	temp_config_file = os.path.join(PROG_DIR, "test_data", organism, organism.lower()+"_config.yaml")
	temp_sequence_file = os.path.join(PROG_DIR, "test_data", organism, organism.lower()+"_sequence.fasta")
	temp_fastq_1_r1_file = os.path.join(PROG_DIR, "test_data", organism, organism.lower()+"_fastq_1_R1.fastq.gz")
	temp_fastq_1_r2_file = os.path.join(PROG_DIR, "test_data", organism, organism.lower()+"_fastq_1_R2.fastq.gz")
	temp_fastq_2_r1_file = os.path.join(PROG_DIR, "test_data", organism, organism.lower()+"_fastq_2_R1.fastq.gz")
	temp_fastq_2_r2_file = os.path.join(PROG_DIR, "test_data", organism, organism.lower()+"_fastq_2_R2.fastq.gz")
	# Print generating message
	print("\n"+"Generating submission test_data", file=sys.stdout)
	# Get combined metadata for all given databases
	database_prefix = {"GENBANK": "gb-", "GISAID": "gs-", "SRA": "sra-", "BIOSAMPLE": "bs-"}
	repeat_columns = ["sample_name", "sequence_name", "collection_date", "organism", "authors", "bioproject", "bs-sample_name"]
	for i in range(len(database)):
		df = pd.read_csv(os.path.join(PROG_DIR, "test_data", organism, organism.lower()+"_"+database[i].lower()+"_metadata.csv"), header = 0, dtype = str, engine = "python", encoding="utf-8", index_col=False, na_filter=False)
		if i == 0:
			combined_metadata = df
			left_match = database_prefix[database[i]] + "sample_name"
		else:
			df = df.drop(columns = [col for col in repeat_columns if col in combined_metadata.columns and col in df.columns])
			combined_metadata = pd.merge(combined_metadata, df, how="left", left_index = True, right_index = True)
	# Write metadata to output directory
	combined_metadata.to_csv(out_metadata_file, index = False)
    # Write config file to output directory
	shutil.copy(temp_config_file, out_config_file)
    # Write fasta file to output directory
	if any([x in ["GENBANK", "GISAID"] for x in database]):
		shutil.copy(temp_sequence_file, out_sequence_file)
    # Write raw reads file to output directory
	if "SRA" in database:
		os.makedirs(out_sra_dir, exist_ok = True)
		shutil.copy(temp_fastq_1_r1_file, out_fastq_1_r1_file)
		shutil.copy(temp_fastq_1_r2_file, out_fastq_1_r2_file)
		shutil.copy(temp_fastq_2_r1_file, out_fastq_2_r1_file)
		shutil.copy(temp_fastq_2_r2_file, out_fastq_2_r2_file)
	print("Files are stored at: "+os.path.join(out_dir), file=sys.stdout)

def download_table2asn(table2asn_dir: str) -> None:
	# Determine which platform to download table2asn
	if platform.system() == "Windows":
		zip_url = "https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/win64.table2asn.zip"
		with urllib.request.urlopen(zip_url) as zip_response:
			with ZipFile(BytesIO(zip_response.read())) as zip_file:
				zip_file.extractall(table2asn_dir)
		return
	elif platform.system() == "Darwin":
		zip_url = "https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/mac.table2asn.gz"
	elif platform.system() == "Linux":
		zip_url = "https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/linux64.table2asn.gz"
	else:
		print("Error: Cannot identify correct system platform. Please download correct Table2asn version for system and place it in script directory.", file=sys.stderr)
		sys.exit(1)
	# Extract table2asn to tmp folder
	try:
		with open(table2asn_dir, "wb") as file:
			with urllib.request.urlopen(zip_url) as zip_response:
				file.write(gzip.decompress(zip_response.read()))
		st = os.stat(table2asn_dir)
		os.chmod(table2asn_dir, st.st_mode | stat.S_IXOTH | stat.S_IRWXU)
	except Exception as error:
		print("Downloading table2asn error", file=sys.stderr)
		print(error, file=sys.stderr)
		sys.exit(1)

# Download xml and write to a file
def download_xml(xml_url: str, output_file: str) -> None:
	r = requests.get(xml_url)
	r.encoding = "UTF-8"
	with open(output_file, "w+") as file:
		file.write(r.text.replace("\xa0", " "))

# Download list of BioSample packages then download the xml for each package
def download_biosample_xml_list() -> None:
	# Download list of all packages
	download_xml(xml_url = (BIOSAMPLE_HTML_PREFIX + BIOSAMPLE_HTML_SUFFIX), output_file = os.path.join(PROG_DIR, "config", "biosample", "biosample_package_list.xml"))
	with open(os.path.join(PROG_DIR, "config", "biosample", "biosample_package_list.xml")) as file:
		line = file.readline()
		i = 0
		while line:
			if "<Name>" in line:
				name = line.replace("<Name>","").replace("</Name>","").strip()
				# Skip hidden template xml on NCBI website
				if "Generic.1.0" in name:
					line = file.readline()
					continue
				print("Downloading Package: " + name)
				try:
					download_xml(xml_url = (BIOSAMPLE_HTML_PREFIX + "/" + name + BIOSAMPLE_HTML_SUFFIX), output_file = os.path.join(PROG_DIR, "config", "biosample", (name + ".xml")))
				except Exception as error:
					print("Error: BioSample package " + name + " failed to download.", file=sys.stderr)
					print(error, file=sys.stderr)
				try:
					biosample_package_to_pandera_schema(os.path.join(PROG_DIR, "config", "biosample", (name + ".xml")), name)
				except Exception as error:
					print("Error: BioSample package " + name + " failed to convert to schema.", file=sys.stderr)
					print(error, file=sys.stderr)
				time.sleep(1)
			line = file.readline()
	tools.update_all_schema_templates()

# Convert downloaded BioSample package xml to Pandera Schema
def biosample_package_to_pandera_schema(xml_file: str, name: str) -> None:
	tree = ET.parse(xml_file)
	root = tree.getroot()
	xmlstr = ET.tostring(root, encoding='utf-8', method='xml')
	# Convert xml to dictionary
	report_dict = xmltodict.parse(xmlstr)
	indentation = "\n\t\t"
	mandatory_group: Dict[str, Any] = dict()
	with open(os.path.join(PROG_DIR, "config", "biosample", (name.replace(".", "_") + ".py")), "w+") as file:
		file.writelines(SCHEMA_HEADER)
		for attribute in report_dict["BioSamplePackages"]["Package"]["Attribute"]:
			# If attribute in reserved words skip
			if attribute["HarmonizedName"] in ["collection_date", "gender_restroom"]:
				continue
			# NCBI canonical field name for submission
			file.write(indentation + "\"bs-" + attribute["HarmonizedName"] + "\": Column(")
			# Pandas datatypes
			indentation += "\t"
			file.write(indentation + "dtype=\"object\",")
			# Validation requirements
			if "Format" in attribute:
				if "@type" in attribute["Format"] and attribute["Format"]["@type"] == "select":
					# For columns with only certain valid values
					valid_values = attribute["Format"]["Description"].strip().split(" | ")
					file.write(indentation + r"checks=Check.str_matches(r\"(?i)(\W|^)(" + ("|".join(valid_values)) + r")(\W|$)\"),")
				else:
					file.write(indentation + "checks=None,")
			else:
				file.write(indentation + "checks=None,")
			# Null fields allowed
			if attribute["@use"] == "mandatory":
				file.write(indentation + "nullable=False,")
			else:
				file.write(indentation + "nullable=True,")
			# Every field must be unique
			file.write(indentation + "unique=False,")
			# Coerce column into specified dtype
			file.write(indentation + "coerce=False,")
			# Column is required for submission
			group_message = ""
			if attribute["@use"] == "mandatory":
				file.write(indentation + "required=True,")
			elif attribute["@use"] == "either_one_mandatory":
				# Collect columns that are required but have different column options
				if attribute["@group_name"] in mandatory_group:
					mandatory_group[attribute["@group_name"]] = mandatory_group[attribute["@group_name"]] + " & df[\"bs-" + attribute["HarmonizedName"] + "\"].isnull()"
				else:
					mandatory_group[attribute["@group_name"]] = "df[\"bs-" + attribute["HarmonizedName"] + "\"].isnull()"
				group_message = "At least one required: Group \\\"" + attribute["@group_name"] + "\\\". "
				file.write(indentation + "required=True,")
			else:
				file.write(indentation + "required=False,")
			# NCBI column description
			if attribute["Description"]:
				# Sanitize description
				attribute_description = attribute["Description"].strip().replace("\"", "\\\"").replace("\n"," ")
				if "gender or physical sex" in attribute_description.lower():
					attribute_description = attribute_description.replace("Gender or physical sex", "Biological sex")
				file.write(indentation + "description=\"" + group_message + attribute_description + "\",")
			# Human readable field name for submission
			file.write(indentation + "title=\"" + attribute["Name"] + "\",")
			# Close attribute
			indentation = indentation[:-1]
			file.write(indentation + "),")
		# Close columns
		file.write(SCHEMA_COLUMNS_FOOTER)
		indentation = indentation[:-1]
		file.write(indentation + "},")
		# Create dataframe wide checks
		if bool(mandatory_group):
			file.write(indentation + "checks=[")
			# Validate columns that are required but have multiple options
			indentation += "\t"
			for key in mandatory_group:
				file.write(indentation + "Check(lambda df: ~(" + mandatory_group[key] + "), ignore_na = False),")
			# Close checks
			indentation = indentation[:-1]
			file.write(indentation + "],")
		else:
			file.write(indentation + "checks=None,")
		file.write(indentation + "index=None,")
		file.write(indentation + "coerce=False,")
		file.write(indentation + "strict=\"filter\",")
		file.write(indentation + "name=\"biosample_package_" + name + "_schema\",")
		file.write(indentation + "ordered=False,")
		file.write(indentation + "unique=None,")
		file.write(indentation + "report_duplicates=\"all\",")
		file.write(indentation + "unique_column_names=True,")
		file.write(indentation + "add_missing_columns=False,")
		file.write(indentation + "title=\"BioSample package " + name + " schema\",")
		file.write(indentation + "description=\"Schema validation for BioSample database using " + name + " package.\",")
		# Close schema
		indentation = indentation[:-1]
		file.write(indentation + ")")
	os.remove(xml_file)

def test_internet_connection(databases: List[str]) -> None:
	error_list = []
	print("Checking network settings...", file=sys.stdout)
	for test, info in TEST_CONNECTIONS.items():
		if info["database"] == "GENERAL" or info["database"] in databases:
			print(f"Checking {test} connection...", file=sys.stdout)
			try:
				query = requests.get(info['website'])
				response = query.status_code
			except Exception as e:
				error_list.append(f"{test} connectivity test failed for '{info['website']}'. Check possible firewall issues. \nException:{e}")
			if response in (200, 204, 301, 302):
				print(f"{test} '{info['website']}' connectivity test ok.", file=sys.stdout)
			else:
				error_list.append(f"{info['error_msg']} Error code received:'{response}'")
	if "NCBI" in databases:
		print("Checking DNS resolution for FTP site...", file=sys.stdout)
		try:
			ip_address = socket.gethostbyname(NCBI_FTP_HOST)
		except Exception as e:
			error_list.append(f"Unable to reach '{NCBI_FTP_HOST}'; possible DNS error. \nException:{e}")
		if not ip_address:
			error_list.append(f"Unable to resolve address for '{NCBI_FTP_HOST}'; check DNS server settings for possible issues.")
		else:
			print(f"DNS resolution test ok. Able to reach ('{NCBI_FTP_HOST} -> {ip_address})", file=sys.stdout)
			print("Checking port status...", file=sys.stdout)
			try:
				ftp = ftplib.FTP()
				ftp.connect(NCBI_FTP_HOST, 21, timeout=10)
				ftp.quit()
				print(f"{NCBI_FTP_HOST} open on port 21.", file=sys.stdout)
			except Exception as e:
				error_list.append(f"Port 21 not open for {NCBI_FTP_HOST}. Check possible firewall/server issues. \nException:{e}")
		try:
			print("Checking Table2asn functionality...", file=sys.stdout)
			table2asn_dir = "/tmp/table2asn"
			download_table2asn(table2asn_dir)
			command = [table2asn_dir, "-version-full-xml"]
			print("Running Table2asn.", file=sys.stdout)
			proc = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd = os.path.join(os.path.dirname(os.path.abspath(__file__))))
			if proc.returncode != 0:
				error_list.append("Table2asn-Error")
				error_list.append(proc.stdout)
				error_list.append(proc.stderr)
		except Exception as e:
			error_list.append(f"Unable to download latest version of Table2asn. \nException:{e}")
	if error_list:
		for error_string in error_list:
			print(error_string, file=sys.stderr)
	else:
		print("No network connection issues detected.", file=sys.stdout)
