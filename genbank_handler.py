#!/usr/bin/env python3

###########################    Description    ##################################
# Functions for NCBI input/output handling
################################################################################

# Python Libraries
import pathlib
import pandas as pd
import os
import subprocess
import sys
from lxml import etree
from datetime import datetime
import shutil
import time
from zipfile import ZipFile
from Bio import SeqIO
from distutils.util import strtobool
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from nameparser import HumanName
from typing import List, Set, Dict, Tuple, Optional, Any

# Local imports
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import setup
import file_handler

API_URL = "https://submit.ncbi.nlm.nih.gov/api/2.0/files/FILE_ID/?format=attachment"

# Main create function for BioSample/SRA
def create_genbank_submission(organism: str, submission_name: str, submission_dir: str, config_dict: Dict[str, Any], metadata: pd.DataFrame, gff_file: Optional[str]):
	create_files(organism=organism, submission_name=submission_name, submission_files_dir=submission_files_dir, config_dict=config_dict, metadata=metadata)
	# If using Table2asn do not generate extra genbank files
	if organism not in ["FLU", "COV"]:
		create_table2asn(submission_name=submission_name, submission_files_dir=submission_files_dir, gff_file=gff_file)
		print("Files are stored at: " + os.path.join(submission_files_dir), file=sys.stdout)
		return
	else:
		# If FTP upload for Genbank, create ZIP file for upload if table2asn is set to False
		create_zip(submission_name=submission_name, submission_files_dir=submission_files_dir)
	# Generate NCBI database submission xml
	xml_str = create_submission_xml(organism=organism, submission_name=submission_name, metadata=metadata, config_dict=config_dict, failed_seqs_auto_removed=True)
	file_handler.save_xml(submission_xml=xml_str, submission_files_dir=submission_files_dir)

# Create GenBank XML
def create_submission_xml(organism: str, submission_name: str, config_dict: Dict[str, Any], metadata: pd.DataFrame, failed_seqs_auto_removed: bool = True) -> bytes:
	# Submission XML header
	root = etree.Element("Submission")
	description = etree.SubElement(root, "Description")
	title = etree.SubElement(description, "Title")
	title.text = config_dict["Description"]["Title"]
	comment = etree.SubElement(description, "Comment")
	comment.text = config_dict["Description"]["Comment"]
	# Description info including organization and contact info
	organization = etree.SubElement(description, "Organization", type=config_dict["Description"]["Organization"]["@type"], role=config_dict["Description"]["Organization"]["@role"])
	org_name = etree.SubElement(organization, "Name")
	org_name.text = config_dict["Description"]["Organization"]["Name"]
	action = etree.SubElement(root, "Action")
	addfiles = etree.SubElement(action, "AddFiles", target_db="GenBank")
	file = etree.SubElement(addfiles, "File", file_path=submission_name + ".zip")
	datatype = etree.SubElement(file, "DataType")
	datatype.text = "genbank-submission-package"
	wizard = etree.SubElement(addfiles, "Attribute", name="wizard")
	if "FLU" in organism:
		wizard.text = "BankIt_influenza_api"
	elif "COV" in organism:
		wizard.text = "BankIt_SARSCoV2_api"
	if failed_seqs_auto_removed == True:
		auto_remove = etree.SubElement(addfiles, "Attribute", name="auto_remove_failed_seqs")
		auto_remove.text = "yes"
	identifier = etree.SubElement(addfiles, "Identifier")
	spuid = etree.SubElement(identifier, "SPUID")
	spuid.text = submission_name
	if "FLU" in organism:
		spuid.set("spuid_namespace", "ncbi-influenza-genbank")
	elif "COV" in organism:
		spuid.set("spuid_namespace", "ncbi-sarscov2-genbank")
	# Pretty print xml
	xml_str = etree.tostring(root, encoding="utf-8", pretty_print=True, xml_declaration=True)
	return xml_str

# Create a authorset file
def create_authorset(config_dict: Dict[str, Any], metadata: pd.DataFrame, submission_name: str, submission_files_dir: str) -> None:
	submitter_first = config_dict["Description"]["Organization"]["Submitter"]["Name"]["First"]
	submitter_last = config_dict["Description"]["Organization"]["Submitter"]["Name"]["Last"]
	submitter_email = config_dict["Description"]["Organization"]["Submitter"]["@email"]
	alt_submitter_email = config_dict["Description"]["Organization"]["Submitter"]["@alt_email"]
	affil = config_dict["Description"]["Organization"]["Address"]["Affil"]
	div = config_dict["Description"]["Organization"]["Address"]["Div"]
	publication_status = metadata["gb-publication_status"].unique()[0]
	publication_title = metadata["gb-publication_title"].unique()[0]
	street = config_dict["Description"]["Organization"]["Address"]["Street"]
	city = config_dict["Description"]["Organization"]["Address"]["City"]
	sub = config_dict["Description"]["Organization"]["Address"]["Sub"]
	country = config_dict["Description"]["Organization"]["Address"]["Country"]
	email = config_dict["Description"]["Organization"]["Address"]["Email"]
	phone = config_dict["Description"]["Organization"]["Address"]["Phone"]
	zip_code = str(config_dict["Description"]["Organization"]["Address"]["Postal_code"])
	# Create authorset file
	with open(os.path.join(submission_files_dir, "authorset.sbt"), "w+") as f:
		f.write("Submit-block ::= {\n")
		f.write("  contact {\n")
		f.write("    contact {\n")
		f.write("      name name {\n")
		f.write("        last \"" + submitter_last + "\",\n")
		f.write("        first \"" + submitter_first + "\",\n")
		f.write("        middle \"\",\n")
		f.write("        initials \"\",\n")
		f.write("        suffix \"\",\n")
		f.write("        title \"\"\n")
		f.write("      },\n")
		f.write("      affil std {\n")
		f.write("        affil \""+ affil + "\",\n")
		f.write("        div \"" + div + "\",\n")
		f.write("        city \"" + city + "\",\n")
		f.write("        sub \"" + sub + "\",\n")
		f.write("        country \"" + country + "\",\n")
		f.write("        street \"" + street + "\",\n")
		f.write("        email \"" + email + "\",\n")
		f.write("        phone \"" + phone + "\",\n")
		f.write("        postal-code \"" + zip_code + "\"\n")
		f.write("      }\n")
		f.write("    }\n")
		f.write("  },\n")
		f.write("  cit {\n")
		f.write("    authors {\n")
		f.write("      names std {\n")
		authors = [HumanName(x.strip()) for x in metadata["authors"].unique()[0].split(";") if x.strip() != ""]
		total_names = len(authors)
		for index, name in enumerate(authors, start = 1):
			f.write("        {\n")
			f.write("          name name {\n")
			f.write("            last \"" + name.last + "\",\n")
			f.write("            first \"" + name.first + "\"")
			if name.middle != "":
				f.write(",\n            middle \"" + name.middle + "\"")
			if name.suffix != "":
				f.write(",\n            suffix \"" + name.suffix + "\"")
			if name.title != "":
				f.write(",\n            title \"" + name.title + "\"")
			f.write("\n          }\n")
			if index == total_names:
				f.write("        }\n")
			else:
				f.write("        },\n")
		f.write("      },\n")
		f.write("      affil std {\n")
		f.write("        affil \"" + affil + "\",\n")
		f.write("        div \"" + div + "\",\n")
		f.write("        city \"" + city + "\",\n")
		f.write("        sub \"" + sub + "\",\n")
		f.write("        country \"" + country + "\",\n")
		f.write("        street \"" + street + "\",\n")
		f.write("        postal-code \"" + zip_code + "\"\n")
		f.write("      }\n")
		f.write("    }\n")
		f.write("  },\n")
		f.write("  subtype new\n")
		f.write("}\n")
		f.write("Seqdesc ::= pub {\n")
		f.write("  pub {\n")
		f.write("    gen {\n")
		f.write("      cit \"" + publication_status + "\",\n")
		f.write("      authors {\n")
		f.write("        names std {\n")
		authors = [HumanName(x.strip()) for x in metadata["authors"].unique()[0].split(";") if x.strip() != ""]
		for index, name in enumerate(authors, start = 1):
			f.write("          {\n")
			f.write("            name name {\n")
			f.write("              last \"" + name.last + "\",\n")
			f.write("              first \"" + name.first + "\"")
			if name.middle != "":
				f.write(",\n              middle \"" + name.middle + "\"")
			if name.suffix != "":
				f.write(",\n              suffix \"" + name.suffix + "\"")
			if name.title != "":
				f.write(",\n              title \"" + name.title + "\"")
			f.write("\n            }\n")
			if index == total_names:
				f.write("          }\n")
			else:
				f.write("          },\n")
		f.write("        }\n")
		f.write("      },\n")
		f.write("      title \"" + publication_title + "\"\n")
		f.write("    }\n")
		f.write("  }\n")
		f.write("}\n")
		if alt_submitter_email is not None and alt_submitter_email.strip() != "":
			f.write("Seqdesc ::= user {\n")
			f.write("  type str \"Submission\",\n")
			f.write("  data {\n")
			f.write("    {\n")
			f.write("      label str \"AdditionalComment\",\n")
			f.write("      data str \"ALT EMAIL: " + alt_submitter_email + "\"\n")
			f.write("    }\n")
			f.write("  }\n")
			f.write("}\n")
		f.write("Seqdesc ::= user {\n")
		f.write("  type str \"Submission\",\n")
		f.write("  data {\n")
		f.write("    {\n")
		f.write("      label str \"AdditionalComment\",\n")
		f.write("      data str \"Submission Title: " + submission_name + "\"\n")
		f.write("    }\n")
		f.write("  }\n")
		f.write("}\n")


# Create a zip file for genbank submission
def create_files(organism: str, config_dict: Dict[str, Any], metadata: pd.DataFrame, submission_name: str, submission_files_dir: str) -> None:
	# Create authorset file
	create_authorset(config_dict=config_dict, metadata=metadata, submission_name=submission_name, submission_files_dir=submission_files_dir)
	file_handler.create_fasta(organism=organism, database=["GENBANK"], metadata=metadata, submission_files_dir=submission_files_dir)
	# Retrieve the source df"
	source_df = metadata.filter(regex="^gb-seq_id$|^src-|^ncbi-spuid$|^ncbi-bioproject$|^organism$|^collection_date$").copy()
	source_df.columns = source_df.columns.str.replace("src-","").str.strip()
	source_df = source_df.rename(columns = {"gb-seq_id":"Sequence_ID", "collection_date":"Collection_date"})
	# Add BioProject if available
	if "ncbi-bioproject" in source_df:
		source_df = source_df.rename(columns={"ncbi-bioproject": "BioProject"})
	# Make sure Sequence_ID stays in first column
	shift_col = source_df.pop("Sequence_ID")
	source_df.insert(0, "Sequence_ID", shift_col)
	source_df.to_csv(os.path.join(submission_files_dir, "source.src"), index=False, sep="\t")
	# Retrieve Structured Comment df
	comment_df = metadata.filter(regex="^cmt-")
	if not comment_df.empty:
		comment_df = metadata.filter(regex="^gb-seq_id$|^cmt-|^organism$|^collection_date$").copy()
		comment_df.columns = comment_df.columns.str.replace("cmt-", "").str.strip()
		comment_df = comment_df.rename(columns = {"gb-seq_id": "SeqID"})
		columns_no_prefix_suffix = list(filter(lambda x: (x not in ["SeqID", "StructuredCommentPrefix", "StructuredCommentSuffix"])==True, comment_df.columns))
		ordered_columns = ["SeqID", "StructuredCommentPrefix"] + columns_no_prefix_suffix + ["StructuredCommentSuffix"]
		comment_df = comment_df.reindex(columns=ordered_columns)
		comment_df.to_csv(os.path.join(submission_files_dir, "comment.cmt"), index=False, sep="\t")

# Create a zip file for genbank submission
def create_zip(submission_name: str, submission_files_dir: str) -> None:
	with ZipFile(os.path.join(submission_files_dir, submission_name + ".zip"), 'w') as zip:
		zip.write(os.path.join(submission_files_dir, "authorset.sbt"), "authorset.sbt")
		zip.write(os.path.join(submission_files_dir, "sequence.fsa"), "sequence.fsa")
		zip.write(os.path.join(submission_files_dir, "source.src"), "source.src")
		if os.path.isfile(os.path.join(submission_files_dir, "comment.cmt")):
			zip.write(os.path.join(submission_files_dir, "comment.cmt"), "comment.cmt")
	# Waiting for the zip file to write
	while not os.path.isfile(os.path.join(submission_files_dir, submission_name + ".zip")):
		time.sleep(10)

# Run Table2asn to generate sqn file for submission
def create_table2asn(submission_name: str, submission_files_dir: str, gff_file: Optional[str]) -> Tuple[str, str]:
	submission_status = "processing"
	submission_id = "Table2asn"
	# Create a temp file to store the downloaded table2asn
	table2asn_dir = "/tmp/table2asn"
	# Download the table2asn
	if os.path.isfile(table2asn_dir) is False:
		print("Downloading Table2asn.", file=sys.stdout)
		setup.download_table2asn(table2asn_dir=table2asn_dir)
	# Command to generate table2asn submission file
	command = [table2asn_dir, "-V","vb","-a","s","-t", os.path.join(submission_files_dir, "authorset.sbt"), "-i", os.path.join(submission_files_dir, "sequence.fsa"), "-src-file", os.path.join(submission_files_dir, "source.src"), "-o", os.path.join(submission_files_dir, submission_name + ".sqn")]
	if os.path.isfile(os.path.join(submission_files_dir, "comment.cmt")):
		command.append("-w")
		command.append( os.path.join(submission_files_dir, "comment.cmt"))
	if gff_file is not None:
		command.append("-f")
		command.append(gff_file)
	print("Running Table2asn.", file=sys.stdout)
	proc = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd = os.path.join(os.path.dirname(os.path.abspath(__file__))))
	if proc.returncode != 0:
		print("Table2asn-Error", file=sys.stderr)
		print(proc.stdout, file=sys.stdout)
		print(proc.stderr, file=sys.stderr)
		sys.exit(1)
	print("Validating Table2asn submission.", file=sys.stdout)
	validation_file = os.path.join(submission_files_dir, submission_name + ".val")
	submission_status = check_table2asn_submission(validation_file=validation_file)
	return submission_id, submission_status

# Check table2asn validation information
def check_table2asn_submission(validation_file: str) -> str:
	# Check if validation file exists
	if os.path.isfile(validation_file) == False:
		return "table2asn-error"
	# If submission has errors reject
	with open(validation_file, "r") as file:
		for line in file:
			if "error:" in line.lower():
				print("Submission has errors after running Table2asn.", file=sys.stderr)
				print("Resolve issues labeled \"Error:\" in table2asn validation file or use send_table2asn function to submit with errors.", file=sys.stderr)
				print("Validation file: " + validation_file, file=sys.stderr)
				return "table2asn-validation-error"
			else:
				return "validated"
	return "table2asn-error"

def process_genbank_report(report_file: str, submission_status_file: str, submission_files_dir: str) -> Tuple[str, str]:
	# Read in report.xml
	tree = ET.parse(report_file)
	root = tree.getroot()
	xmlstr = ET.tostring(root, encoding='utf8', method='xml')
	# Read in submission status csv
	status_submission_df = pd.read_csv(submission_status_file, header = 0, dtype = str, engine = "python", encoding="utf-8", index_col=False)
	status_submission_df = status_submission_df.fillna("")
	status_submission_df.columns = status_submission_df.columns.str.strip()
	# Convert xml to dictionary
	report_dict = xmltodict.parse(xmlstr)
	# Get submission status
	try:
		# Get submission status and id from report.xml
		submission_status = report_dict["SubmissionStatus"]["@status"]
	except:
		submission_status = "submitted"
	# Get submission id
	try:
		submission_id = report_dict["SubmissionStatus"]["@submission_id"]
	except:
		submission_id = "pending"
	# CHECK SUBMISSION STATUS ON THE REPORT
	try:
		if report_dict["SubmissionStatus"]["Action"]["@status"] == "processed-ok":
			try:
				for item in report_dict["SubmissionStatus"]["Action"]["Response"]:
					if "File" in item:
						filename_dict = item["File"]
						break
				for file in filename_dict:
					file_name = file["@file_path"]
					file_path = file["@file_id"]
					r = requests.get(API_URL.replace("FILE_ID", file_path), allow_redirects=True)
					open(os.path.join(submission_files_dir, file_name), 'wb').write(r.content)
					# Waiting for the file to write
					while not os.path.exists(os.path.join(submission_files_dir, file_name)):
						time.sleep(10)
					if file_name == "AccessionReport.tsv":
						accession_report_df = pd.read_csv(os.path.join(submission_files_dir, file_name), header = 0, sep = "\t", dtype = str, engine = "python", encoding="utf-8", index_col=False)
						for index, row in accession_report_df.iterrows():
							status_submission_df.loc[status_submission_df["ncbi-sequence_name"] == row["Sequence ID"].strip(), "genbank-status"] = "processed-ok"
							status_submission_df.loc[status_submission_df["ncbi-sequence_name"] == row["Sequence ID"].strip(), "genbank-accession"] = row["#Accession"]
							status_submission_df.loc[status_submission_df["ncbi-sequence_name"] == row["Sequence ID"].strip(), "genbank-message"] = row["Release Date"]
						status_submission_df.to_csv(submission_status_file, header = True, index = False)
			except:
				pass
	except:
		pass
	return submission_status, submission_id

# Check if it has BioSample and BioProject accession number (update status report)
def update_genbank_files(database: List[str], organism: str, submission_files_dir: str, submission_status_file: str) -> None:
	# Read in the submission status report
	status_df = pd.read_csv(submission_status_file, header = 0, dtype = str, engine = "python", encoding="utf-8", index_col=False)
	# Read in genbank source file
	if os.path.isfile(os.path.join(submission_files_dir, "source.src")):
		source_df = pd.read_csv(os.path.join(submission_files_dir, "source.src"), header = 0, sep="\t", dtype = str, engine = "python", encoding="utf-8", index_col=False)
	else:
		print("Error: submission source file does not exist at "+os.path.join(submission_files_dir, "source.src"), file=sys.stderr)
		sys.exit(1)
	# Read in genbank comment file
	if os.path.isfile(os.path.join(submission_files_dir, "comment.cmt")):
		cmt_df = pd.read_csv(os.path.join(submission_files_dir, "comment.cmt"), header = 0, dtype = str, engine = "python", sep="\t", encoding = "utf-8", index_col = False)
	# Retrieve accession info
	src_accessions = dict()
	cmt_accessions = dict()
	# Pull accessions only if field has valid info
	if ("BIOSAMPLE" in database) and ("biosample_accession" in status_df) and (status_df["biosample_accession"].isna().all() == False):
		src_accessions["biosample_accession"] = "BioSample"
	if ("SRA" in database) and ("sra_accession" in status_df) and (status_df["sra_accession"].isna().all() == False):
		src_accessions["sra_accession"] = "SRA"
	if ("GISAID" in database) and ("gisaid_accession_epi_isl_id" in status_df) and (status_df["gisaid_accession_epi_isl_id"].isna().all() == False):
		cmt_accessions["gisaid_accession_epi_isl_id"] = "EPI_ISOLATE_ID"
	if ("GISAID" in database) and ("gisaid_accession_epi_id" in status_df) and (status_df["gisaid_accession_epi_id"].isna().all() == False):
		cmt_accessions["gisaid_accession_epi_id"] = "EPI_SEQUENCE_ID"
	# Update source_df
	if len(src_accessions) > 0:
		# If accession columns exist drop to overwrite
		source_df = source_df.drop(columns=src_accessions.values(), errors="ignore")
		src_accessions["ncbi-sequence_name"] = "Sequence_ID"
		# Merge and rewrite source file
		src_accessions_df = status_df[src_accessions.keys()].copy()
		src_accessions_df = src_accessions_df.rename(columns=src_accessions)
		source_df = pd.merge(source_df, src_accessions_df, how="left", on="Sequence_ID")
		source_df.to_csv(os.path.join(submission_files_dir, "source.src"), index=False, sep="\t")
	# Update CMT file
	if len(cmt_accessions) > 0:
		if os.path.isfile(os.path.join(submission_files_dir, "comment.cmt")):
			cmt_df = pd.read_csv(os.path.join(submission_files_dir, "comment.cmt"), header = 0, dtype = str, engine = "python", sep="\t", encoding = "utf-8", index_col = False)
			# If accession columns exist drop to overwrite
			cmt_df = cmt_df.drop(columns=cmt_accessions.values(), errors="ignore")
			cmt_accessions["ncbi-sequence_name"] = "SeqID"
			# merge fields
			cmt_accessions_df = status_df[cmt_accessions.keys()].copy()
			cmt_accessions_df = cmt_accessions_df.rename(columns=cmt_accessions)
			cmt_df = pd.merge(cmt_df, cmt_accessions_df, how="left", on="SeqID")
		else:
			# If cmt field doesn't exist and must to write accessions then create it
			cmt_df = status_df[cmt_accessions.keys()].copy()
			cmt_df = cmt_df.rename(columns=cmt_accessions)
			if "FLU" in organism:
				cmt_df["StructuredCommentPrefix"] = "FluData"
				cmt_df["StructuredCommentSuffix"] = "FluData"
			elif "COV" in organism:
				cmt_df["StructuredCommentPrefix"] = "Assembly-Data"
				cmt_df["StructuredCommentSuffix"] = "Assembly-Data"
		# Correct order of cmt file columns
		cmt_start = ["SeqID", "StructuredCommentPrefix"]
		cmt_end = ["StructuredCommentSuffix"]
		if "EPI_ISOLATE_ID" in cmt_df:
			cmt_start.append("EPI_ISOLATE_ID")
		if "EPI_SEQUENCE_ID" in cmt_df:
			cmt_end.insert(0, "EPI_SEQUENCE_ID")
		columns_no_prefix_suffix = list(filter(lambda x: (x not in ["SeqID", "StructuredCommentPrefix", "StructuredCommentSuffix", "EPI_ISOLATE_ID", "EPI_SEQUENCE_ID"])==True, cmt_df.columns))
		ordered_columns = cmt_start + columns_no_prefix_suffix + cmt_end
		cmt_df = cmt_df.reindex(columns=ordered_columns)
		cmt_df.to_csv(os.path.join(submission_files_dir, "comment.cmt"), index=False, sep="\t")
