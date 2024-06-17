#!/usr/bin/env python3

###########################	Description	##################################
# Functions for NCBI input/output handling
################################################################################

import pandas as pd
from typing import Set, Dict, Any, Tuple, List

# Check raw reads files listed in metadata file
def check_raw_read_files(submission_name: str, submission_dir: str, metadata: pd.DataFrame) -> Set[str]:
	# Check raw reads files if SRA is provided
	raw_reads_path_default = os.path.join(submission_dir, submission_name, "raw_reads")
	# Separate samples stored in local and cloud
	local_df = metadata[metadata["sra-file_location"] == "local"]
	validated_files = set()
	invalid_raw_files = []
	for index, row in local_df.iterrows():
		# If multiple files check each one
		for file in row["sra-file_name"].split(","):
			file = file.strip()
			# If full path don't use rawreads folder
			if os.path.isabs(file):
				file_path = file
			else:
				# If just file name use rawreads path
				file_path = os.path.join(raw_reads_path, file)
			if os.path.isfile(file_path) == False:
				invalid_raw_files.append(f"Error: Raw read files for {row['sra-sample_name']} does not exist at: {file_path}")
			else:
				validated_files.add(file_path)
	if invalid_raw_files:
		for error_msg in invalid_raw_files:
			print(error_msg, file=sys.stderr)
		sys.exit(1)
	return validated_files

# Create files for optional manual submission to repositories biosample and sra
def create_manual_submission_files(database: str, submission_files_dir: str, metadata: pd.DataFrame, config_dict: Dict[str, Any]) -> None:
	if "SRA" in database:
		metadata_regex = "^ncbi-spuid|^sra-|^organism$|^collection_date$"
		rename_columns = {"sra-library_name":"sra-library_ID"}
		drop_columns = ["ncbi-spuid", "ncbi-spuid_namespace", "sra-file_location", "sra-loader", "sra-file_name"]
		column_ordered = ["sample_name","library_ID"]
		prefix = "sra-"
		# Create SRA specific fields
		metadata["sra-title"] = config_dict["Description"]["Title"]
		metadata = metadata["sra-file_name"].str.split(",", expand=True).add_prefix("sra-filename").join(metadata)
		filename_cols = [col for col in metadata.columns.tolist() if col.startswith("sra-filename")]
		# Correct index for filename column
		for col in filename_cols:
			# Remove 0 index
			if col == "sra-filename0":
				rename_columns[col] = "sra-filename"
			else:
				index = str(int(col.split("filename")[-1]) + 1)
				rename_columns[col] = "sra-filename" + index
	elif "BIOSAMPLE" in database:
		metadata_regex = "^ncbi-|^bs-|^organism$|^collection_date$"
		rename_columns = {"bs-description":"sample_title","ncbi-bioproject":"bioproject_accession"}
		drop_columns = ["ncbi-spuid", "ncbi-spuid_namespace", "bs-package"]
		column_ordered = ["sample_name"]
		prefix = "bs-"
	else:
		print("Error: create_manual_submission_files function only for databases SRA/BioSample. Not " + database + ".", file=sys.stderr)
		sys.exit(1)
	# Filter only needed columns
	database_df = metadata.filter(regex=metadata_regex).copy()
	database_df = database_df.drop_duplicates()
	# Rename columns for manual submission
	database_df = database_df.rename(columns={key: value for key, value in rename_columns.items() if key in database_df})
	# Remove columns not required for manual submission
	database_df = database_df.drop(columns=[col for col in drop_columns if col in database_df.columns])
	# Remove database prefix
	database_df.columns = database_df.columns.str.replace(prefix,"")
	# Order columns
	columns_no_order = list(filter(lambda x: (x not in column_ordered)==True, database_df.columns))
	column_ordered = column_ordered + columns_no_order
	print(database_df.columns)
	database_df = database_df.reindex(columns=column_ordered)
	database_df.to_csv(os.path.join(submission_files_dir, "metadata.tsv"), index=False, sep="\t")

# Create submission XML
def create_submission_xml(organism: str, database: str, submission_name: str, config_dict: Dict[str, Any], metadata: pd.DataFrame, failed_seqs_auto_removed: bool = True) -> bytes:
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
	contact = etree.SubElement(organization, "Contact", email=config_dict["Description"]["Organization"]["Submitter"]["@email"])
	name = etree.SubElement(contact, "Name")
	first_name = etree.SubElement(name, "First")
	first_name.text = config_dict["Description"]["Organization"]["Submitter"]["Name"]["First"]
	last_name = etree.SubElement(name, "Last")
	last_name.text = config_dict["Description"]["Organization"]["Submitter"]["Name"]["Last"]
	# XML actions
	if "BIOSAMPLE" in database:
		database_df = metadata.filter(regex="^ncbi-|^bs-|^organism$|^collection_date$").copy()
		database_df = database_df.drop_duplicates()
		for index, row in database_df.iterrows():
			action = etree.SubElement(root, "Action")
			add_data = etree.SubElement(action, "AddData", target_db="BioSample")
			data = etree.SubElement(add_data, "Data", content_type="xml")
			xmlcontent = etree.SubElement(data, "XmlContent")
			biosample = etree.SubElement(xmlcontent, "BioSample", schema_version="2.0")
			sampleid = etree.SubElement(biosample, "SampleId")
			spuid = etree.SubElement(sampleid, "SPUID", spuid_namespace=row["ncbi-spuid_namespace"])
			spuid.text = row["ncbi-spuid"]
			descriptor = etree.SubElement(biosample, "Descriptor")
			title = etree.SubElement(descriptor, "Title")
			title.text = row["bs-description"]
			organism = etree.SubElement(biosample, "Organism")
			organismname = etree.SubElement(organism, "OrganismName")
			organismname.text = row["organism"]
			if pd.notnull(row["ncbi-bioproject"]) and row["ncbi-bioproject"].strip() != "":
				bioproject = etree.SubElement(biosample, "BioProject")
				primaryid = etree.SubElement(bioproject, "PrimaryId", db="BioProject")
				primaryid.text = row["ncbi-bioproject"]
			package = etree.SubElement(biosample, "Package")
			package.text = row["bs-package"]
			# Attributes
			attributes = etree.SubElement(biosample, "Attributes")
			# Remove columns with bs-prefix that are not attributes
			biosample_cols = [col for col in database_df.columns.tolist() if (col.startswith('bs-')) and (col not in ["bs-package","bs-description"])]
			for col in biosample_cols:
				attribute = etree.SubElement(attributes, "Attribute", attribute_name=col.replace("bs-",""))
				attribute.text = row[col]
			# Add collection date to Attributes
			attribute = etree.SubElement(attributes, "Attribute", attribute_name="collection_date")
			attribute.text = row["collection_date"]
			identifier = etree.SubElement(add_data, "Identifier")
			spuid = etree.SubElement(identifier, "SPUID", spuid_namespace=row["ncbi-spuid_namespace"] + "_bs")
			spuid.text = row["ncbi-spuid"]
	if "SRA" in database:
		database_df = metadata.filter(regex="^ncbi-|^sra-|^organism$|^collection_date$").copy()
		database_df = database_df.drop_duplicates()
		for index, row in database_df.iterrows():
			action = etree.SubElement(root, "Action")
			addfiles = etree.SubElement(action, "AddFiles", target_db="SRA")
			for sra_file in row["sra-file_name"].split(","):
				if row["sra-file_location"].strip().lower() == "cloud":
					file = etree.SubElement(addfiles, "File", cloud_url = sra_file.strip())
				elif row["sra-file_location"].strip().lower() == "local":
					file = etree.SubElement(addfiles, "File", file_path = os.path.basename(sra_file.strip()))
				else:
					print("Error: Metadata field file_location must be either cloud or local. Field currently contains: " + row["sra-file_location"].strip().lower(), file=sys.stderr)
					sys.exit(1)
				datatype = etree.SubElement(file, "DataType")
				datatype.text = "generic-data"
			# Remove columns with sra- prefix that are not attributes
			sra_cols = [col for col in database_df.columns.tolist() if (col.startswith('sra-')) and (col not in ["sra-file_location","sra-file_name"])]
			for col in sra_cols:
				attribute = etree.SubElement(addfiles, "Attribute", name=col.replace("sra-",""))
				attribute.text = row[col]
			if pd.notnull(row["ncbi-bioproject"]) and row["ncbi-bioproject"].strip() != "":
				attribute_ref_id = etree.SubElement(addfiles, "AttributeRefId", name="BioProject")
				refid = etree.SubElement(attribute_ref_id, "RefId")
				primaryid = etree.SubElement(refid, "PrimaryId")
				primaryid.text = row["ncbi-bioproject"]
			if metadata.columns.str.contains("bs-").any():
				attribute_ref_id = etree.SubElement(addfiles, "AttributeRefId", name="BioSample")
				refid = etree.SubElement(attribute_ref_id, "RefId")
				spuid = etree.SubElement(refid, "SPUID", spuid_namespace=row["ncbi-spuid_namespace"] + "_bs")
				spuid.text = row["ncbi-spuid"]
			identifier = etree.SubElement(addfiles, "Identifier")
			spuid = etree.SubElement(identifier, "SPUID", spuid_namespace=row["ncbi-spuid_namespace"] + "_sra")
			spuid.text = row["ncbi-spuid"]
	# Pretty print xml
	xml_str = etree.tostring(root, encoding="utf-8", pretty_print=True, xml_declaration=True)
	return xml_str

# Create list of raw read paths inside sra submission folder
def create_raw_reads_list(submission_files_dir: str, raw_files_list: List[str]) -> None:
	with open(os.path.join(submission_files_dir, "raw_reads_location.txt"), "w+") as file:
		for line in raw_files_list:
			file.write(line + "\n")

# Main create function for BioSample/SRA
def create_biosample_sra_submission(organism: str, database: str, submission_name: str, submission_dir: str, config_dict: Dict[str, Any], metadata: pd.DataFrame):
	if database == "SRA":
		# Validate and write raw reads location
		raw_files_list = check_raw_read_files(submission_name=submission_name, submission_dir=submission_dir, metadata=metadata)
		create_raw_reads_list(raw_files_list)
	create_manual_submission_files(database, submission_files_dir, metadata, config_dict)
	xml_str = create_submission_xml(organism=organism, submission_name=submission_name, metadata=metadata, config_dict=config_dict, failed_seqs_auto_removed=True)
	file_handler.save_xml(xml_str, submission_files_dir)

# Read xml report and get status of the submission
def process_biosample_sra_report(report_file: str, submission_status_file: str) -> Tuple[str, str]:
	# # Read in report.xml
	tree = ET.parse(report_file)
	root = tree.getroot()
	xmlstr = ET.tostring(root, encoding='utf8', method='xml')
	# Read in submission status csv
	status_submission_df = pd.read_csv(submission_status_file, header = 0, dtype = str, engine = "python", encoding="utf-8", index_col=False)
	status_submission_df = status_submission_df.fillna("")
	status_submission_df.columns = status_submission_df.columns.str.strip()
	# Create sample ids that match ids in report.xml
	status_submission_df['report_sample_id'] = [re.sub(r"/", "_", x.lower()) for x in status_submission_df["ncbi-sample_name"].values.tolist()]
	# Convert xml to dictionary
	report_dict = xmltodict.parse(xmlstr)
	# Get submission status from report.xml
	try:
		submission_status = report_dict["SubmissionStatus"]["@status"]
	except:
		submission_status = "submitted"
	# Get submission id from report.xml
	try:
		submission_id = report_dict["SubmissionStatus"]["@submission_id"]
	except:
		submission_id = "pending"
	# Get status of individual samples from report.xml
	try:
		report_action_dict = report_dict["SubmissionStatus"]["Action"]
		# CHECK SUBMISSION ACTIONS ON THE REPORT
		try:
			for i in range(len(report_action_dict)):
				rp_spuid = re.sub(submission_id+"-", "", report_action_dict[i]["@action_id"])
				rp_spuid = re.sub(r"-sra", "", rp_spuid)
				rp_spuid = re.sub(r"-bs", "", rp_spuid)
				rp_target_db = report_action_dict[i]["@target_db"]
				# CHECK THE RESPONSE SECTION OF THE REPORT
				try:
					if (rp_target_db.upper() == "BIOSAMPLE") and (rp_spuid in status_submission_df["report_sample_id"].values.tolist()):
						df_partial = status_submission_df.loc[status_submission_df["report_sample_id"] == rp_spuid]
						status_submission_df.loc[df_partial.index.values, 'biosample_status'] = report_action_dict[i]["@status"]
						try:
							status_submission_df.loc[df_partial.index.values, 'biosample_message'] = report_action_dict[i]["Response"][0]["Message"]["#text"]
						except:
							pass
						if report_action_dict[i]["@status"] == "processed-ok":
							status_submission_df.loc[df_partial.index.values, 'biosample_accession'] = report_action_dict[i]["Response"][0]["Object"]["@accession"]
							status_submission_df.loc[df_partial.index.values , 'bioproject'] = report_action_dict[i]["Response"][0]["Object"]["Meta"]["BioProject"]["#text"]
					elif (rp_target_db.upper() == "SRA") and (rp_spuid in status_submission_df["report_sample_id"].values.tolist()):
						df_partial = status_submission_df.loc[status_submission_df["report_sample_id"] == rp_spuid]
						status_submission_df.loc[df_partial.index.values, 'sra_status'] = report_action_dict[i]["@status"]
						try:
							status_submission_df.loc[df_partial.index.values, 'sra_message'] = report_action_dict[i]["Response"][0]["Message"]["#text"]
						except:
							pass
						if report_action_dict[i]["@status"] == "processed-ok":
							status_submission_df.loc[df_partial.index.values , 'sra_accession'] = report_action_dict[i]["Response"][0]["Object"]["@accession"]
							status_submission_df.loc[df_partial.index.values , 'bioproject'] = report_action_dict[i]["Response"][0]["Object"]["Meta"]["BioProject"]["#text"]
				# CHECK THE RESPONSE SECTION OF THE REPORT
				except:
					try:
						if (rp_target_db.upper() == "BIOSAMPLE") and (rp_spuid in status_submission_df["report_sample_id"].values.tolist()):
							df_partial = status_submission_df.loc[status_submission_df["report_sample_id"] == rp_spuid]
							status_submission_df.loc[df_partial.index.values, 'biosample_status'] = report_action_dict[i]["@status"]
							status_submission_df.loc[df_partial.index.values , 'bioproject'] = report_action_dict[i]["Response"]["Object"]["Meta"]["BioProject"]["#text"]
							try:
								status_submission_df.loc[df_partial.index.values, 'biosample_message'] = report_action_dict[i]["Response"]["Message"]["#text"]
							except:
								pass
							if report_action_dict[i]["@status"] == "processed-ok":
								status_submission_df.loc[df_partial.index.values, 'biosample_accession'] = report_action_dict[i]["Response"]["Object"]["@accession"]
						elif (rp_target_db.upper() == "SRA") and (rp_spuid in status_submission_df["report_sample_id"].values.tolist()):
							df_partial = status_submission_df.loc[status_submission_df["report_sample_id"] == rp_spuid]
							status_submission_df.loc[df_partial.index.values, 'sra_status'] = report_action_dict[i]["@status"]
							try:
								status_submission_df.loc[df_partial.index.values, 'sra_message'] = report_action_dict[i]["Response"]["Message"]["#text"]
							except:
								pass
							if report_action_dict[i]["@status"] == "processed-ok":
								status_submission_df.loc[df_partial.index.values , 'sra_accession'] = report_action_dict[i]["Response"]["Object"]["@accession"]
								status_submission_df.loc[df_partial.index.values , 'bioproject'] = report_action_dict[i]["Response"]["Object"]["Meta"]["BioProject"]["#text"]
					except:
						pass
		except:
			pass
	except:
		pass
	# Save a copy of submission status df without report_sample_id
	status_submission_df = status_submission_df.drop(["report_sample_id"], errors='ignore')
	status_submission_df.to_csv(submission_status_file, header = True, index = False)
	return submission_status, submission_id
