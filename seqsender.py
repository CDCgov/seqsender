#!/usr/bin/env python3

import ftplib
import yaml
import os
import sys
import subprocess
import xmltodict
import pandas as pd
from datetime import datetime
import re
import shutil
from zipfile import ZipFile
import argparse
import time
import json

# Get working and program directory
work_dir = os.getcwd()
prog_dir = os.path.dirname(os.path.abspath(__file__))

# Define version of seqsender
version = "1.0 (Beta)"

# Define current time
STARTTIME = datetime.now()

# Get execution time
def get_execution_time():
    print(f"Total runtime (HRS:MIN:SECS): {str(datetime.now() - STARTTIME)}\n", file=sys.stdout)

# Define NCBI database
ncbi_database = ["BioSample","SRA","BioSample_SRA","Genbank"]

# Define biosample package 
biosample_packages = {"FLU": "Pathogen.cl.1.0", 
					  "COV": "SARS-CoV-2.cl.1.0"}

# Defing genbank wizard
genbank_wizard = {"FLU": "BankIt_influenza_api", 
				  "COV": "BankIt_SARSCoV2_api"}

# Read in the gisaid flu metadata fields
with open(os.path.join(prog_dir, "config", "gisaid_flu_config.yaml"), "r") as f:
	gisaid_flu = yaml.safe_load(f)
	f.close()

# Read in the gisaid covid metadata fields
with open(os.path.join(prog_dir, "config", "gisaid_cov_config.yaml"), "r") as f:
	gisaid_cov = yaml.safe_load(f)
	f.close()				  

# Define required column names for metadata file
sra_required_colnames = ['spuid', 'spuid_namespace', 'bioproject', 'file_location', 'file_path', 'datatype']
sra_required_attributes = ['sra-instrument_model', 'sra-library_name', 'sra-library_strategy', 'sra-library_source', 'sra-library_selection', 'sra-library_layout', 'sra-library_construction_protocol']
biosample_required_colnames = ["spuid", "spuid_namespace", "bioproject", "description_title", "organism_name"]
bioSample_required_attributes = {"FLU": ["bs-collected_by", "bs-collection_date", "bs-geo_loc_name", "bs-host", "bs-host_disease", "bs-strain", "bs-isolation_source", "bs-lat_lon"],
								 "COV": ["bs-collected_by", "bs-collection_date", "bs-geo_loc_name", "bs-host", "bs-host_disease", "bs-isolate", "bs-isolation_source"]}
genbank_required_colnames = {"FLU": ["Sequence_ID", "src-strain", "src-collection_date", "src-country", "src-host", "src-serotype", "src-isolation_source"],
						     "COV": ["Sequence_ID", "src-organism", "src-isolate", "src-collection_date", "src-country", "src-host", "src-isolation_source"]}
all_required_colnames = {"FLU": {"BioSample": biosample_required_colnames + bioSample_required_attributes["FLU"], 
						 		 "SRA": sra_required_colnames + sra_required_attributes,
						 		 "BioSample_SRA": [*set(biosample_required_colnames + bioSample_required_attributes["FLU"] + sra_required_colnames + sra_required_attributes)],
						 		 "Genbank": genbank_required_colnames["FLU"],
						 		 "GISAID": list(filter(lambda x: (("**" in x) or ("***" in x))==True, gisaid_flu["EpiFlu"].keys()))
					 			},
			             "COV": {"BioSample": biosample_required_colnames + bioSample_required_attributes["COV"], 
						 		 "SRA": sra_required_colnames + sra_required_attributes,
						 		 "BioSample_SRA": [*set(biosample_required_colnames + bioSample_required_attributes["COV"] + sra_required_colnames + sra_required_attributes)],
						 		 "Genbank": genbank_required_colnames["COV"],
						 		 "GISAID": list(filter(lambda x: (("**" in x) or ("***" in x))==True, gisaid_cov["EpiCoV"].keys()))
					 			}
						}

print(all_required_colnames["FLU"]["GISAID"]); print(all_required_colnames["COV"]["GISAID"]);

# Create gisaid-flu Segment Ids
Segment_Ids = ["HA", "HE", "MP", "NA", "NP", "NS", "P3", "PA", "PB1", "PB2"]

# Define required the sample colunm names to check if samples are unique
sample_colnames = {"FLU": {"BioSample": "spuid",
						   "SRA": "spuid",
						   "BioSample_SRA": "spuid",
						   "Genbank": "Sequence_ID",
						   "GISAID": ["Seq_Id ("+x+")" for x in Segment_Ids]
						  },
			 	   "COV": {"BioSample": "spuid",
						   "SRA": "spuid",
						   "BioSample_SRA": "spuid",
						   "Genbank": "Sequence_ID",								 
						   "GISAID": "covv_virus_name"
						  }
				  }

# NCBI FTP Server
NCBI_FTP_HOST = "ftp-private.ncbi.nlm.nih.gov"

# NCBI API URL
NCBI_API_URL = "https://submit.ncbi.nlm.nih.gov/api/2.0/files/FILE_ID/?format=attachment"

# Read in config file
def get_config(config_file):
	if os.path.isfile(config_file) == False:
		print("Error: Config file does not exist at: " + config_file, file=sys.stderr)
		sys.exit(1)
	else:
		with open(config_file, "r") as f:
			config_dict = yaml.safe_load(f)
			f.close()
		if type(config_dict) is dict:
			try:
				config_dict = config_dict["Submission"]
				return config_dict
			except:
				print("Error: there is no Submission information in config file.", file=sys.stderr) 
				sys.exit(1)					
		else:
			print("Error: Config file is incorrect. Config file must be a dictionary.", file=sys.stderr)
			sys.exit(1)

# Read in token file
def get_token(token_file, organism):
	if os.path.isfile(token_file) == False:
		print("Error: Token file does not exist at: " + token_file, file=sys.stderr)
		sys.exit(1)
	else:
		with open(token_file, "r") as f:
			token_dict = yaml.safe_load(f)
			f.close()			
		if type(token_dict) is dict:
			try:
				token_dict = token_dict[organism.lower()]
				return token_dict
			except:
				print("Error: there is no " + organism.lower() + " authenticated information in token file.", file=sys.stderr) 
				sys.exit(1)	
		else:
			print("Error: Token file is incorrect. Token file must be a dictionary.", file=sys.stderr)
			sys.exit(1)		

# Read in metadata file
def get_metadata(organism, database, metadata_file):
	# Check if file exists
	if os.path.isfile(metadata_file) == False:
		print("Error: Metadata file does not exist at: " + metadata_file, file=sys.stderr)
		sys.exit(1)
	else:
		metadata = pd.read_csv(metadata_file, header = 0, dtype = str, engine = "python", encoding="utf-8", index_col=False)
		metadata = metadata.fillna("")
		metadata.columns = metadata.columns.str.strip() 
		# In cov gisaid, fields marked by "**" are required. If a value is missing in those fields, use the term "unknown"
		gisaid_unknown_colnames = []
		if organism == "COV" and database == "GISAID":
			gisaid_unknown_colnames += list(filter(lambda x: (("**" in x) and ("***" not in x))==True, all_required_colnames[organism][database]))		
		# Obtain the required column names without the asterisks 
		required_colnames = [str(x).replace("*", "") for x in  all_required_colnames[organism][database]]
		# Check if required columns are existed in metadata file
		if not set(required_colnames).issubset(set(metadata.columns)):
			failed_rq_colnames = list(filter(lambda x: (x in metadata.columns)==False, required_colnames))
			print("Error: Metadata file must have the following required column names: " + ", ".join(failed_rq_colnames) + ".", file=sys.stderr)
			sys.exit(1)
		# In flu gisaid, Isolate_Id and Segment_Ids must be empty String
		if organism == "FLU" and database == "GISAID":
			# Isolate_Id and Segment_Ids must be empty String
			if any(metadata["Isolate_Id"] != "") and any(metadata["Segment_Ids"] != ""):
				print("Error: 'Isolate_Id' and 'Segment_Ids' in metadata file must be empty String.", file=sys.stderr)
				sys.exit(1)
			# In flu gisaid, submission must has at least one segment name
			segment_names = list(filter(lambda x: all(metadata["Seq_Id ("+x+")"]==""), Segment_Ids))
			if len(segment_names) == len(Segment_Ids):
				print("Error: Metadata file must has at least one segment name: " + ", ".join(["Seq_Id ("+s+")" for s in Segment_Ids]) + " that is not empty.", file=sys.stderr)
				sys.exit(1)	
		# In Genbank, Sequence_ID field cannot contain spaces
		if database == "Genbank" and len(list(filter(lambda x: (" " in str(x))==True, metadata["Sequence_ID"].values.tolist()))) > 0:
			print("Error: 'Sequence_ID' field in metadata file cannot contain spaces.", file=sys.stderr)
			sys.exit(1)
		# Make sure required fields do not contain any empty values
		for name in required_colnames:
			# Make sure required fields in gisaid with empty values are filled with 'Unknown'	
			if (name in gisaid_unknown_colnames) and any(metadata[name]==""):
				print("Error: The required '" + name + "' field in metatdata file cannot contain any empty values. Please fill empty values with 'Unknown'", file=sys.stderr)
				sys.exit(1)
			# Make sure required fields do not contain any empty values
			elif (name not in gisaid_unknown_colnames) and pd.isna(metadata[name]).any():
				print("Error: The required '" + name + "' field in metatdata file cannot contain any empty values.", file=sys.stderr)
				sys.exit(1)
		# Make sure the column contains sample names are unique
		unique_colnames = sample_colnames[organism][database]
		for name in unique_colnames:
			unique_name = set(metadata[name])
			if len(unique_name) != metadata.shape[0]:
				print("Error: The required '" + name + "' field in metatdata file must have samples that are unique.", file=sys.stderr)
				sys.exit(1)
		return metadata

# Check user credentials information
def get_credentials(config_file, database):
	# Get config file
	config_dict = get_config(config_file)		
	# Extract submission information
	try:
		submission_dict = {k.title():v for k,v in config_dict[database].items()}
	except:
		print("Error: there is no Submission > " + database + " information in config file.", file=sys.stderr)
		sys.exit(1)
	else:
		# Create an empty dict to store credentials
		cred_dict = {}
		# Check username		
		if "Username" not in submission_dict.keys():
			print("Error: there is no Submission > " + database + " > Username information in config file.", file=sys.stderr)
			sys.exit(1)	
		elif ("Username" in submission_dict.keys()) and ((submission_dict["Username"] is None) or (submission_dict["Username"] == "")):
			print("Error: Submission > " + database + " > Username in the config file cannot be empty.", file=sys.stderr)
			sys.exit(1)	
		elif ("Username" in submission_dict.keys()) and ((submission_dict["Username"] is not None) and (submission_dict["Username"] != "")):
			cred_dict["Username"] = submission_dict["Username"]		
		# Check password		
		if "Password" not in submission_dict.keys():
			print("Error: there is no Submission > " + database + " > Password information in config file.", file=sys.stderr)
			sys.exit(1)		
		elif ("Password" in submission_dict.keys()) and ((submission_dict["Password"] is None) or (submission_dict["Password"] == "")):
			print("Error: Submission > " + database + " > Password in the config file cannot be empty.", file=sys.stderr)
			sys.exit(1)	
		elif ("Password" in submission_dict.keys()) and ((submission_dict["Password"] is not None) and (submission_dict["Password"] != "")):
			cred_dict["Password"] = submission_dict["Password"]	
		# Check client-id if database is GISAID		
		if (database == "GISAID") and ("Client-Id" not in submission_dict.keys()):
			print("Error: there is no Submission > " + database + " > Client-Id information in config file.", file=sys.stderr)
			sys.exit(1)	
		elif (database == "GISAID") and ("Client-Id" in submission_dict.keys() and ((submission_dict["Client-Id"] is None) or (submission_dict["Client-Id"] == ""))):
			print("Error: Submission > " + database + " > Client-Id in the config file cannot be empty.", file=sys.stderr)
			sys.exit(1)	
		elif (database == "GISAID") and ("Client-Id" in submission_dict.keys() and ((submission_dict["Client-Id"] is not None) and (submission_dict["Client-Id"] != ""))):
			cred_dict["Client-Id"] = submission_dict["Client-Id"]	
		return cred_dict

# Authenticate user credentials
def authenticate(database, organism, config_file):	
	if database == "GISAID":
		credentials = get_credentials(config_file=config_file, database="GISAID")		
		if organism == "FLU":
			command = "python3 /MIRA/epiflu_cli/__main__.py authenticate --token gisaid.flu.authtoken --username %s --password %s --client_id %s --log gisaid-logfileFlu.log --force > /dev/null" % (credentials["Username"], credentials["Password"], credentials["Client-Id"])
		elif organism == "COV":
			command = "python3 /MIRA/epicov_cli/__main__.py authenticate --token gisaid.cov.authtoken --username %s --password %s --client_id %s --log gisaid-logfileCoV.log --force > /dev/null" % (credentials["Username"], credentials["Password"], credentials["Client-Id"])
		# Check status of the command
		if os.system(command) == 0:
			print("\n" + "User authenticated", file=sys.stdout)
		else:
			print("\n" + "Error: " + e, file=sys.stderr)
			sys.exit(1)
	else:
		credentials = get_credentials(config_file=config_file, database="NCBI")		
		# Login into NCBI FTP Server
		try:
			ftp = ftplib.FTP(NCBI_FTP_HOST)
			ftp.login(user=credentials["Username"], passwd = credentials["Password"])
			if organism == "FLU":
				token_file = {"flu": {"username": credentials["Username"], "password": credentials["Password"]}}
			elif organism == "COV":
				token_file = {"cov": {"username": credentials["Username"], "password": credentials["Password"]}}				
			# Save an authentication file 
			with open(os.path.join(work_dir, "ncbi."+organism.lower()+".authtoken"), "w") as f:
				yaml.dump(token_file, f)
				f.close()
			# Output authenticated message
			print("\n" + "User authenticated", file=sys.stdout)
		except ftplib.all_errors as e:
			print("\n" + "Authentication failed", file=sys.stderr)
			print("\n" + 'Error:' + e, file=sys.stderr)
			sys.exit(1)

# Create description section in submission.xml
def create_submission_description(submission_dict, database):
	# Check if config has description information
	if "Description" not in submission_dict.keys():
		print("Error: there is no Submission > " + database + " > Description information in config file.", file=sys.stderr)
		sys.exit(1)	
	# Check if description has a title (required)
	try:
		title = submission_dict["Description"]["Title"]
	except:
		print("Error: there is no Submission > " + database + " > Description > Title in the config file.", file=sys.stderr)
		sys.exit(1)	
	else:
		if (title is None) or (title == ""):
			print("Error: Submission > " + database + " > Description > Title in the config file cannot be empty.", file=sys.stderr)
			sys.exit(1)	
	# Check if description has organization information (required) for NCBI
	if database == "NCBI":
		try:
			organization = submission_dict["Description"]["Organization"]
		except:
			print("Error: there is no Submission > " + database + " > Description > Organization in the config file.", file=sys.stderr)
			sys.exit(1)	
		else:
			if (organization["@role"] is None) or (organization["@role"] == ""):
				print("Error: Submission > " + database + " > Description > Organization > @role in the config file cannot be empty.", file=sys.stderr)
				sys.exit(1)	
			if (organization["@type"] is None) or (organization["@type"] == ""):
				print("Error: Submission > " + database + " > Description > Organization > @type in the config file cannot be empty.", file=sys.stderr)
				sys.exit(1)	
			if (organization["Name"] is None) or (organization["Name"] == ""):
				print("Error: Submission > " + database + " > Description > Organization > Name in the config file cannot be empty.", file=sys.stderr)
				sys.exit(1)
			if (organization["Contact"] is None) or (organization["Contact"] == ""):
				print("Error: Submission > " + database + " > Description > Organization > Contact in the config file cannot be empty.", file=sys.stderr)
				sys.exit(1)
			else:
				contact = organization["Contact"]
				if(contact["@email"] is None) or (contact["@email"] == ""):
					print("Error: Submission > " + database + " > Description > Organization > Contact > @email in the config file cannot be empty.", file=sys.stderr)
					sys.exit(1)			
				if(contact["Name"] is None) or (contact["Name"] == ""):
					print("Error: Submission > " + database + " > Description > Organization > Contact > Name in the config file cannot be empty.", file=sys.stderr)
					sys.exit(1)	
				else:
					contact_name = contact["Name"]
					if(contact_name["First"] is None) or (contact_name["First"] == ""):
						print("Error: Submission > " + database + " > Description > Organization > Contact > Name > First in the config file cannot be empty.", file=sys.stderr)
						sys.exit(1)	
					if(contact_name["Last"] is None) or (contact_name["Last"] == ""):
						print("Error: Submission > " + database + " > Description > Organization > Contact > Name > Last in the config file cannot be empty.", file=sys.stderr)
						sys.exit(1)			
	return submission_dict["Description"]

# Create action section in submission.xml for BioSample
def create_biosample_submission(metadata, package):
	# Create a dict to combine action for each submission
	combined_action_dict = {'AddData': [""]*metadata.shape[0]}
	## Retrieve the attributes df
	attributes_df = metadata.filter(regex="^bs-")
	## Obtain the names of the attributes without the bs-
	attributes_df.columns = [re.sub(r"^bs-", "", str(x)) for x in attributes_df.columns]
	attributes_df.columns = attributes_df.columns.str.strip() 
	attributes_colnames = attributes_df.columns.values.tolist()
	# Create an action for each subission in metadata file
	for index, row in metadata.iterrows():
		# Read in the action config file
		action_config_path = os.path.join(prog_dir, "config", "biosample_action_config.yaml")
		action_config_dict = get_config(config_file=action_config_path)	
		# Fill in the spuid and spuid_namespace information
		action_config_dict["Action"]["AddData"]["Data"]["XmlContent"]["BioSample"]["SampleId"]["SPUID"]["#text"] = row["spuid"]
		action_config_dict["Action"]["AddData"]["Data"]["XmlContent"]["BioSample"]["SampleId"]["SPUID"]["@spuid_namespace"] = row["spuid_namespace"]
		action_config_dict["Action"]["AddData"]["Identifier"]["SPUID"]["#text"] = row["spuid"]
		action_config_dict["Action"]["AddData"]["Identifier"]["SPUID"]["@spuid_namespace"] = row["spuid_namespace"]
		## Fill in the description title and organism name information
		action_config_dict["Action"]["AddData"]["Data"]["XmlContent"]["BioSample"]["Descriptor"]['Title'] = row["description_title"]
		action_config_dict["Action"]["AddData"]["Data"]["XmlContent"]["BioSample"]["Organism"]["OrganismName"] = row["organism_name"]
		## Fill in bioproject information
		action_config_dict["Action"]["AddData"]["Data"]["XmlContent"]["BioSample"]["BioProject"]["PrimaryId"]["@db"] = row["bioproject"]
		## Fill in package information
		action_config_dict["Action"]["AddData"]["Data"]["XmlContent"]["BioSample"]["Package"] = package
		## Fill in attributes information
		combined_attribute_dict = dict()
		combined_attribute_dict.setdefault("Attribute", [""]*len(attributes_colnames))
		for i in range(len(attributes_colnames)):
			attr_dict = {}
			attr_dict["#text"] = attributes_df.loc[index, attributes_colnames[i]]
			attr_dict["@attribute_name"] = attributes_colnames[i]
			combined_attribute_dict["Attribute"][i] = attr_dict
		## Update attributes information
		action_config_dict["Action"]["AddData"]["Data"]["XmlContent"]["BioSample"]["Attributes"] = combined_attribute_dict
		## Update action information
		combined_action_dict['AddData'][index] = action_config_dict["Action"]["AddData"]
	return combined_action_dict

# Check fasta files listed in metadata file
def check_fasta_files(metadata, fasta_path):
	# Check fasta local location is local or stored on cloud
	failed_fasta_location = list(filter(lambda x: (x in ["local", "cloud"])==False, metadata["file_location"]))
	if len(failed_fasta_location) > 0:
		print("Error: the value of file_location in metatdata file can only be local/cloud.", file=sys.stderr)
		sys.exit(1)	
	# Separate samples stored in local and cloud 
	local_fasta_files = []; cloud_fasta_files = [];
	for index, row in metadata.iterrows():
		# Create an empty list to store all fasta files
		if row["file_location"] == "local":
			local_fasta_files += [os.path.join(fasta_path, g.strip()) for g in row["file_path"].split(",")]
		elif row["file_location"] == "cloud":
			cloud_fasta_files += [os.path.join(fasta_path, g.strip()) for g in row["file_path"].split(",")]
	# Get a list of fasta paths that does not exist locally
	if len(local_fasta_files) > 0:
		failed_fasta = list(filter(lambda x: os.path.isfile(x)==False, local_fasta_files))
		if len(failed_fasta) > 0:
			print("Error: fasta files do not exist at " + ", ".join(failed_fasta) + ".", file=sys.stderr)
			sys.exit(1)	
	# Get a list of fasta paths that does not exist on aws cloud
	if len(cloud_fasta_files) > 0:
		failed_fasta = list(filter(lambda x: os.path.isfile(x)==False, local_fasta_files))
		if len(failed_fasta) > 0:
			print("Error: fasta files do not exist at " + ", ".join(failed_fasta) + ".", file=sys.stderr)
			sys.exit(1)	

# Create action section in submission.xml for SRA
def create_sra_submission(metadata, fasta_path):	
	# Check fasta files listed in metadata file exist local or on cloud given fasta path
	check_fasta_files(metadata=metadata, fasta_path=fasta_path)
	# Create a dict to combine action for each submission
	combined_action_dict = {"AddFiles": [""]*metadata.shape[0]}
	## Retrieve the attributes df
	attributes_df = metadata.filter(regex="^sra-")
	## Obtain the names of the attributes without sra-
	attributes_df.columns = [re.sub(r'^sra-', '', str(x)) for x in attributes_df.columns]
	attributes_df.columns = attributes_df.columns.str.strip() 
	attributes_colnames = attributes_df.columns.values.tolist()
	# Create an action for each subission in metadata file
	for index, row in metadata.iterrows():
		# Read in the action config file
		action_config_path = os.path.join(prog_dir, "config", "sra_action_config.yaml")
		action_config_dict = get_config(config_file=action_config_path)	
		# Fill in the spuid and spuid_namespace information
		action_config_dict["Action"]["AddFiles"]["Identifier"]["SPUID"]["#text"] = row["spuid"]
		action_config_dict["Action"]["AddFiles"]["Identifier"]["SPUID"]["@spuid_namespace"] = row["spuid_namespace"]
		## Fill in bioproject information
		action_config_dict["Action"]["AddFiles"]["AttributeRefId"][0]["@name"] = "BioProject"
		action_config_dict["Action"]["AddFiles"]["AttributeRefId"][0]["RefId"]["PrimaryId"]["#text"] = row["bioproject"]
		action_config_dict["Action"]["AddFiles"]["AttributeRefId"][0]["RefId"]["PrimaryId"]["@db"] = "BioProject"
		## Fill in biosample_accession information
		if metadata.columns.isin(["biosample_accession"]).any() and row["biosample_accession"] != "":
			action_config_dict["Action"]["AddFiles"]["AttributeRefId"][1]["@name"] = "BioSample"
			action_config_dict["Action"]["AddFiles"]["AttributeRefId"][1]["RefId"]["PrimaryId"]["#text"] = row["biosample_accession"]
			action_config_dict["Action"]["AddFiles"]["AttributeRefId"][1]["RefId"]["PrimaryId"]["@db"] = "BioSample"
			del action_config_dict["Action"]["AddFiles"]["AttributeRefId"][1]["RefId"]["SPUID"]
		else:
			action_config_dict["Action"]["AddFiles"]["AttributeRefId"][1]["@name"] = "BioSample"
			action_config_dict["Action"]["AddFiles"]["AttributeRefId"][1]["RefId"]["SPUID"]["#text"] = row["spuid"]
			action_config_dict["Action"]["AddFiles"]["AttributeRefId"][1]["RefId"]["SPUID"]["@spuid_namespace"] = row["spuid_namespace"]
			del action_config_dict["Action"]["AddFiles"]["AttributeRefId"][1]["RefId"]["PrimaryId"]
		## Get fasta file information
		fasta = [file.strip() for file in row["file_path"].split(",")]
		## Fill in fasta file information
		combined_fasta_dict = dict()
		combined_fasta_dict.setdefault("File", [""]*len(fasta))
		for i in range(len(fasta)):
			fasta_dict = {}
			if row["file_location"] == "local":
				fasta_dict["@file_path"] = fasta[i]
			elif row["file_location"] == "cloud":
				fasta_dict["@cloud_url"] = fasta[i]
			fasta_dict["DataType"] = row["datatype"]
			combined_fasta_dict["File"][i] = fasta_dict
		action_config_dict["Action"]["AddFiles"]["File"] = combined_fasta_dict["File"]
		## Fill in attributes information
		combined_attribute_dict = dict()
		combined_attribute_dict.setdefault("Attribute", [""]*len(attributes_colnames))
		for i in range(len(attributes_colnames)):
			attr_dict = {}
			attr_dict["#text"] = attributes_df.loc[index, attributes_colnames[i]]
			attr_dict["@name"] = attributes_colnames[i]
			combined_attribute_dict["Attribute"][i] = attr_dict
		## Update attributes information
		action_config_dict["Action"]["AddFiles"]["Attribute"] = combined_attribute_dict["Attribute"]
		## Update action information
		combined_action_dict["AddFiles"][index] = action_config_dict["Action"]["AddFiles"]
	return combined_action_dict

# Check sample names in metadata file are listed in fasta file
def check_fasta_samples(metadata, fasta_file, sample_colnames):
	for index, row in metadata.iterrows():
		for name in sample_colnames:
			command = os.system("grep -c -w >%s %s" % (row[name], fasta_file))
			if os.system(command) == 0:
				print("Error: Cannot find " + name + " = " + row[sample] + " in " + fasta_file + ". Make sure to start a sample with >[sample_name].", file=sys.stderr)
				sys.exit(1)
			elif os.system(command) > 1:
				print("Error: " + name + " = " + row[sample] + " is duplicated in " + fasta_file + ".", file=sys.stderr)
				sys.exit(1)

# Create a zip file for genbank submission
def create_genbank_zip(organism, metadata, fasta_file, authorset_file, submission_name, submission_dir, comment=False):
	# Check each row to make sure each sample name exists in fasta file
	check_fasta_samples(metadata=metadata, fasta_file=fasta_file, sample_colnames=sample_colnames[organism]["Genbank"])
	# Copy fasta to output directory
	shutil.copyfile(fasta_file, os.path.join(submission_dir, "sequence.fsa"))
	# Copy authorset to output directory
	shutil.copyfile(authorset_file, os.path.join(submission_dir, "authorset.sbt"))
	# Retrieve the source df
	source_df = metadata.filter(regex="Sequence_ID|^src-")
	source_df.columns = [re.sub(r'^src-', '', str(x)) for x in source_df.columns]
	source_df.columns = source_df.columns.str.strip() 
	source_df.to_csv(os.path.join(submission_dir, "source.src"), index=False, sep="\t")
	# Retrieve comment df
	comment_df = metadata.filter(regex="^cmt-")
	if comment_df.shape[0] > 0:
		comment_df = metadata.filter(regex="Sequence_ID|^cmt-")
		comment_df.columns = [re.sub(r'^cmt-', '', str(x)) for x in comment_df.columns]
		comment_df.columns = [re.sub(r'Sequence_ID', 'SeqID', str(x)) for x in comment_df.columns]
		#comment_df.columns = [re.sub(r'-', ' ', str(x)) for x in comment_df.columns]
		comment_df.columns = comment_df.columns.str.strip() 
		comment_df.to_csv(os.path.join(submission_dir, "comment.cmt"), index=False, sep="\t")
		comment = True
	# Create a zip file for genbank submission
	with ZipFile(os.path.join(submission_dir, submission_name + ".zip"), 'w') as zip:
		zip.write(os.path.join(submission_dir, "authorset.sbt"), "authorset.sbt")
		zip.write(os.path.join(submission_dir, "sequence.fsa"), "sequence.fsa")
		zip.write(os.path.join(submission_dir, "source.src"), "source.src")
		if comment == True:
			zip.write(os.path.join(submission_dir, "comment.cmt"), "comment.cmt")
	# Waiting for the zip file to write
	while not os.path.exists(os.path.join(submission_dir, submission_name + ".zip")):
		time.sleep(10)

# Create action section in submission.xml for Genbank
def create_genbank_submission(organism, metadata, fasta_file, authorset_file, submission_name, submission_dir):
	# Create a zip file for genbank submission
	create_genbank_zip(metadata=metadata, fasta_file=fasta_file, authorset_file=authorset_file, submission_name=submission_name, submission_dir=submission_dir)
	# Read in the action config file
	action_config_path = os.path.join(prog_dir, "config", "genbank_action_config.yaml")
	action_config_dict = get_config(action_config_path)["Action"]
	# Fill in the spuid and spuid_namespace information
	action_config_dict["AddFiles"]["Identifier"]["SPUID"]["#text"] = submission_name
	action_config_dict["AddFiles"]["Identifier"]["SPUID"]["@spuid_namespace"] = "ncbi-"+organism.lower()+"-genbank"
	# Fill in genbank wizard information
	action_config_dict["AddFiles"]["Attribute"][0]["#text"] = genbank_wizard[organism]
	# Fill zip file information
	action_config_dict["AddFiles"]["File"]["@file_path"] = submission_name + ".zip"
	return action_config_dict

def create_submission_xml(organism, database, config_file, metadata_file, fasta_file=None, fasta_path=None, authorset_file=None, test=False):
	# Get config file
	config_dict = get_config(config_file=config_file)
	# Get metadata file
	metadata = get_metadata(organism=organism, database=database, metadata_file=metadata_file)
	# Check if fasta file is provided
	if fasta_file is not None and os.path.isfile(fasta_file) == False:
		print("Error: fasta file does not exist at: \n" + fasta_file, file=sys.stderr)
		sys.exit(1)	
	# Check if fasta path is provided
	if fasta_path is not None and os.path.isdir(fasta_path) == False:
		print("Error: fasta path does not exist at: \n" + fasta_path, file=sys.stderr)
		sys.exit(1)	
	# Check if authorset file is provided
	if authorset_file is not None and  os.path.isfile(authorset_file) == False:
		print("Error: authorset file does not exist at: " + authorset_file, file=sys.stderr)
		sys.exit(1)	
	# Extract NCBI submission information
	try:
		submission_dict = {k.title():v for k,v in config_dict["NCBI"].items()}
	except:
		print("Error: there is no Submission > NCBI information in config file.", file=sys.stderr)
		sys.exit(1)
	else:
		# Extract submission description information
		description_dict = create_submission_description(submission_dict=submission_dict)
		# Extract description title to store as directory name to dump submission files
		submission_name = description_dict["Title"]
		# Create an output directory
		submission_dir = os.path.join(work_dir, "submission"+str(datetime.now().strftime("%y%m%d")), submission_name, database, organism)
		# Check if output directory exists	
		if os.path.exists(submission_dir) == False:
			os.makedirs(submission_dir)	
		# Extract submission action information
		if database == "BioSample":
			# Create biosample submission dict
			bs_submission_dict = create_biosample_submission(metadata=metadata, package=biosample_packages[organism])
			# Create final submission file
			submission_xml = {"Submission": {"Description": description_dict, "Action1": bs_submission_dict}}				
		elif database == "SRA":
			# Create sra submission dict
			sra_submission_dict = create_sra_submission(metadata=metadata, fasta_path=fasta_path)
			# Create final submission file
			submission_xml = {"Submission": {"Description": description_dict, "Action2": sra_submission_dict}}	
		elif database == "BioSample_SRA":
			# Create biosample submission dict
			bs_submission_dict = create_biosample_submission(metadata=metadata, package=package)
			# Create sra submission dict
			sra_submission_dict = create_sra_submission(metadata=metadata, fasta_path=fasta_path)
			# Combined biosample and sra as one final submission file
			submission_xml = {"Submission": {"Description": description_dict, "Action1": bs_submission_dict, "Action2": sra_submission_dict}}				
		elif database == "Genbank":
			# Create genebank submission dict
			genbank_submission_dict = create_genbank_submission(organism=organism, metadata=metadata, fasta_file=fasta_file, authorset_file=authorset_file, submission_name=submission_name, submission_dir=submission_dir)		
			# Create final submission file
			submission_xml = {"Submission": {"Description": description_dict, "Action": genbank_submission_dict}}							
		# Create submission string
		xmlstr = xmltodict.unparse(submission_xml, pretty=True)
		# Save string as submission.xml
		with open(os.path.join(submission_dir, "submission.xml"), "w") as f:
			f.write(xmlstr)
			f.close()
		# Waiting for the xml file to write
		while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
			time.sleep(10)			
		# Replace Action1, Action2 with just Action in submission.xml
		update_xml_cmd = "bash %s/update_submission_xml.sh %s" % (prog_dir, os.path.join(submission_dir, "submission.xml"))
		# Check status of the bash command
		if os.system(update_xml_cmd) == 0:
			print("\n" + "submission.xml is created", file=sys.stdout)
			print("\n" + "File is stored at: " + os.path.join(submission_dir, "submission.xml"), file=sys.stdout)
			return {"submission_name": submission_name, "submission_dir": submission_dir}
		else:
			print("\n" + "Creating submission.xml failed", file=sys.stderr)
			sys.exit(1)			

# Submit to NCBI
def submit_ncbi(organism, database, config_file, metadata_file, fasta_file=None, fasta_path=None, authorset_file=None, test=False, overwrite=False):
	# Create submission.xml
	submission_dict = create_submission_xml(organism=organism, database=database, config_file=config_file, metadata_file=metadata_file, fasta_file=fasta_file, authorset_file=authorset_file)				
	# Extract submission name
	submission_name = submission_dict["submission_name"]
	# Extract submission directory
	submission_dir = submission_dict["submission_dir"]
	# Extract user credentials
	credentials = get_credentials(config_file=config_file, database="NCBI")
	# Extract metadata file
	metadata = get_metadata(organism=organism, database=database, metadata_file=metadata_file)
	# Create a submission ready file
	open(os.path.join(work_dir, "submit.ready"), 'w+').close()
	# Submit sequences to NCBI via FTP Server
	try:
		# Login into NCBI FTP Server
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=credentials["Username"], passwd=credentials["Password"])
		# print out message
		print("\n" + "User authenticated", file=sys.stdout)
		# Submit to test or production ftp server
		if test == True:
			ftp_dir = "Test"
		else:
			ftp_dir = "Production"
		# Create ftp directory if it does not exist
		if ftp_dir not in ftp.nlst():
			ftp.mkd(ftp_dir)
		# CD to to ftp folder
		ftp.cwd(ftp_dir)
		# Create submission directory if it does not exist
		if submission_name not in ftp.nlst():
			ftp.mkd(submission_name)
		# CD to submission folder
		ftp.cwd(submission_name)
		# Check if report.xml exists
		if "report.xml" in ftp.nlst() and overwrite == False:
			print("\n" + "Submission report exists pulling down.", file=sys.stdout)
			report_file = open(os.path.join(submission_dir, "report.xml"), 'wb')
			ftp.retrbinary('RETR report.xml', report_file.write, 262144)
			report_file.close()
		else:
			print("\n" + "Submitting to NCBI-" + database + "-" + organism, file=sys.stdout)
			res = ftp.storlines("STOR " + "submission.xml", os.path.join(submission_dir, "submission.xml"), 'rb')
			if not res.startswith('226 Transfer complete'):
				print('Submission.xml upload failed.', file=sys.stderr)
				sys.exit(1)
			if database in ["BioSample_SRA", "SRA"]:
				for index, row in metadata.iterrows():
					if row["file_location"] == "local":
						local_fasta_files = [os.path.join(fasta_path, g.strip()) for g in row["file_path"].split(",")]
						for file in local_fasta_files:
							res = ftp.storbinary("STOR " + os.path.basename(file.strip()), open(file.strip(), 'rb'))
							if not res.startswith('226 Transfer complete'):
								print('SRA file upload failed. Try again.', file=sys.stderr)
								sys.exit(1)
			elif database == "Genbank":
				res = ftp.storbinary("STOR " + submission_name + ".zip", open(os.path.join(submission_dir, submission_name + ".zip"), 'rb'))
				if not res.startswith('226 Transfer complete'):
					print("Uploading " + os.path.join(submission_dir, submission_name + ".zip") + " failed.", file=sys.stderr)
					sys.exit(1)
			res = ftp.storlines("STOR " + "submit.ready", open(os.path.join(submission_dir, "submit.ready"), 'rb'))
			if not res.startswith('226 Transfer complete'):
				print('submit.ready upload failed.', file=sys.stderr)
				sys.exit(1)
	except ftplib.all_errors as e:
		print("\n" + 'Error:', e, file=sys.stderr)
		sys.exit(1)

# Submit to GISAID
def submit_gisaid(organism, config_file, metadata_file, fasta_file, test=False):
	# Get config file
	config_dict = get_config(config_file=config_file)
	# Get metadata file
	metadata = get_metadata(organism=organism, database="GISAID", metadata_file=metadata_file)
	# Extract GISAID submission information
	try:
		submission_dict = {k.title():v for k,v in config_dict["GISAID"].items()}
	except:
		print("Error: there is no Submission > NCBI information in config file.", file=sys.stderr)
		sys.exit(1)
	else:
		# Extract submission description information
		description_dict = create_submission_description(submission_dict=submission_dict)	
		# Extract description title to store as directory name to dump submission files
		submission_name = description_dict["Title"]
		# Create an output directory
		submission_dir = os.path.join(work_dir, "submission"+str(datetime.now().strftime("%y%m%d")), submission_name, "GISAID", organism)
		# Get fasta file
		if os.path.isfile(fasta_file) == False:
			print("Error: fasta file does not exist at: \n" + fasta_file, file=sys.stderr)
			sys.exit(1)
		# Check each row to make sure each sample name exists in fasta file
		check_fasta_samples(metadata=metadata, fasta_file=fasta_file, sample_colnames=sample_colnames[organism]["GISAID"])
		# Make sure authentication token is valid
		authenticate(database="GISAID", organism=organism, config_file=config_file)
		# Upload the sequences to GISAID
		if test == True:
			print("Performing test submission with Client-Id: " + cid, file=sys.stdout)
			print("If this is not a test Client-Id interrupt the submission immediately.", file=sys.stdout)
			time.sleep(10)
		# Output message
		print("\n" + "Submitting to " + database + "-" + organism, file=sys.stdout)
		# Create log file
		log_file = os.path.join(submission_dir, "gisaid-logfile"+ organism.lower() + ".log")
		# Create failed file
		failed_file = os.path.join(submission_dir, "gisaid-"+ organism.lower() + "-failed.csv")
		# Try submission two times before erroring out
		attempts = 1
		complete = False
		while attempts < 3 and complete == False:
			if organism == "FLU":
				command = subprocess.run("python3 /MIRA/epiflu_cli/__main__.py upload --token gisaid.flu.authtoken --metadata %s --fasta %s --log %s --failed %s --dateformat YYYYMMDD --debug" % (metadata_file, fasta_file, log_file, failed_file),
					env = os.environ.copy(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)
			elif organism == "COV":
				command = subprocess.run("python3 /MIRA/epicov_cli/__main__.py upload --token gisaid.cov.authtoken --metadata %s --fasta %s --log %s --failed %s --dateformat YYYYMMDD --debug" % (metadata_file, fasta_file, log_file, failed_file),
					env = os.environ.copy(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)					
			# Check if uploading is successful
			if command.returncode == 0:
				complete = True
			else:
				print(command.stdout, file=sys.stdout)
				print(command.stderr, file=sys.stdout)
				attempts += 1
				print("Gisaid submission attempt: " + str(attempts), file=sys.stdout)
		# Check if log file is produced
		if complete == True:
			while not os.path.exists(os.path.join(submission_dir, "gisaid-logfile"+ organism.lower() + ".log")):
				time.sleep(10)
			# Read in the log file
			gisaid_submitted_total, gisaid_failed_total = read_gisaid_log(log_file=log_file, failed_file=failed_file)
			# Update submission log
			create_submission_log(database=database, organism=organism, submission_name=submission_name, submission_status="Submitted", gisaid_submitted_total=gisaid_submitted_total, gisaid_failed_total=gisaid_failed_total)
		else:
			print("Submission errored out.", file=sys.stderr)
			sys.exit(1)

# Create submission log csv
def create_submission_log(database, organism, submission_name, submission_status, submisison_type, gisaid_submitted_total="", gisaid_failed_total=""):
	curr_time = datetime.now()
	if os.path.isfile(os.path.join(work_dir, "submission_log.csv")):
		df = pd.read_csv(os.path.join(work_dir, "submission_log.csv"), header = 0, dtype = str, engine = "python", encoding="utf-8", index_col=False)
	else:
		df = pd.DataFrame(columns = ["database","organism","submission_name","submission_date","submission_status", "submisison_type", "GISAID_submitted_total", "GISAID_failed_total"])
	# Check if row exists in log to update instead of write new
	df_partial = df.loc[df["database"]==database, df["organism"]==organism, df['submission_name']==submission_name, df['submisison_type']==submisison_type]
	# Fill in the log with new values
	if df_partial.shape[0] > 0:
		df.loc[df_partial.index.values, 'submission_date'] = curr_time.strftime("%Y-%m-%d")
		df.loc[df_partial.index.values, 'submission_status'] = submission_status
		df.loc[df_partial.index.values, 'GISAID_submitted_total'] = gisaid_submitted_total
		df.loc[df_partial.index.values, 'GISAID_failed_total'] = gisaid_failed_total
	else:
		new_entry = {'database': database,
					 'organism': organism,
					 'submission_name': submission_name,
					 'submission_date': curr_time.strftime("%Y-%m-%d"),
					 'submission_status': submission_status,
					 'submisison_type': submisison_type,
					 'GISAID_submitted_total': gisaid_submitted_total,
					 'GISAID_failed_total': gisaid_failed_total
					}
		df = df.append(new_entry, ignore_index = True)
	df.to_csv(os.path.join(work_dir, "submission_log.csv"), header = True, index = False)		

# # Update submission log
# def get_submision_status(database, organism, submission_name, submission_dir, submisison_type):
# 	# Check if submission log exists
# 	if os.path.isfile(os.path.join(work_dir, "submission_log.csv")):
# 		df = pd.read_csv(os.path.join(work_dir, "submission_log.csv"), header = 0, dtype = str, engine = "python", encoding="utf-8", index_col=False)
# 	else:
# 		print("Error: Either a submission has not been made or submission_log.csv has been moved.", file=sys.stderr)
# 		sys.exit(1)
# 	# Check if row exists in log to update instead of write new file
# 	df_partial = df.loc[df["database"]==database, df["organism"]==organism, df['submission_name']==submission_name, df['submisison_type']==submisison_type]
# 	# Fill in the log with new values
# 	if df_partial.shape[0] > 0:
# 		# Obtain the authentication file
# 		if database == "GISAID":
#
# # Read output log from gisaid submission script
# def read_gisaid_log(organism, log_file, failed_file):
# 	if os.path.isfile(log_file) == False:
# 		print("Error: gisaid log file does not exist at: " + log_file, file=sys.stderr)
# 		sys.exit(1)
# 	else:        
# 		with open(log_file) as f:
#             data = json.load(f)
#             f.close()
#         number_submitted = 0
#         number_failed = 0
#         already_submitted = []
#         for i in data:
#             # Sequence successfully uploaded
#             if i["code"] == "upload_count":
#                 number_submitted = int(i["msg"].strip().split("uploaded: ")[1])
#             # Sequence failed upload
#             elif i["code"] == "failed_count":
#                 number_failed = int(i["msg"].strip().split("failed: ")[1])
#             # Correct number of successfully uploaded for if a sequence fails for already existing
#             elif (i["code"] == "validation_error") and ("\"covv_virus_name\": \"already exists\"" in i["msg"]):
#                 already_submitted.append(i["msg"].split("; validation_error;")[0])
#             elif i["code"] == "epi_isl_id":
#                 continue
#         # Clean failed log
# 		if number_failed == 0 or number_failed == len(already_submitted):
# 			if os.path.exists(failed_file):
# 				print("No failed sequences.\nCleaning up files.", file=sys.stdout)
# 				os.remove(failed_file)
# 		else:
# 				df = pd.read_csv(failed_file, header = 0, dtype = str, engine = "python", encoding="utf-8", index_col=False)
# 				clean_df = df[~df["covv_virus_name"].isin(already_submitted)]
# 				clean_df.to_csv(faile_file, header = True, index = False)
# 				print("Error: Sequences failed please check: " + os.path.join(config_dict["general"]["submission_directory"], unique_name, "gisaid", unique_name + ".log"))
#         number_failed = number_failed - len(already_submitted)
#         number_submitted = number_submitted + len(already_submitted)
#         return str(number_submitted), str(number_failed)
#     else:
#         return "There is no log "

# # Cleans failed meta log if some of the submissions are just already submitted or
# # Removes file if it is empty
# def clean_failed_log(unique_name, number_failed, already_submitted):
#     if number_failed == 0 or number_failed == len(already_submitted):
#         if os.path.exists(os.path.join(config_dict["general"]["submission_directory"], unique_name, "gisaid", unique_name + "_failed_meta.csv")):
#             print("No failed sequences.\nCleaning up files.")
#             os.remove(os.path.join(config_dict["general"]["submission_directory"], unique_name, "gisaid", unique_name + "_failed_meta.csv"))
#     else:
#         df = pd.read_csv(os.path.join(config_dict["general"]["submission_directory"], unique_name, "gisaid", unique_name + "_failed_meta.csv"), header = 0, dtype = str)
#         clean_df = df[~df.covv_virus_name.isin(already_submitted)]
#         clean_df.to_csv(os.path.join(config_dict["general"]["submission_directory"], unique_name, "gisaid", unique_name + "_failed_meta.csv"), header = True, index = False)
#         print("Error: Sequences failed please check: " + os.path.join(config_dict["general"]["submission_directory"], unique_name, "gisaid", unique_name + ".log"))

def args_parser():
	"""
	Argument parser setup and build.
	"""
	# Create submission.xml 
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
									 description='Generate submission.xml for submissions')
	organism_parser = argparse.ArgumentParser(add_help=False)
	config_file_parser = argparse.ArgumentParser(add_help=False)
	database_parser = argparse.ArgumentParser(add_help=False)
	metadata_file_parser = argparse.ArgumentParser(add_help=False)
	fasta_file_parser = argparse.ArgumentParser(add_help=False)
	fasta_path_parser = argparse.ArgumentParser(add_help=False)
	authorset_file_parser = argparse.ArgumentParser(add_help=False)
	test_parser = argparse.ArgumentParser(add_help=False)
	overwrite_parser = argparse.ArgumentParser(add_help=False)

	database_parser.add_argument("--database",
		help="NCBI/GISAID",
		required=True)
	organism_parser.add_argument("--organism",
		help="FLU/COV",
		required=True)		
	config_file_parser.add_argument("--config_file",
		help="Config file",
		required=True)		
	metadata_file_parser.add_argument("--metadata_file",
		help="Metadata file",
		required=True)
	fasta_file_parser.add_argument("--fasta_file",
		help="fasta file",
		required=True)	
	fasta_path_parser.add_argument("--fasta_path",
		help="fasta path",
		required=True)
	authorset_file_parser.add_argument("--authorset_file",
		help="authorser file",
		required=True)
	test_parser.add_argument("--test",
		required=False,
		default="Production",
		action="store_const",
		const="Test")
	overwrite_parser.add_argument("--overwrite",
		help='Overwrite existing submission on NCBI',
		required=False,
		default=False,
		action="store_true")

	subparser_modules = parser.add_subparsers(dest='command')

	auth_module = subparser_modules.add_parser(
		'authenticate', 
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		help='(Interactively) ask for user credentials, generate an authentication-token and store it in the auth-file',
		parents=[database_parser, organism_parser, config_file_parser]
	)

	biosample_module = subparser_modules.add_parser(
		'biosample', 
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		help='Creates submission files and begins automated process of submitting to public databases.',
		parents=[organism_parser, config_file_parser, metadata_file_parser, test_parser, overwrite_parser]
	)

	sra_module = subparser_modules.add_parser(
		'sra', 
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		help='Creates submission files and begins automated process of submitting to public databases.',
		parents=[organism_parser, config_file_parser, metadata_file_parser, fasta_path_parser, test_parser, overwrite_parser]
	)

	biosample_sra_module = subparser_modules.add_parser(
		'biosample_sra', 
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		help='Creates submission files and begins automated process of submitting to public databases.',
		parents=[organism_parser, config_file_parser, metadata_file_parser, fasta_path_parser, test_parser, overwrite_parser]
	)

	genbank_module = subparser_modules.add_parser(
		'genbank', 
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		help='Creates submission files and begins automated process of submitting to public databases.',
		parents=[organism_parser, config_file_parser, metadata_file_parser, fasta_file_parser, authorset_file_parser, test_parser, overwrite_parser]
	)

	gisaid_module = subparser_modules.add_parser(
		'gisaid', 
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		help='Creates submission files and begins automated process of submitting to public databases.',
		parents=[organism_parser, config_file_parser, metadata_file_parser, fasta_file_parser]
	)

	verion_module = subparser_modules.add_parser(
		"version",
		help = "Show version and exit."
	)

	return(parser)

def main():
	"""The main routine"""
	parser = args_parser()
	args = parser.parse_args()

	if args.command == "authenticate":
		authenticate(config_file=args.config_file, database=args.database, organism=args.organism)
		get_execution_time()
	elif args.command == "biosample":
		submission_res = create_submission_xml(organism=args.organism, database="BioSample", config_file=args.config_file, metadata_file=args.metadata_file)
		#submit_ncbi(database=args.database, organism=args.organism, config_file=args.config_file, submission_dict=submission_res["submission_dict"], submission_dir=submission_res["submission_dir"], submission_name=submission_res["submission_name"], test=args.test, overwrite=args.overwrite)
		get_execution_time()
	elif args.command == "sra" or args.command == "biosample_sra":
		submission_dict = create_submission_xml(organism=args.organism, database="BioSample_SRA", config_file=args.config_file, metadata_file=args.metadata_file, fasta_path=args.fasta_path)
		#submit_ncbi(database=args.database, organism=args.organism, config_file=args.config_file, submission_dict=submission_res["submission_dict"], submission_dir=submission_res["submission_dir"], submission_name=submission_res["submission_name"], test=args.test, overwrite=args.overwrite)
		get_execution_time()
	elif args.command == "genbank":
		submission_dict = create_submission_xml(organism=args.organism, database="Genbank", config_file=args.config_file, metadata_file=args.metadata_file, fasta_file=args.fasta_file, authorset_file=args.authorset_file)				
		#submit_ncbi(database=args.database, organism=args.organism, config_file=args.config_file, submission_dict=submission_res["submission_dict"], submission_dir=submission_res["submission_dir"], submission_name=submission_res["submission_name"], test=args.test, overwrite=args.overwrite)
		get_execution_time()
	elif args.command == "gisaid":
		submit_gisaid(organism=args.organism, metadata_file=args.metadata_file, fasta_file=args.fasta_file)
		get_execution_time()
	elif args.command == "version":
	    print("Version: " + version, file=sys.stdout)

if __name__ == "__main__":
	main()
