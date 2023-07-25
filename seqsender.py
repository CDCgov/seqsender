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
#prog_dir = os.path.dirname(os.path.abspath(__file__))
prog_dir = "/seqsender"

# Define version of seqsender
version = "1.0 (Beta)"

# Define current time
STARTTIME = datetime.now()

# Get execution time
def get_execution_time():
    print(f"\nTotal runtime (HRS:MIN:SECS): {str(datetime.now() - STARTTIME)}")

# Create a list of database dictionary
database_dict = {
	"BIOSAMPLE": "BioSample",
	"SRA": "SRA",
	"BIOSAMPLE_SRA": "BioSample_SRA",
	"GENBANK": "Genbank",
	"GISAID": "GISAID"	
}

# Define biosample package 
biosample_packages = {"FLU": "Pathogen.cl.1.0", 
					  "COV": "SARS-CoV-2.cl.1.0"}

# Defing genbank wizard
genbank_wizard = {"FLU": "BankIt_influenza_api", 
				  "COV": "BankIt_SARSCoV2_api"}

# Define gisaid log file
gisaid_logfile = {"FLU": "gisaid.flu.log",
				  "COV": "gisaid.cov.log"}

# Define gisaid failed csv file
gisaid_failedfile = {"FLU": "gisaid.flu.failed.csv",
				     "COV": "gisaid.cov.failed.csv"}
				  
# Define required column names for metadata file
with open(os.path.join(prog_dir, "config", "all_required_colnames.yaml"), "r") as f:
	all_required_colnames = yaml.safe_load(f)	
	f.close()					

# Create gisaid-flu Segment Ids
Segment_Ids = ["HA", "HE", "MP", "NA", "NP", "NS", "P3", "PA", "PB1", "PB2"]

# Define required the sample colunm names to check if samples are unique
sample_colnames = {"FLU": {"BioSample": ["spuid"],
						   "SRA": ["spuid"],
						   "BioSample_SRA": ["spuid"],
						   "Genbank": ["Sequence_ID"],
						   "GISAID": ["Seq_Id ("+x+")" for x in Segment_Ids]
						  },
			 	   "COV": {"BioSample": ["spuid"],
						   "SRA": ["spuid"],
						   "BioSample_SRA": ["spuid"],
						   "Genbank": ["Sequence_ID"],								 
						   "GISAID": ["covv_virus_name"]
						  }
				  }

# NCBI FTP Server
NCBI_FTP_HOST = "ftp-private.ncbi.nlm.nih.gov"

# NCBI API URL
NCBI_API_URL = "https://submit.ncbi.nlm.nih.gov/api/2.0/files/FILE_ID/?format=attachment"

# Read in config file
def get_config(config_file, database):
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
			except:
				print("Error: there is no Submission information in config file.", file=sys.stderr) 
				sys.exit(1)
			else:
				try:
					config_dict = {k.title():v for k,v in config_dict[database].items()}
					return config_dict
				except:
					print("Error: there is no Submission > " + database + " information in config file.", file=sys.stderr)
					sys.exit(1)
		else:
			print("Error: Config file is incorrect. File must be in a valid yaml format.", file=sys.stderr)
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
		# Make sure required fields in gisaid with empty values are filled with 'Unknown'	
		for name in required_colnames:
			if (name in gisaid_unknown_colnames) and any(metadata[name]==""):
				print("Error: The required '" + name + "' field in metatdata file cannot contain any empty values. Please fill empty values with 'Unknown'", file=sys.stderr)
				sys.exit(1)
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
		# Make sure the column contains sample names are unique
		if database != "GISAID" and organism == "FLU":
			unique_colnames = sample_colnames[organism][database]
			for name in unique_colnames:
				unique_names = set(metadata[name])
				if len(unique_names) != metadata.shape[0]:
					print("Error: The required '" + name + "' field in metatdata file must have samples that are unique.", file=sys.stderr)
					sys.exit(1)
		return metadata

# Check user credentials information
def get_credentials(config_file, database):
	# Get config file
	config_dict = get_config(config_file=config_file, database=database)	
	# Create an empty dict to store credentials
	cred_dict = {}
	# Check username		
	if "Username" not in config_dict.keys():
		print("Error: there is no Submission > " + database + " > Username information in config file.", file=sys.stderr)
		sys.exit(1)	
	elif ("Username" in config_dict.keys()) and ((config_dict["Username"] is None) or (config_dict["Username"] == "")):
		print("Error: Submission > " + database + " > Username in the config file cannot be empty.", file=sys.stderr)
		sys.exit(1)	
	elif ("Username" in config_dict.keys()) and ((config_dict["Username"] is not None) and (config_dict["Username"] != "")):
		cred_dict["Username"] = config_dict["Username"]		
	# Check password		
	if "Password" not in config_dict.keys():
		print("Error: there is no Submission > " + database + " > Password information in config file.", file=sys.stderr)
		sys.exit(1)		
	elif ("Password" in config_dict.keys()) and ((config_dict["Password"] is None) or (config_dict["Password"] == "")):
		print("Error: Submission > " + database + " > Password in the config file cannot be empty.", file=sys.stderr)
		sys.exit(1)	
	elif ("Password" in config_dict.keys()) and ((config_dict["Password"] is not None) and (config_dict["Password"] != "")):
		cred_dict["Password"] = config_dict["Password"]	
	# Check client-id if database is GISAID		
	if (database == "GISAID") and ("Client-Id" not in config_dict.keys()):
		print("Error: there is no Submission > " + database + " > Client-Id information in config file.", file=sys.stderr)
		sys.exit(1)	
	elif (database == "GISAID") and ("Client-Id" in config_dict.keys() and ((config_dict["Client-Id"] is None) or (config_dict["Client-Id"] == ""))):
		print("Error: Submission > " + database + " > Client-Id in the config file cannot be empty.", file=sys.stderr)
		sys.exit(1)	
	elif (database == "GISAID") and ("Client-Id" in config_dict.keys() and ((config_dict["Client-Id"] is not None) and (config_dict["Client-Id"] != ""))):
		cred_dict["Client-Id"] = config_dict["Client-Id"]	
	return cred_dict

# Authenticate user credentials
def authenticate(database, organism, config_file):	
	credentials = get_credentials(config_file=config_file, database=database)
	if database == "GISAID":
		log_file = os.path.join("/tmp", gisaid_logfile[organism])
		# If log file exists, removes it
		if os.path.isfile(log_file) == True:
			os.remove(log_file)
		if organism == "FLU":
			command = "python3 %s/epiflu_cli/__main__.py authenticate --token gisaid.flu.authtoken --username %s --password %s --client_id %s --log %s --force > /dev/null" % (prog_dir, credentials["Username"], credentials["Password"], credentials["Client-Id"], log_file)
		elif organism == "COV":
			command = "python3 %s/epicov_cli/__main__.py authenticate --token gisaid.cov.authtoken --username %s --password %s --client_id %s --log %s --force > /dev/null" % (prog_dir, credentials["Username"], credentials["Password"], credentials["Client-Id"], log_file)
		# Check status of the command
		if os.system(command) == 0:
			print("\n" + "User authenticated", file=sys.stdout)
		else:
			print("\n" + "Authentication failed", file=sys.stderr)
			###### Add a check for authenticate date #####
			sys.exit(1)
	elif database == "NCBI":
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
				json.dump(token_file, f)
				f.close()
			# Output authenticated message
			print("\n" + "User authenticated", file=sys.stdout)
		except ftplib.all_errors as e:
			print("\n" + "Authentication failed", file=sys.stderr)
			print("\n" + "Error: " + str(e), file=sys.stderr)
			sys.exit(1)
	return credentials

# Create description section in submission.xml
def create_submission_description(config_dict, database):
	# Check if config has description information
	if "Description" not in config_dict.keys():
		print("Error: there is no Submission > " + database + " > Description information in config file.", file=sys.stderr)
		sys.exit(1)	
	# Check if description has a title (required)
	try:
		title = config_dict["Description"]["Title"]
	except:
		print("Error: there is no Submission > " + database + " > Description > Title in the config file.", file=sys.stderr)
		sys.exit(1)	
	else:
		# Make sure title is not empty or none
		if (title is None) or (title == ""):
			print("Error: Submission > " + database + " > Description > Title in the config file cannot be empty (required).", file=sys.stderr)
			sys.exit(1)	
	# Check if description has organization information (required) for NCBI
	if database == "NCBI":
		try:
			organization = config_dict["Description"]["Organization"]
		except:
			print("Error: there is no Submission > " + database + " > Description > Organization in the config file (required).", file=sys.stderr)
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
	return config_dict["Description"]

# Create action section in submission.xml for BioSample
def create_biosample_submission(metadata, package):
	# Create a dict to combine action for each submission
	combined_action_dict = {"AddData": [""]*metadata.shape[0]}
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
		action_config_dict = get_config(config_file=action_config_path, database="NCBI")	
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
def check_rawread_files(metadata, raw_reads_path):
	# Check fasta local location is local or stored on cloud
	failed_file_location= list(filter(lambda x: (x in ["local", "cloud"])==False, metadata["file_location"]))
	if len(failed_file_location) > 0:
		print("Error: the value of file_location in metatdata file can only be local/cloud.", file=sys.stderr)
		sys.exit(1)	
	# Separate samples stored in local and cloud 
	local_rawread_files = []; cloud_rawread_files = [];
	for index, row in metadata.iterrows():
		# Create an empty list to store all fasta files
		if row["file_location"] == "local":
			local_rawread_files += [os.path.join(raw_reads_path, g.strip()) for g in row["file_path"].split(",")]
		elif row["file_location"] == "cloud":
			cloud_rawread_files += [os.path.join(raw_reads_path, g.strip()) for g in row["file_path"].split(",")]
	# Get a list of fasta paths that does not exist locally
	if len(local_rawread_files) > 0:
		failed_fasta = list(filter(lambda x: os.path.exists(x)==False, local_rawread_files))
		if len(failed_fasta) > 0:
			print("Error: Raw reads files do not exist at " + ", ".join(failed_fasta) + ".", file=sys.stderr)
			sys.exit(1)	
	# Get a list of fasta paths that does not exist on aws cloud
	if len(cloud_rawread_files) > 0:
		failed_fasta = list(filter(lambda x: (subprocess.run("aws s3 ls %s" % x, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).returncode)==1, cloud_rawread_files))
		if len(failed_fasta) > 0:
			print("Error: Raw reads files do not exist at " + ", ".join(failed_fasta) + ".", file=sys.stderr)
			sys.exit(1)	

# Create action section in submission.xml for SRA
def create_sra_submission(metadata, raw_reads_path):	
	# Check fasta files listed in metadata file exist local or on cloud given fasta path
	check_rawread_files(metadata=metadata, raw_reads_path=raw_reads_path)
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
		action_config_dict = get_config(config_file=action_config_path, database="NCBI")	
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
			if row[name] != "":
				command = subprocess.run("grep -c -w '>%s' %s" % (row[name], fasta_file), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
				if command.returncode == 0 and command.stdout is not None:
					if int(command.stdout[:-1]) == 0:
						print("Error: Cannot find " + name + " = " + row[name] + " in " + fasta_file + ". Make sure to start a sample with >sample_name.", file=sys.stderr)
						sys.exit(1)
					elif int(command.stdout[:-1]) > 1:
						print("Error: " + name + " = " + row[name] + " is duplicated in " + fasta_file + ".", file=sys.stderr)
						sys.exit(1)
				else:
					print("Error: Cannot find " + name + " = " + row[name] + " in " + fasta_file + ". Make sure to start a sample with >sample_name.", file=sys.stderr)
					sys.exit(1)

# Create a zip file for genbank submission
def create_genbank_zip(organism, database, metadata, fasta_file, authorset_file, submission_name, submission_dir, comment=False):
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
def create_genbank_submission(organism, database, metadata, fasta_file, authorset_file, submission_name, submission_dir):
	# Check each row to make sure each sample name exists in fasta file
	check_fasta_samples(metadata=metadata, fasta_file=fasta_file, sample_colnames=sample_colnames[organism][database])
	# Create a zip file for genbank submission
	create_genbank_zip(organism=organism, database=database, metadata=metadata, fasta_file=fasta_file, authorset_file=authorset_file, submission_name=submission_name, submission_dir=submission_dir)
	# Read in the action config file
	action_config_path = os.path.join(prog_dir, "config", "genbank_action_config.yaml")
	action_config_dict = get_config(action_config_path, database="NCBI")
	# Fill in the spuid and spuid_namespace information
	action_config_dict["Action"]["AddFiles"]["Identifier"]["SPUID"]["#text"] = submission_name
	action_config_dict["Action"]["AddFiles"]["Identifier"]["SPUID"]["@spuid_namespace"] = "ncbi-"+organism.lower()+"-genbank"
	# Fill in genbank wizard information
	action_config_dict["Action"]["AddFiles"]["Attribute"][0]["#text"] = genbank_wizard[organism]
	# Fill zip file information
	action_config_dict["Action"]["AddFiles"]["File"]["@file_path"] = submission_name + ".zip"
	return action_config_dict["Action"]

def create_submission_xml(organism, database, config_file, metadata_file, fasta_file=None, raw_reads_path=None, authorset_file=None):
	# Get config file
	config_dict = get_config(config_file=config_file, database="NCBI")
	# Get metadata file
	metadata = get_metadata(organism=organism, database=database, metadata_file=metadata_file)
	# Check if fasta file is provided
	if fasta_file is not None and os.path.isfile(fasta_file) == False:
		print("Error: fasta file does not exist at: " + fasta_file, file=sys.stderr)
		sys.exit(1)	
	# Check if fasta path is provided
	if raw_reads_path is not None and os.path.isdir(raw_reads_path) == False:
		print("Error: fasta path does not exist at: " + raw_reads_path, file=sys.stderr)
		sys.exit(1)	
	# Check if authorset file is provided
	if authorset_file is not None and  os.path.isfile(authorset_file) == False:
		print("Error: authorset file does not exist at: " + authorset_file, file=sys.stderr)
		sys.exit(1)	
	# Extract submission description information
	description_dict = create_submission_description(config_dict=config_dict, database="NCBI")
	# Extract description title to store as directory name to dump submission files
	submission_name = description_dict["Title"]
	# Create an output directory
	submission_dir = os.path.join(work_dir, "submission", submission_name, database, organism)
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
		sra_submission_dict = create_sra_submission(metadata=metadata, raw_reads_path=raw_reads_path)
		# Create final submission file
		submission_xml = {"Submission": {"Description": description_dict, "Action2": sra_submission_dict}}	
	elif database == "BioSample_SRA":
		# Create biosample submission dict
		bs_submission_dict = create_biosample_submission(metadata=metadata, package=biosample_packages[organism])
		# Create sra submission dict
		sra_submission_dict = create_sra_submission(metadata=metadata, raw_reads_path=raw_reads_path)
		# Combined biosample and sra as one final submission file
		submission_xml = {"Submission": {"Description": description_dict, "Action1": bs_submission_dict, "Action2": sra_submission_dict}}				
	elif database == "Genbank":
		# Create genebank submission dict
		genbank_submission_dict = create_genbank_submission(organism=organism, database=database, metadata=metadata, fasta_file=fasta_file, authorset_file=authorset_file, submission_name=submission_name, submission_dir=submission_dir)		
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
def submit_ncbi(database, organism, config_file, metadata_file, fasta_file=None, raw_reads_path=None, authorset_file=None, test=True):
	# Create submission.xml
	submission_dict = create_submission_xml(database=database, organism=organism, config_file=config_file, metadata_file=metadata_file, fasta_file=fasta_file, raw_reads_path=raw_reads_path, authorset_file=authorset_file)				
	# Extract submission name
	submission_name = submission_dict["submission_name"]
	# Extract submission directory
	submission_dir = submission_dict["submission_dir"]
	# Extract user credentials
	credentials = get_credentials(config_file=config_file, database="NCBI")
	# Extract metadata file
	metadata = get_metadata(organism=organism, database=database, metadata_file=metadata_file)
	# Create an empty submit.ready file
	open(os.path.join(work_dir, "submit.ready"), 'w+').close()
	# Submit sequences to NCBI via FTP Server
	try:
		# Login into NCBI FTP Server
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=credentials["Username"], passwd=credentials["Password"])
		print("\n" + "User authenticated", file=sys.stdout)
		# Submit to test or production ftp server
		if test == True:
			submission_type = "Test"
		else:
			submission_type = "Production"
		print("\n"+"Submitting to FTP " + submission_type + " server", file=sys.stdout)
		# Create test/production directory if it does not exist
		if submission_type not in ftp.nlst():
			ftp.mkd(submission_type)
		# CD to to test/production folder
		ftp.cwd(submission_type)
		# Create submission directory if it does not exist
		if submission_name not in ftp.nlst():
			ftp.mkd(submission_name)
		# CD to submission folder
		ftp.cwd(submission_name)	
		print("\n" + "Submitting to NCBI-" + database + "-" + organism + " database", file=sys.stdout)
		res = ftp.storlines("STOR " + "submission.xml", open(os.path.join(submission_dir, "submission.xml"), 'rb'))
		if not res.startswith('226 Transfer complete'):
			print('Submission.xml upload failed.', file=sys.stderr)
			sys.exit(1)
		if database in ["BioSample_SRA", "SRA"]:
			for index, row in metadata.iterrows():
				if row["file_location"] == "local":
					local_rawread_files = [os.path.join(raw_reads_path, g.strip()) for g in row["file_path"].split(",")]
					for file in local_rawread_files:
						res = ftp.storbinary("STOR " + os.path.basename(file.strip()), open(file.strip(), 'rb'))
						if not res.startswith('226 Transfer complete'):
							print('SRA file upload failed. Try again.', file=sys.stderr)
							sys.exit(1)
		elif database == "Genbank":
			res = ftp.storbinary("STOR " + submission_name + ".zip", open(os.path.join(submission_dir, submission_name + ".zip"), 'rb'))
			if not res.startswith('226 Transfer complete'):
				print("Uploading " + os.path.join(submission_dir, submission_name + ".zip") + " failed.", file=sys.stderr)
				sys.exit(1)
		res = ftp.storlines("STOR " + "submit.ready", open(os.path.join(work_dir, "submit.ready"), 'rb'))
		if not res.startswith('226 Transfer complete'):
			print('submit.ready upload failed.', file=sys.stderr)
			sys.exit(1)
		# Update submission log
		create_submission_log(database=database, organism=organism, submission_name=submission_name, submission_type=submission_type, submission_status="Submitted", submission_id="", submitted_total=metadata.shape[0], failed_total=0, submission_dir=submission_dir)	
	except ftplib.all_errors as e:
		print("\n" + 'Error:' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to GISAID
def submit_gisaid(organism, database, config_file, metadata_file, fasta_file, test=True):
	# Get config file
	config_dict = get_config(config_file=config_file, database=database)
	# Get metadata file
	metadata = get_metadata(organism=organism, database=database, metadata_file=metadata_file)
	# Get fasta file
	if os.path.isfile(fasta_file) == False:
		print("Error: fasta file does not exist at: \n" + fasta_file, file=sys.stderr)
		sys.exit(1)
	# Extract submission description information
	description_dict = create_submission_description(config_dict=config_dict, database=database)	
	# Extract description title to store as directory name to dump submission files
	submission_name = description_dict["Title"]
	# Create an output directory
	submission_dir = os.path.join(work_dir, "submission", submission_name, database, organism)
	# Check if output directory exists	
	if os.path.exists(submission_dir) == False:
		os.makedirs(submission_dir)	
	# Check each row to make sure each sample name exists in fasta file
	check_fasta_samples(metadata=metadata, fasta_file=fasta_file, sample_colnames=sample_colnames[organism][database])
	# Authenticate to obtain the authentication token 
	credentials = authenticate(database=database, organism=organism, config_file=config_file)
	# Upload the sequences to GISAID
	if test == True:
		submission_type = "Test"
	else:
		submission_type = "Production"
	# Output message
	print("\n"+"Performing a " + submission_type + " submission with Client-Id: " + credentials["Client-Id"], file=sys.stdout)
	print("If this is not a " + submission_type + " submission, interrupts submission immediately.", file=sys.stdout)
	time.sleep(10)
	# Create log file
	log_file = os.path.join(submission_dir, gisaid_logfile[organism])
	# If log file exists, removes it
	if os.path.isfile(log_file) == True:
		os.remove(log_file)
	# Create failed file
	failed_file = os.path.join(submission_dir, gisaid_failedfile[organism])
	# Try submission two times before erroring out
	print("\n" + "Submitting to " + database + "-" + organism, file=sys.stdout)
	attempts = 1
	complete = False
	while attempts < 3 and complete == False:
		if organism == "FLU":
			command = subprocess.run("python3 %s/epiflu_cli/__main__.py upload --token gisaid.flu.authtoken --metadata %s --fasta %s --log %s --debug" % (prog_dir, metadata_file, fasta_file, log_file),
				env = os.environ.copy(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)
		elif organism == "COV":
			command = subprocess.run("python3 %s/epicov_cli/__main__.py upload --token gisaid.cov.authtoken --metadata %s --fasta %s --log %s --failed %s" % (prog_dir, metadata_file, fasta_file, log_file, failed_file),
				env = os.environ.copy(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)					
		# Check if uploading is successful
		if command.returncode == 0:
			complete = True
		else:
			print(command.stdout, file=sys.stdout)
			print(command.stderr, file=sys.stdout)
			attempts += 1
			print("GISAID submission attempt: " + str(attempts), file=sys.stdout)
	# Check if log file is produced
	if complete == True:
		while not os.path.exists(os.path.join(submission_dir, log_file)):
			time.sleep(10)
		# Read in the log file
		submitted_total, failed_total = read_gisaid_log(organism=organism, database=database, log_file=log_file, submission_dir=submission_dir)
		# Update submission log
		create_submission_log(database=database, organism=organism, submission_name=submission_name, submission_type=submission_type, submission_status="Submitted", submission_id="", submitted_total=submitted_total, failed_total=failed_total, submission_dir=submission_dir)
	else:
		print("Submission errored out.", file=sys.stderr)
		print("Please check failed log file at: " + log_file, file=sys.stdout)
		sys.exit(1)

# Read output log from gisaid submission script
def read_gisaid_log(organism, database, log_file, submission_dir):
	if os.path.isfile(log_file) == False:
		print("Error: GISAID log file does not exist at: " + log_file, file=sys.stderr)
		print("Try to re-upload the sequences again.", file=sys.stderr)
		sys.exit(1)
	else:
		# Read in log file
		with open(log_file) as f:
			gisaid_log = json.load(f)
			f.close()
		# Create variables to obtain total submitted and total failed
		number_submitted = 0
		number_failed = 0
		sample_id = []
		gisaid_status = []
		gisaid_accession = []
		gisaid_message = []
		for i in range(len(gisaid_log)):
			if gisaid_log[i]["code"].upper() == "epi_isl_id".upper():
				msg = [g.strip() for g in gisaid_log[i]["msg"].split(";")]
				sample_id.append(msg[0])
				gisaid_accession.append(msg[-1])
				gisaid_status.append("proccessed-ok")
				gisaid_message.append("")
			elif gisaid_log[i]["code"].upper() == "upload_count".upper():
				msg = [g.strip() for g in gisaid_log[i]["msg"].split("submissions uploaded:")]
				number_submitted += int(msg[-1])
			elif gisaid_log[i]["code"].upper() == "failed_count".upper():
				msg = [g.strip() for g in gisaid_log[i]["msg"].split("submissions failed:")] 
				number_failed += int(msg[-1])
			elif gisaid_log[i]["code"].upper() == "validation_error".upper():
				msg = [g.strip() for g in gisaid_log[i]["msg"].split(";")] 
				sample_id.append(msg[0])
				gisaid_status.append("proccessed-error")
				gisaid_message.append(gisaid_log[i]["msg"])
			elif gisaid_log[i]["code"].upper() == "upload_error".upper():
				msg = [g.strip() for g in gisaid_log[i]["msg"].split(";")] 
				sample_id.append(msg[0])
				gisaid_status.append("proccessed-error")
				gisaid_message.append(gisaid_log[i]["msg"])
		# Create the submission status df for GISAID
		gisaid_submission_df = pd.DataFrame().assign(
			sample_id = sample_id, 
			gisaid_status = gisaid_status,
			gisaid_accession = gisaid_accession,
			gisaid_message = gisaid_message
		)
		# Save a copy to submission status df
		gisaid_submission_df.to_csv(os.path.join(submisison_dir, "submission_report_status.csv"), header = True, index = False)
		if number_failed > 0:
			print("Warnings: Some sequences failed to upload to GISAID", file=sys.stdout)
			print("Please check failed log file at: " + log_file, file=sys.stdout)
		return str(number_submitted), str(number_failed)
		
# Create submission log csv
def create_submission_log(database, organism, submission_name, submission_type, submission_status, submission_id, submitted_total, failed_total, submission_dir):
	if os.path.isfile(os.path.join(work_dir, "submission_log.csv")) == True:
		df = pd.read_csv(os.path.join(work_dir, "submission_log.csv"), header = 0, dtype = str, engine = "python", encoding="utf-8", index_col=False)
	else:
		df = pd.DataFrame(columns = ["database", "organism", "submission_name", "submission_type", "submission_date", "submission_status", "submission_id", "submitted_total", "failed_total", "submission_dir", "updated_date"])
	# Fill in the log if exist, otherwise create new
	df_partial = df.loc[(df["database"] == database) & (df["organism"] == organism) & (df['submission_name'] == submission_name) & (df['submission_type'] == submission_type)]
	# Update values or create new enty
	if df_partial.shape[0] > 0:
		df.loc[df_partial.index.values, 'submission_status'] = submission_status
		df.loc[df_partial.index.values, 'submission_id'] = submission_id
		df.loc[df_partial.index.values, 'submitted_total'] = submitted_total
		df.loc[df_partial.index.values, 'failed_total'] = failed_total
		df.loc[df_partial.index.values, 'submission_dir'] = submission_dir
		df.loc[df_partial.index.values, 'updated_date'] = datetime.now().strftime("%Y-%m-%d")		
	else:
		new_entry = {'database': database,
					 'organism': organism,
					 'submission_name': submission_name,
					 'submission_type': submission_type,
					 'submission_date': datetime.now().strftime("%Y-%m-%d"),
					 'submission_status': submission_status,
					 'submission_id': submission_id,
					 'submitted_total': submitted_total,
					 'failed_total': failed_total,
					 'submission_dir': submission_dir,
					 'updated_date': datetime.now().strftime("%Y-%m-%d")
					}
		df.loc[len(df)] = new_entry
	df.to_csv(os.path.join(work_dir, "submission_log.csv"), header = True, index = False)		

# Update submission log
def check_submission_status(database, organism, submission_name, test=True):
	# Check if submission log exists
	if os.path.isfile(os.path.join(work_dir, "submission_log.csv")):
		df = pd.read_csv(os.path.join(work_dir, "submission_log.csv"), header = 0, dtype = str, engine = "python", encoding="utf-8", index_col=False)
	else:
		print("Error: Either a submission has not been made or submission_log.csv has been moved.", file=sys.stderr)
		sys.exit(1)
	# Get the submission type: test/production
	if test == True:
		submission_type = "Test"
	else:
		submission_type = "Production"	
	# Check if submission name is in the submission log
	df_partial = df.loc[(df["database"] == database) & (df["organism"] == organism) & (df['submission_name'] == submission_name) & (df['submission_type'] == submission_type)]
	if df_partial.shape[0] > 0:
		# Obtain submission directory
		submission_dir = df_partial["submission_dir"].values[0]
		# Check if output directory exists	
		if os.path.exists(submission_dir) == False:
			print("Error: Submission directiory does not exit for '"+submission_name+"' at: "+submission_dir, file=sys.stderr)
			print("Error: Either a submission has not been made or folder has been moved.", file=sys.stderr)
			sys.exit(1)
		# Obtain the authentication file
		if database == "GISAID":
			# Get gisaid submission total and failed total
			submitted_total = df_partial["submitted_total"].values[0]
			failed_total = df_partial["failed_total"].values[0]
			# Output message
			print("Submission name: " + submission_name, file=sys.stdout)
			print("Submission database: " + database, file=sys.stdout)
			print("Submission organism: " + organism, file=sys.stdout)
			print("Submission type: " + submission_type, file=sys.stdout)
			print("Submitted total: " + submitted_total, file=sys.stdout)
			print("Failed total: " + failed_total, file=sys.stdout)
		else:
			report_file = get_process_report(database=database, organism=organism, submission_name=submission_name, submission_dir=submission_dir, submission_type=submission_type)
			# Obtain submission.xml path
			submission_file = os.path.join(submission_dir, "submission.xml")
			# Check if submission_xml exists	
			if os.path.exists(submission_file) == False:
				print("Error: submission.xml does not exit for "+submission_name+" at: "+submission_dir, file=sys.stderr)
				print("Error: Either a submission has not been made or submission.xml has been moved.", file=sys.stderr)
				sys.exit(1)
			# Processing the report and output status of the submission
			submission_status, submission_id, submitted_total, failed_total = process_report(database=database, report_file=report_file, submission_file=submission_file)
			# Update submission log
			create_submission_log(database=database, organism=organism, submission_name=submission_name, submission_status=submission_status, submission_id=submission_id, submitted_total=submitted_total, failed_total=failed_total, submission_dir=submission_dir)
			# Output message
			print("Submission name: " + submission_name, file=sys.stdout)
			print("Submission database: " + database, file=sys.stdout)
			print("Submission organism: " + organism, file=sys.stdout)
			print("Submission type: " + submission_type, file=sys.stdout)
			print("Submission id: " + submission_id, file=sys.stdout)
			print("Submission status: " + submission_status, file=sys.stdout)
			print("Submitted total: " + submitted_total, file=sys.stdout)
			print("Failed total: " + failed_total, file=sys.stdout)
	else:
		print("\n" + "Error: " + submission_name + " is not in the submisison log file")

# Get authentication token
def get_token(token_file, organism):
	if os.path.exists(token_file) == False:
		print("Error: Authentication token does not exit at: "+token_file, file=sys.stderr)
		print("Error: Either a submission has not been made or token has been moved.", file=sys.stderr)
		print("Please re-authenticate the database to obtain the token again", file=sys.stderr)
		sys.exit(1)
	else:
		with open(token_file) as f:
			token_dict = json.load(f)
			f.close()    
		# Extract user credentials
		if type(token_dict) is dict:
			try:
				token_dict = {k.title():v for k,v in token_dict[organism.lower()].items()}
			except:
				print("Error: there is no " + organism.lower() + " information in token file.", file=sys.stderr)
				sys.exit(1)
			else:
				# Check username		
				if "Username" not in token_dict.keys():
					print("Error: there is no " + organism.lower() + " > username information in token file.", file=sys.stderr)
					sys.exit(1)	
				elif ("Username" in token_dict.keys()) and ((token_dict["Username"] is None) or (token_dict["Username"] == "")):
					print("Error: " + organism + " > username in the token file cannot be empty.", file=sys.stderr)
					sys.exit(1)		
				# Check password		
				if "Password" not in token_dict.keys():
					print("Error: there is no" + organism.lower() + " > password information in token file.", file=sys.stderr)
					sys.exit(1)		
				elif ("Password" in token_dict.keys()) and ((token_dict["Password"] is None) or (token_dict["Password"] == "")):
					print("Error: " + organism.lower() + " > password in the token file cannot be empty.", file=sys.stderr)
					sys.exit(1)
				return token_dict
		else:
			print("Error: Token file is incorrect. Token file must be in a valid json format.", file=sys.stderr)
			sys.exit(1)

# Process Biosample Report
def get_process_report(database, organism, submission_name, submission_dir, submission_type):
	token_file = os.path.join(work_dir, "ncbi."+organism.lower()+".authtoken")
	credentials = get_token(token_file=token_file, organism=organism)
	# Login into NCBI FTP Server
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=credentials["Username"], passwd=credentials["Password"])
		print("\n" + "User authenticated", file=sys.stdout)
		print("\n" + "Checking "+ submission_type + " submission server", file=sys.stdout)
		# Create test/production directory if it does not exist
		if submission_type not in ftp.nlst():
			ftp.mkd(submission_type)
		# CD to to test/production folder
		ftp.cwd(submission_type)
		# Check if submission name exists
		if submission_name not in ftp.nlst():
			print("There is no submission with the name of '"+ submission_name + "' exists on "+submission_type+" server", file=sys.stderr)
			print("Please try to submit the submission again", file=sys.stderr)
			sys.exit(1)
		# CD to submission folder
		ftp.cwd(submission_name)
		# Check if report.xml exists
		if "report.xml" in ftp.nlst():
			print("\n" + "Pulling down report.xml", file=sys.stdout)
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
				f.close()
			return report_file
		else:
			print("There is no report.xml exists on the FTP server.", file=sys.stdout)
			print("Please try again later.", file=sys.stdout)
			sys.exit(1)
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Read xml report and get status of the submission
def process_report(database, report_file, submission_file):
	# Read in report.xml
	with open(report_file, "r", encoding="utf-8") as f:
		report_xml = f.read()
		f.close()
	# Read in submission.xml
	with open(submission_file, "r", encoding="utf-8") as f:
		submission_xml = f.read()
		f.close()	
	# Convert submission.xml to dictionary	
	submission_dict = xmltodict.parse(submission_xml)
	# Get a list of actions from submission.xml
	submission_action_dict = submission_dict["Submission"]["Action"]
	# Create fields to store submission values for BioSample
	bs_spuid = [""]*len(submission_action_dict)
	biosample_status = [""]*len(submission_action_dict)
	biosample_accession = [""]*len(submission_action_dict)
	biosample_message = [""]*len(submission_action_dict)
	# Create fields to store submission values for SRA
	sra_spuid = [""]*len(submission_action_dict)
	sra_status = [""]*len(submission_action_dict)
	sra_accession = [""]*len(submission_action_dict)
	sra_message = [""]*len(submission_action_dict)
	# Create fields to store submission values for genbank
	genbank_spuid = [""]*len(submission_action_dict)
	genbank_status = [""]*len(submission_action_dict)
	genbank_accession = [""]*len(submission_action_dict)
	genbank_message = [""]*len(submission_action_dict)
	for i in range(len(submission_action_dict)):
		if database == "BioSample":
			bs_spuid[i] = submission_action_dict[i]["AddData"]["Identifier"]["SPUID"]["#text"]
		elif database == "SRA":
			sra_spuid[i] = submission_action_dict[i]["AddFiles"]["Identifier"]["SPUID"]["#text"]
		elif database == "BioSample_SRA":
			try:
				target_db = submission_action_dict["AddData"]["@target_db"]
			except:
				try:
					target_db = submission_action_dict["AddFiles"]["@target_db"]
				except:
					continue
				else:					
					sra_spuid[i] = submission_action_dict[i]["AddFiles"]["Identifier"]["SPUID"]["#text"]
			else:
				bs_spuid[i] = submission_action_dict[i]["AddData"]["Identifier"]["SPUID"]["#text"]
		elif database == "Genbank":
			genbank_spuid[i] = submission_action_dict[i]["AddFiles"]["Identifier"]["SPUID"]["#text"]
	# create the submission df for biosample
	bs_submission_df = pd.DataFrame().assign(
		spuid = bs_spuid, 
		biosample_status = biosample_status,
		biosample_accession = biosample_accession,
		biosample_message = biosample_message
	) 
	# create the submission df for sra
	sra_submission_df = pd.DataFrame().assign(
		spuid = sra_spuid, 
		sra_status = sra_status,
		sra_accession = sra_accession,
		sra_message = sra_message
	) 
	genbank_submission_df = pd.DataFrame().assign(
		Sequence_ID = genbank_spuid, 
		genbank_status = genbank_status,
		genbank_accession = genbank_accession,
		genbank_message = genbank_message
	)
	# Convert xml to dictionary	
	report_dict = xmltodict.parse(report_xml)
	# Get submission status and id from report.xml
	submission_status = report_dict["SubmissionStatus"]["@status"]
	submission_id = report_dict["SubmissionStatus"]["@submission_id"]
	# Get a list of actions from report.xml
	report_action_dict = report_dict["SubmissionStatus"]["Action"]
	for i in range(len(report_action_dict)):
		rp_spuid = re.sub(submission_id+"-", "", report_action_dict[i]["@action_id"])
		rp_spuid = re.sub("_", "/", rp_spuid)
		rp_target_db = report_action_dict["@target_db"]
		if rp_spuid in submission_df["spuid"].lower():
			if rp_target_db == "BioSample":
				bs_submision_df.loc[bs_submision_df["spuid"].lower() == rp_spuid, 'biosample_status'] = report_action_dict[i]["@status"]
				bs_submision_df.loc[bs_submision_df["spuid"].lower() == rp_spuid, 'biosample_message'] = report_action_dict[i]["Response"][0]["Message"]["#text"]
				if report_action_dict[i]["@status"] == "processed-ok":	
					bs_submision_df.loc[bs_submision_df["spuid"].lower() == rp_spuid, 'biosample_accession'] = report_action_dict[i]["Response"][0]["Object"]["@accession"]
			elif rp_target_db == "SRA":
				sra_submision_df.loc[sra_submision_df["spuid"].lower() == rp_spuid, 'sra_status'] = report_action_dict[i]["@status"]
				sra_submision_df.loc[sra_submision_df["spuid"].lower() == rp_spuid, 'sra_message'] = report_action_dict[i]["Response"][0]["Message"]["#text"]
				if report_action_dict[i]["@status"] == "processed-ok":	
					sra_submision_df.loc[sra_submision_df["spuid"].lower == rp_spuid, 'sra_accession'] = report_action_dict[i]["Response"][0]["Object"]["@accession"]
			elif rp_target_db == "GenBank":
				genbank_submision_df.loc[genbank_submision_df["spuid"].lower() == rp_spuid, 'genbank_status'] = report_action_dict[i]["@status"]
				genbank_submision_df.loc[genbank_submision_df["spuid"].lower() == rp_spuid, 'genbank_message'] = report_action_dict[i]["Response"][0]["Message"]["#text"]
				if report_action_dict[i]["@status"] == "processed-ok":	
					genbank_submision_df.loc[genbank_submision_df["spuid"].lower == rp_spuid, 'genbank_accession'] = report_action_dict[i]["Response"][0]["Object"]["@accession"]				
	# Create final submission status df
	if database == "BioSample":
		final_submission_df = bs_submission_df
	elif database == "SRA":
		final_submission_df = sra_submission_df
	elif database == "BioSample_SRA":
		final_submission_df = bs_submission_df.merge(sra_submission_df, left_on='spuid', right_on='spuid') 
	elif database == "Genbank":
		final_submission_df = genbank_submision_df
	# Save final submission status df
	final_submission_df.to_csv(os.path.join(os.path.dirname(submission_file), "submission_report_status.csv"), header = True, index = False)
	# Get number of submitted total and failed total
	submitted_total = final_submission_df.shape[0]
	failed_total = final_submission_df.loc[(final_submission_df["biosample_status"]=="processed-error") | (final_submission_df["sra_status"]=="processed-error") | (final_submission_df["genbank_status"]=="processed-error")].shape[0]
	return submission_status, submission_id, str(number_submitted), str(number_failed)

# Create template for testings
def create_zip_template(database, organism):
	out_dir = os.path.join(work_dir, "template")
	# Check if output directory exists	
	if os.path.exists(out_dir) == False:
		os.makedirs(out_dir)	
	# Create a zip file for genbank submission
	if database == "GISAID":
		with ZipFile(os.path.join(out_dir, database+"-"+organism+"-template.zip"), 'w') as zip:
			zip.write(os.path.join(prog_dir, "data", "metadata", database.lower()+"-"+organism.lower()+"-metadata.csv"), "metadata.csv")
			zip.write(os.path.join(prog_dir, "data", "fasta", database, "-" + organism.lower()+".fasta"), "sequence.fasta")
			zip.write(os.path.join(prog_dir, "data", "config", "default-config.yaml"), "config.yaml")
	elif database == "BioSample":
		with ZipFile(os.path.join(out_dir, database+"-"+organism+"-template.zip"), 'w') as zip:
			zip.write(os.path.join(prog_dir, "data", "metadata", database.lower()+"-"+organism.lower()+"-metadata.csv"), "metadata.csv")
			zip.write(os.path.join(prog_dir, "data", "config", "default-config.yaml"), "config.yaml")
	elif database in ["SRA", "BioSample_SRA"]:
		with ZipFile(os.path.join(out_dir, database+"-"+organism+"-template.zip"), 'w') as zip:
			zip.write(os.path.join(prog_dir, "data", "metadata", database.lower()+"-"+organism.lower()+"-metadata.csv"), "metadata.csv")
			zip.write(os.path.join(prog_dir, "data", "fasta", "sample_1_fastq_R1.fastq.gz"), "sample_1_fastq_R1.fastq.gz")
			zip.write(os.path.join(prog_dir, "data", "fasta", "sample_1_fastq_R2.fastq.gz"), "sample_1_fastq_R2.fastq.gz")
			zip.write(os.path.join(prog_dir, "data", "config", "default-config.yaml"), "config.yaml")
	elif database == "Genbank":
		with ZipFile(os.path.join(out_dir, database+"-"+organism+"-template.zip"), 'w') as zip:
			zip.write(os.path.join(prog_dir, "data", "metadata", database.lower()+"-"+organism.lower()+"-metadata.csv"), "metadata.csv")
			zip.write(os.path.join(prog_dir, "data", "fasta", database, "-" + organism.lower()+".fsa"), "sequence.fsa")
			zip.write(os.path.join(prog_dir, "data", "config", "authorset.sbt"), "authorset.sbt")
			zip.write(os.path.join(prog_dir, "data", "config", "default-config.yaml"), "config.yaml")
	# Waiting for the zip file to write
	while not os.path.exists(os.path.join(out_dir, database+"-"+organism+"-template.zip")):
		time.sleep(10)
	# Print message
	print("\n"+"Generating a zip template for "+database+"-"+organism, file=sys.stdout)
	print("File is stored at: "+os.path.join(out_dir, database+"-"+organism+"-template.zip"), file=sys.stdout)

def args_parser():
	"""
	Argument parser setup and build.
	"""
	# Create submission.xml 
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
									 description='Generate submission.xml for submissions')
	organism_parser = argparse.ArgumentParser(add_help=False)
	config_file_parser = argparse.ArgumentParser(add_help=False)
	auth_database_parser = argparse.ArgumentParser(add_help=False)
	metadata_file_parser = argparse.ArgumentParser(add_help=False)
	fasta_file_parser = argparse.ArgumentParser(add_help=False)
	raw_reads_path_parser = argparse.ArgumentParser(add_help=False)
	authorset_file_parser = argparse.ArgumentParser(add_help=False)
	test_parser = argparse.ArgumentParser(add_help=False)
	check_submission_database_parser = argparse.ArgumentParser(add_help=False)
	submission_name_parser = argparse.ArgumentParser(add_help=False)

	auth_database_parser.add_argument("--database",
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
		help="Fasta file",
		required=True)	
	raw_reads_path_parser.add_argument("--raw_reads_path",
		help="Fasta path",
		required=True)
	authorset_file_parser.add_argument("--authorset_file",
		help="Authorset file",
		required=True)
	test_parser.add_argument("--test",
		help="Whether to perform test submission.",
		required=False,
		action="store_const",
		default=False,
		const=True)
	check_submission_database_parser.add_argument("--database",
		help="BioSample/SRA/BioSample_SRA/Genbank/GISAID",
		required=True)
	submission_name_parser.add_argument("--submission_name",
		help='Name of submission',
		required=True)

	subparser_modules = parser.add_subparsers(dest='command')

	auth_module = subparser_modules.add_parser(
		'authenticate', 
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		help='(Interactively) ask for user credentials, generate an authentication-token and store it in the auth-file',
		parents=[auth_database_parser, organism_parser, config_file_parser]
	)

	biosample_module = subparser_modules.add_parser(
		'biosample', 
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		help='Creates submission files and begins automated process of submitting to public databases.',
		parents=[organism_parser, config_file_parser, metadata_file_parser, test_parser]
	)

	sra_module = subparser_modules.add_parser(
		'sra', 
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		help='Creates submission files and begins automated process of submitting to public databases.',
		parents=[organism_parser, config_file_parser, metadata_file_parser, raw_reads_path_parser, test_parser]
	)

	biosample_sra_module = subparser_modules.add_parser(
		'biosample_sra', 
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		help='Creates submission files and begins automated process of submitting to public databases.',
		parents=[organism_parser, config_file_parser, metadata_file_parser, raw_reads_path_parser, test_parser]
	)

	genbank_module = subparser_modules.add_parser(
		'genbank', 
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		help='Creates submission files and begins automated process of submitting to public databases.',
		parents=[organism_parser, config_file_parser, metadata_file_parser, fasta_file_parser, authorset_file_parser, test_parser]
	)

	gisaid_module = subparser_modules.add_parser(
		'gisaid', 
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		help='Creates submission files and begins automated process of submitting to public databases.',
		parents=[organism_parser, config_file_parser, metadata_file_parser, fasta_file_parser, test_parser]
	)

	update_module = subparser_modules.add_parser(
		'check_submission_status', 
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		help='Check existing process of submissions using submission log.',
		parents=[check_submission_database_parser, organism_parser, submission_name_parser, test_parser]
	)

	template_module = subparser_modules.add_parser(
		'template', 
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		help='Return a zip file containning config file, metadata file, fasta files, etc., that can be used to make a test submission',
		parents=[check_submission_database_parser, organism_parser]
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

	# Make sure organism input all caps
	if args.organism.upper() not in ["FLU", "COV"]:
		print("--organism <> options are FLU/COV")
		sys.exit(1)

	if args.command == "authenticate":
		if args.database.upper() not in ["NCBI", "GISAID"]:
			print("--database <> options are NCBI/GISAID")
			sys.exit(1)
		authenticate(organism=args.organism.upper(), database=args.database.upper(), config_file=args.config_file)  
		get_execution_time()
	elif args.command == "biosample":
		submit_ncbi(organism=args.organism.upper(), database="BioSample", config_file=args.config_file, metadata_file=args.metadata_file, test=args.test)
		get_execution_time()
	elif args.command == "sra":
		submit_ncbi(organism=args.organism.upper(), database="SRA", config_file=args.config_file, metadata_file=args.metadata_file, raw_reads_path=args.raw_reads_path, test=args.test)
		get_execution_time()
	elif args.command == "biosample_sra":
		submit_ncbi(organism=args.organism.upper(), database="BioSample_SRA", config_file=args.config_file, metadata_file=args.metadata_file, raw_reads_path=args.raw_reads_path, test=args.test)
		get_execution_time()
	elif args.command == "genbank":
		submit_ncbi(organism=args.organism.upper(), database="Genbank", config_file=args.config_file, metadata_file=args.metadata_file, fasta_file=args.fasta_file, authorset_file=args.authorset_file, test=args.test)				
		get_execution_time()
	elif args.command == "gisaid":
		submit_gisaid(organism=args.organism.upper(), database="GISAID", config_file=args.config_file, metadata_file=args.metadata_file, fasta_file=args.fasta_file, test=args.test)
		get_execution_time()
	elif args.command == "check_submission_status":
		if args.database.upper() not in database_dict.keys():
			print("--database <> options are BioSample/SRA/BioSample_SRA/Genbank/GISAID")
			sys.exit(1)
		check_submission_status(organism=args.organism.upper(), database=database_dict[args.database.upper()], submission_name=args.submission_name, test=args.test)
		get_execution_time()
	elif args.command == "template":
		if args.database.upper() not in database_dict.keys():
			print("--database <> options are BioSample/SRA/BioSample_SRA/Genbank/GISAID")
			sys.exit(1)
		create_zip_template(organism=args.organism.upper(), database=database_dict[args.database.upper()])
		get_execution_time()
	elif args.command == "version":
	    print("Version: " + version, file=sys.stdout)

if __name__ == "__main__":
	main()
