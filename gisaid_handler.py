#!/usr/bin/env python3

###########################    Description    ##################################
# Functions for GISAID input/output handling
################################################################################

import shutil
import subprocess
from typing import Dict, Any
import pandas as pd

# Create directory and files for GISAID submission
def create_gisaid_files(organism: str, database: str, submission_name: str, submission_dir: str, config_dict: Dict[str, Any], metadata: pd.DataFrame) -> None:
	# Get column names for gisaid submission only
	gisaid_df = metadata.filter(regex="^gs-|^collection_date$|^authors").copy()
	gisaid_df.columns = gisaid_df.columns.str.replace("gs-","").str.strip()
	#Add required gisaid fields
	if organism in ["COV", "POX", "ARBO"]:
		gisaid_df = gisaid_df.rename(columns = {"authors": "Authors", "collection_date": "Collection_Date"})
		gisaid_df = gisaid_df.add_prefix((organism.lower() + "_"))
		gisaid_df["submitter"] = config_dict["Username"]
		gisaid_df["fn"] = ""
		first_cols = ["submitter", "fn", (organism.lower() + "_virus_name")]
	elif "FLU" in organism:
		gisaid_df["Isolate_Id"] = ""
		gisaid_df["Segment_Ids"] = ""
		# Rename column names
		gisaid_df = gisaid_df.rename(columns = {"authors": "Authors", "collection_date": "Collection_Date"})
		# Pivot FLU segment names from long form to wide form
		gisaid_df["segment"] = "Seq_Id (" + gisaid_df["segment"].astype(str) + ")"
		group_df = gisaid_df.pivot(index="Isolate_Name", columns="segment", values="sample_name").reset_index()
		gisaid_df = gisaid_df.drop(columns=["sample_name","segment"])
		gisaid_df = gisaid_df.drop_duplicates(keep="first")
		gisaid_df = gisaid_df.merge(group_df, on="Isolate_Name", how="inner", validate="1:1")
		first_cols = ["Isolate_Id","Segment_Ids","Isolate_Name"]
	# Restructure column order
	last_cols = [col for col in gisaid_df.columns if col not in first_cols]
	gisaid_df = gisaid_df[first_cols + last_cols]
	# Create submission files
	gisaid_df.to_csv(os.path.join(submission_files_dir, "metadata.csv"), index=False, sep=",")
	shutil.copy(os.path.join(submission_files_dir, "metadata.csv"), os.path.join(submission_files_dir, "orig_metadata.csv"))
	file_handler.create_fasta(organism=organism, database=["GISAID"], metadata=metadata, submission_files_dir=submission_files_dir)
	shutil.copy(os.path.join(submission_files_dir, "sequence.fsa"), os.path.join(submission_files_dir, "orig_sequence.fsa"))
	print("\n"+"Creating submission files for " + database, file=sys.stdout)
	print("Files are stored at: " + os.path.join(submission_files_dir), file=sys.stdout)

# Read output log from gisaid submission script
def read_gisaid_log(log_file: str, submission_status_file: str) -> pd.DataFrame:
	if os.path.isfile(log_file) is False:
		print("Error: GISAID log file does not exist at: "+log_file, file=sys.stderr)
		print("Error: Either a submission has not been made or log file has been moved.", file=sys.stderr)
		print("Try to re-upload the sequences again.", file=sys.stderr)
		sys.exit(1)
	if os.path.isfile(submission_status_file) is False:
		print("Error: GISAID submission status file does not exist at: "+submission_status_file, file=sys.stderr)
		print("Error: Either a submission has not been made or file has been moved.", file=sys.stderr)
		print("Try to re-upload the sequences again.", file=sys.stderr)
		sys.exit(1)
	# Read in submission status csv
	submission_status = pd.read_csv(submission_status_file, header = 0, dtype = str, engine = "python", encoding="utf-8", index_col=False)
	submission_status = submission_status.fillna("")
	# Read in log file
	with open(log_file) as f:
		while True:
			line = f.readline()
			if not line:
				break
			else:
				# Get sample or segment accession
				if "epi_isl".upper() in line.upper():
					column_name = "gs-sample_name"
					sample_name = list(set(filter(lambda x: (x.upper() in line.upper())==True, submission_status[column_name])))
					accession_id = "epi_isl_id"
					accession = re.search("EPI_ISL_[1-9]+", line)
				elif "epi_id".upper() in line.upper():
					column_name = "gs-sequence_name"
					sample_name = list(set(filter(lambda x: (x.upper() in line.upper())==True, submission_status[column_name])))
					accession_id = "epi_id"
					accession = re.search("EPI[1-9]+", line)
				else:
					continue
				# Get the accession number only
				if accession is not None:
					start = accession.span()[0]
					end = accession.span()[1]
					accession_number = line[start:end]
					sample_message = submission_status.loc[submission_status[column_name].isin(sample_name), "gisaid_message"].astype(str)
					submission_status.loc[submission_status[column_name].isin(sample_name), ("gisaid_accession_" + accession_id)] = accession_number
					submission_status.loc[submission_status[column_name].isin(sample_name), "gisaid_message"] = sample_message + line
				else:
					continue
	# Save submission status df
	submission_status.to_csv(submission_status_file, header = True, index = False)
	not_submitted = submission_status[~submission_status["gisaid_accession_epi_isl_id"].str.contains("EPI", na=False)].copy()
	return not_submitted[["gs-sample_name"]]


# Submit to GISAID
def submit_gisaid(organism: str, database: str, submission_dir: str, submission_name: str, config_dict: Dict[str, Any], gisaid_cli: str, submission_status_file: str, submission_type: str) -> str:
	# Get the directory that stores all submission files
	submission_files_dir = os.path.join(submission_dir, submission_name, "submission_files", database)
	# Gather all required files
	metadata = os.path.join(submission_files_dir, "metadata.csv")
	orig_metadata = os.path.join(submission_files_dir, "orig_metadata.csv")
	fasta = os.path.join(submission_files_dir, "sequence.fsa")
	orig_fasta = os.path.join(submission_files_dir, "orig_sequence.fsa")
	# Extract user credentials (e.g. username, password, client-id)
	process.check_credentials(config_dict=config_dict, database="GISAID")
	# Output message
	print("\n"+"Uploading submission files to GISAID-"+organism, file=sys.stdout)
	print("Performing a '" + submission_type + "' submission with Client-Id: " + config_dict["Client-Id"], file=sys.stdout)
	print("If this is not a '" + submission_type + "' submission, interrupts submission immediately.", file=sys.stdout)
	# Set number of attempt to 3 if erroring out occurs
	attempts = 1
	# Submit to GISAID
	while attempts <= 3:
		print("\n"+"Submission attempt: " + str(attempts), file=sys.stdout)
		# Create a log submission for each attempt
		log_file = os.path.join(submission_files_dir, "gisaid_upload_log_" + str(attempts) + ".txt")
		# If log file exists, removes it
		if os.path.isfile(log_file) == True:
			os.remove(log_file)
		# Upload submission
		command = subprocess.run([gisaid_cli, "upload", "--username", config_dict["Username"], "--password", config_dict["Password"], "--clientid", config_dict["Client-Id"], "--metadata", metadata, "--fasta", fasta, "--log", log_file, "--debug"],
			cwd=submission_files_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		# Check if uploading is successful
		if command.returncode != 0:
			print("Error: upload command error", file=sys.stderr)
			print(command.stdout)
			print(command.stderr)
			sys.exit(1)
		# Check if log file exists
		while not os.path.exists(log_file):
			time.sleep(10)
		# Check submission log to see if all samples are uploaded successfully
		not_submitted_df = process.read_gisaid_log(log_file=log_file, submission_status_file=submission_status_file)
		# If submission completed, no more attempts
		if not_submitted_df.empty:
			print("Uploading successfully", file=sys.stdout)
			print("Status report is stored at: " + submission_status_file, file=sys.stdout)
			print("Log file is stored at: " + submission_files_dir + "/gisaid_upload_log_attempt_" + str(attempts) +  ".txt", file=sys.stdout)
			return "processed-ok"
		else:
			# If submission is not completed, try again
			metadata_df = pd.read_csv(metadata, header = 0, dtype = str, engine = "python", encoding="utf-8", index_col=False)
			if "FLU" in organism:
				column_name = "Isolate_Name"
			elif "COV" in organism:
				column_name = "virus_name"
			metadata_df = metadata_df.merge(not_submitted_df, how="inner", left_on=column_name, right_on="gs-sample_name")
			fasta_names = metadata_df["gs-sequence_name"].tolist()
			metadata_df = metadata_df.drop(columns=["gs-sample_name", "gs-sequence_name"])
			metadata_df.to_csv(orig_metadata, header = True, index = False)
			fasta_dict = []
			with open(orig_fasta, "r") as fsa:
				records = SeqIO.parse(fsa, "fasta")
				for record in records:
					if record.id in fasta_names:
						fasta_dict.append(record)
			with open(fasta, "w+") as fasta_file:
				SeqIO.write(fasta_dict, fasta_file, "fasta")
			attempts += 1
	if not not_submitted_df.empty:
		print("Error: " + str(len(not_submitted_df.index)) + " sample(s) failed to upload to GISAID", file=sys.stderr)
		print("Please check status report at: " + submission_status_file, file=sys.stdout)
		print("Please check log file at: " + submission_files_dir + "/gisaid_upload_log_attempt_{1,2,3}.txt", file=sys.stderr)
		return "Error-Submission-Incomplete"
	else:
		return "processed-ok"

def update_gisaid_files(organism: str, submission_files_dir: str, submission_status_file: str) -> None:
	# Read in the submission status report
	status_df = pd.read_csv(submission_status_file, header = 0, dtype = str, engine = "python", encoding="utf-8", index_col=False)
	# Gather all required files
	metadata = os.path.join(submission_files_dir, "metadata.csv")
	orig_metadata = os.path.join(submission_files_dir, "orig_metadata.csv")
	fasta = os.path.join(submission_files_dir, "sequence.fsa")
	orig_fasta = os.path.join(submission_files_dir, "orig_sequence.fsa")
	# Filter out genbank that has accession number
	genbank_status_df = status_df[status_df["genbank-status"].str.contains("processed-ok", na=False)].copy()
	gisaid_status_df = genbank_status_df[["gs-sample_name", "gs-sequence_name"]]
	# Add required gisaid fields
	metadata_df = pd.read_csv(metadata, header = 0, dtype = str, engine = "python", encoding="utf-8", index_col=False)
	if "FLU" in organism:
		column_name = "Isolate_Name"
	elif "COV" in organism:
		column_name = "virus_name"
	metadata_df = metadata_df.merge(gisaid_status_df, how="inner", left_on=column_name, right_on="gs-sample_name")
	fasta_names = metadata_df["gs-sequence_name"].tolist()
	metadata_df = metadata_df.drop(columns=["gs-sample_name", "gs-sequence_name"])
	metadata_df.to_csv(orig_metadata, header = True, index = False)
	fasta_dict = []
	with open(orig_fasta, "r") as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			if record.id in fasta_names:
				fasta_dict.append(record)
	with open(fasta, "w+") as fasta_file:
		SeqIO.write(fasta_dict, fasta_file, "fasta")
