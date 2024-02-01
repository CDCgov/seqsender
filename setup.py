
# Python Libraries
import pathlib
import sys
from zipfile import ZipFile
import ftplib
import os
import json
import subprocess
import pandas as pd
import shutil
import platform
from urllib.request import urlopen
import gzip
import stat
from io import BytesIO

# Local imports
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import create
import process
import report
import seqsender
import submit

# Get program directory
PROG_DIR = os.path.dirname(os.path.abspath(__file__))

# Create example templates for testing
def create_zip_template(organism, database, submission_dir, submission_name):
	# Create output directory 
	submission_dir = os.path.abspath(submission_dir)
	out_dir = os.path.join(submission_dir, submission_name)
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
	# Create a list of template files to output
	temp_config_file = os.path.join(PROG_DIR, "template", organism, organism.lower()+"_config.yaml")	
	temp_sequence_file = os.path.join(PROG_DIR, "template", organism, organism.lower()+"_sequence.fasta")	
	temp_fastq_1_r1_file = os.path.join(PROG_DIR, "template", organism, organism.lower()+"_fastq_1_R1.fastq.gz")	
	temp_fastq_1_r2_file = os.path.join(PROG_DIR, "template", organism, organism.lower()+"_fastq_1_R2.fastq.gz")
	temp_fastq_2_r1_file = os.path.join(PROG_DIR, "template", organism, organism.lower()+"_fastq_2_R1.fastq.gz")	
	temp_fastq_2_r2_file = os.path.join(PROG_DIR, "template", organism, organism.lower()+"_fastq_2_R2.fastq.gz")
	# Print generating message
	print("\n"+"Generating submission template", file=sys.stdout)
	# Get combined metadata for all given databases
	for i in range(len(database)):
		df = pd.read_csv(os.path.join(PROG_DIR, "template", organism, organism.lower()+"_"+database[i].lower()+"_metadata.csv"), header = 0, dtype = str, engine = "python", encoding="utf-8", index_col=False, na_filter=False)
		if i == 0:
			combined_metadata = df
		else:
			combined_metadata = pd.merge(combined_metadata, df, how='left')
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

def download_table2asn(table2asn_dir):
	# Determine which platform to download table2asn
	if platform.system() == "Windows":
		zip_url = "https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/win64.table2asn.zip"
		with urlopen(zip_url) as zip_response:
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
			with urlopen(zip_url) as zip_response:
				file.write(gzip.decompress(zip_response.read()))
		st = os.stat(table2asn_dir)
		os.chmod(table2asn_dir, st.st_mode | stat.S_IXOTH | stat.S_IRWXU)
	except Exception as error:
		print("Downloading table2asn error", file=sys.stderr)
		print(error, file=sys.stderr)
		sys.exit(1)

