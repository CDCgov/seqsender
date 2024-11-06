
###########################    Description    ##################################
################################################################################

# Python Libraries
import paramiko
import base64
import os
import sys
import time
import subprocess
import pandas as pd
import smtplib
import xml.etree.ElementTree as ET
import xmltodict
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.application import MIMEApplication
from typing import List, Set, Dict, Any, Union, Tuple, Optional
from settings import NCBI_FTP_HOST, NCBI_FTP_PORT, NCBI_FTP_PUBLIC_KEY, TABLE2ASN_EMAIL

# Local imports
import tools

# Process NCBI Report file
def get_ncbi_report(database: str, submission_name: str, submission_dir: str, config_dict: Dict[str, Any], submission_type: str) -> Optional[str]:
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Login into NCBI FTP Server
	try:
		client, sftp = ncbi_login(config_dict=config_dict)
		sftp_navigate_to_folder(sftp=sftp, folder_name=ncbi_submission_name, submission_type=submission_type)
		# Check if report.xml exists
		if "report.xml" in sftp.listdir():
			print("Downloading report.xml", file=sys.stdout)
			report_file = os.path.join(submission_dir, "report.xml")
			sftp.get("report.xml", report_file)
			return report_file
		else:
			print("The report.xml has not yet been generated.", file=sys.stdout)
			return None
		sftp.close()
		client.close()
	except paramiko.ssh_exception.AuthenticationException:
		print("Error: Authenication failed. Please check your NCBI username and password.", file=sys.stderr)
		sys.exit(1)
	except paramiko.ssh_exception.BadHostKeyException:
		print(f"Error: Bad host key. The host key could not be verified for {NCBI_FTP_HOST}.", file=sys.stderr)
		sys.exit(1)

# Create an empty submit.ready file if it not exists
def create_submit_ready_file(sftp, submission_dir: str) -> None:
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		sftp.put(submit_ready_file)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.", file=sys.stdout)
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)

def ncbi_login(config_dict: Dict[str, Any]):
	client = paramiko.SSHClient()
	client.set_missing_host_key_policy(paramiko.RejectPolicy)
	host_key = paramiko.RSAKey(data = base64.decodebytes(NCBI_FTP_PUBLIC_KEY))
	client.get_host_keys().add(NCBI_FTP_HOST, "ssh-rsa", host_key)
	client.connect(hostname = NCBI_FTP_HOST, port = NCBI_FTP_PORT, username = config_dict["Username"], password = config_dict["Password"], look_for_keys = False, allow_agent = False)
	sftp = client.open_sftp()
	directories = sftp.listdir()
	return client, sftp

def sftp_upload_file(sftp, upload_file: str, upload_name: Optional[str] = None) -> None:
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	sftp.put(upload_file, upload_name)

def sftp_navigate_to_folder(sftp, folder_name: str, submission_type: str, make_folder=False) -> None:
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	directories = sftp.listdir()
	if submission_type in directories:
		sftp.chdir(submission_type)
	elif submission_type not in directories and "submit" not in directories:
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		sftp.chdir("submit")
		directories = sftp.listdir()
		if submission_type in directories:
			sftp.chdir(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	directories = sftp.listdir()
	if not make_folder and folder_name not in directories:
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in directories:
		sftp.mkdir(f"{folder_name}")
	sftp.chdir(folder_name)

def upload_raw_reads(sftp, submission_dir: str, submission_name: str) -> None:
	raw_reads_files = os.path.join(submission_dir, "raw_reads_location.txt")
	if os.path.isfile(raw_reads_files) is False:
		print(f"Error: Submission {submission_name} is missing raw reads file at {raw_reads_files}", file=sys.stderr)
		sys.exit(1)
	# Upload SRA files
	with open(raw_reads_files, "r") as file:
		for line in file:
			line = line.strip()
			if line is None or line == "":
				continue
			elif os.path.isfile(line):
				sftp_upload_file(sftp=sftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)

# Submit to NCBI
def submit_ncbi(database: str, submission_name: str, submission_dir: str, config_dict: Dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.", file=sys.stdout)
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		client, sftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server", file=sys.stdout)
		print(f"Submission name: {ncbi_submission_name}", file=sys.stdout)
		sftp_navigate_to_folder(sftp=sftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'", file=sys.stdout)
		# Upload submission xml
		sftp_upload_file(sftp=sftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			upload_raw_reads(sftp=sftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			sftp_upload_file(sftp=sftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		create_submit_ready_file(sftp=sftp, submission_dir=submission_dir)
		sftp.close()
		client.close()
	except paramiko.ssh_exception.AuthenticationException:
		print("Error: Authenication failed. Please check your NCBI username and password.", file=sys.stderr)
		sys.exit(1)
	except paramiko.ssh_exception.BadHostKeyException:
		print(f"Error: Bad host key. The host key could not be verified for {NCBI_FTP_HOST}.", file=sys.stderr)
		sys.exit(1)

# Send table2asn file through email
def email_table2asn(submission_name: str, submission_dir: str, config_dict: Dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.", file=sys.stdout)
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.", file=sys.stdout)
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"])
		msg['From'] = from_email
		msg['To'] = ", ".join(to_email)
		if len(cc_email) != 0:
			msg['Cc'] = ", ".join(cc_email)
		with open(sqn_file, 'rb') as file_input:
			part = MIMEApplication(file_input.read(), Name=submission_name + ".sqn")
		part['Content-Disposition'] = "attachment; filename=" + submission_name + ".sqn"
		msg.attach(part)
		s = smtplib.SMTP('localhost')
		s.sendmail(from_email, to_email, msg.as_string())
		submission_status = "PROCESSED"
	except Exception as e:
		print("Error: Unable to send mail automatically. If unable to email, submission can be made manually using the sqn file.", file=sys.stderr)
		print(f"sqn_file:{sqn_file}", file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

def standardize_submission_status(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "submitted" in submission_status:
		return "SUBMITTED"
	elif "created" in submission_status:
		return "CREATED"
	elif "queued" in submission_status:
		return "QUEUED"
	elif "processing" in submission_status:
		return "PROCESSING"
	elif "failed" in submission_status:
		return "FAILED"
	elif "processed-ok" in submission_status:
		return "PROCESSED"
	elif "processed-error" in submission_status:
		return "ERROR"
	elif "deleted" in submission_status:
		return "DELETED"
	elif "waiting" in submission_status:
		return "WAITING"
	elif "retried" in submission_status:
		return "RETRIED"
	else:
		return "ERROR"

def process_report_header(report_file: str) -> Tuple[Dict[str, Any], str, str]:
	# Read in report.xml
	tree = ET.parse(report_file)
	root = tree.getroot()
	xmlstr = ET.tostring(root, encoding='utf8', method='xml')
	# Convert xml to dictionary
	report_dict = xmltodict.parse(xmlstr)
	# Get submission status
	try:
		# Get submission status and id from report.xml
		submission_status = report_dict["SubmissionStatus"]["@status"]
	except KeyError:
		submission_status = "SUBMITTED"
	submission_status = standardize_submission_status(submission_status=submission_status)
	# Get submission id
	try:
		submission_id = report_dict["SubmissionStatus"]["@submission_id"]
	except KeyError:
		submission_id = "PENDING"
	return report_dict, submission_status, submission_id
