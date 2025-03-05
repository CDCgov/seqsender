
###########################    Description    ##################################
################################################################################

# Python Libraries
import ftplib
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
from settings import NCBI_FTP_HOST, TABLE2ASN_EMAIL
from logging_handler import CONFIGURED_LOGGER as logger
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
		ftp = ncbi_login(config_dict=config_dict)
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type)
		# Check if report.xml exists
		if "report.xml" in ftp.nlst():
			logger.info("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			logger.success("Download successful.")
			return report_file
		else:
			logger.info("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		logger.error(f"FTP error:\n{e}")
		sys.exit(1)

# Create an empty submit.ready file if it not exists
def create_submit_ready_file(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			logger.error("'submit.ready' upload failed.")
			sys.exit(1)
	except Exception as e:
		if str(e).startswith('Error:550 submit.ready: Permission denied'):
			logger.info("The submission has already been made and is currently processing.")
		else:
			logger.error(f"Unable to upload submit.ready file.\n{e}")
			sys.exit(1)
	return ftp

def ncbi_login(config_dict: Dict[str, Any]):
	logger.debug("Logging into NCBI FTP...")
	ftp = ftplib.FTP(NCBI_FTP_HOST)
	ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	logger.debug("Login successful.")
	return ftp

def ftp_upload_file(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	logger.debug(f"Uploading file: {upload_file}")
	res = ftp.storbinary(f"STOR {upload_name}", open(upload_file, "rb"))
	if not res.startswith('226 Transfer complete'):
		logger.error(f"Uploading {upload_file} failed.")
		sys.exit(1)
	logger.debug(f"'{upload_name}' uploaded successfully.")
	return ftp

def ftp_navigate_to_folder(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		logger.debug(f"Located '{submission_type}' folder.")
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		logger.error("Cannot find submission folder on NCBI FTP site.")
		sys.exit(1)
	else:
		logger.debug("Located 'submit' folder.")
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			logger.debug(f"Located '{submission_type}' folder.")
			ftp.cwd(submission_type)
		else:
			logger.error("Cannot find submission folder on NCBI FTP site.")
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		logger.error("Cannot find submission folder on NCBI FTP site.")
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	logger.debug(f"Submission folder '{folder_name}' created.")
	ftp.cwd(folder_name)
	return ftp

def upload_raw_reads(ftp, submission_dir: str, submission_name: str):
	logger.info("Uploading SRA raw reads...")
	raw_reads_files = os.path.join(submission_dir, "raw_reads_location.txt")
	if os.path.isfile(raw_reads_files) is False:
		logger.error(f"Submission '{submission_name}' is missing raw reads file at: {raw_reads_files}")
		sys.exit(1)
	# Upload SRA files
	with open(raw_reads_files, "r") as file:
		for line in file:
			line = line.strip()
			if line is None or line == "":
				continue
			elif os.path.isfile(line):
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				logger.error(f"Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: {line}")
				sys.exit(1)
	logger.success("SRA raw reads successfully uploaded.")
	return ftp

# Submit to NCBI
def submit_ncbi(database: str, submission_name: str, submission_dir: str, config_dict: Dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	logger.warning(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		logger.info(f"Preparing submission: {ncbi_submission_name}")
		ftp = ncbi_login(config_dict)
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		logger.info(f"Submitting: '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
		logger.success(f"Submission '{submission_name}' uploaded successfully.")
	except ftplib.all_errors as e:
		logger.error(f"FTP error {e}")
		sys.exit(1)

# Send table2asn file through email
def email_table2asn(submission_name: str, submission_dir: str, config_dict: Dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	logger.debug(f"SQN file to email for submission: {sqn_file}")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			logger.warning(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			logger.warning(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			logger.error(f"Submission type '{submission_type}' is not a valid option.")
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
		logger.error(f"Unable to send email automatically. If unable to email, submission can be made manually using the sqn file: {sqn_file}")
		logger.critical(f"Email error:\n{e}")
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
