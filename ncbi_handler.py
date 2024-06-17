
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
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.application import MIMEApplication
from typing import List, Set, Dict, Any, Union

# Local imports
import process

FTP_HOST = "ftp-private.ncbi.nlm.nih.gov"
TABLE2ASN_EMAIL = "gb-admin@ncbi.nlm.nih.gov"

# Process NCBI Report file
def get_ncbi_process_report(database: str, submission_name: str, submission_files_dir: str, config_dict: Dict[str, Any], submission_type: str) -> Union[str, None]:
	# Check user credentials
	process.check_credentials(config_dict=config_dict, database="NCBI")
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Login into NCBI FTP Server
	try:
		ftp = ftplib.FTP(FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
		# Check if submit folder exists
		if "submit" in ftp.nlst():
			ftp.cwd("submit")
			# If submit folder exists check if Production/Test folder exists
			if submission_type not in ftp.nlst():
				print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
				sys.exit(1)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
		ftp.cwd(submission_type)
		# Check if submission name exists
		if ncbi_submission_name not in ftp.nlst():
			print("There is no submission with the name of '"+ ncbi_submission_name +"' on NCBI FTP server.", file=sys.stderr)
			print("Please try the submission again.", file=sys.stderr)
			sys.exit(1)
		# CD to submission folder
		ftp.cwd(ncbi_submission_name)
		# Check if report.xml exists
		if "report.xml" in ftp.nlst():
			print("Pulling down report.xml", file=sys.stdout)
			report_file = os.path.join(submission_files_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.", file=sys.stdout)
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def submit_ncbi(database: str, submission_name: str, submission_dir: str, config_dict: Dict[str, Any], submission_type: str) -> None:
	# Get the directory that stores all submission files
	submission_files_dir = os.path.join(submission_dir, submission_name, "submission_files", database)
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	process.check_credentials(config_dict=config_dict, database="NCBI")
	# Create an empty submit.ready file if it not exists
	submit_ready_file = os.path.join(submission_dir, "submit.ready")
	open(submit_ready_file, 'w+').close()
	# Submit sequences to NCBI via FTP Server
	try:
		# Login into NCBI FTP Server
		ftp = ftplib.FTP(FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
		print("\n"+"Uploading submission files to NCBI-"+database, file=sys.stdout)
		print("Performing a '" + submission_type + "' submission", file=sys.stdout)
		print("If this is not a '" + submission_type + "' submission, interrupts submission immediately.", file=sys.stdout)
		print("\n"+"Connecting to NCBI FTP Server", file=sys.stdout)
		print("Submission name: " + submission_name, file=sys.stdout)
		# Check FTP folder structure either /submit/Production/ or /Production/
		if submission_type not in ftp.nlst():
			# Check if submit folder exists
			if "submit" in ftp.nlst():
				ftp.cwd("submit")
				# If submit folder exists check if Production/Test folder exists
				if submission_type not in ftp.nlst():
					print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
					sys.exit(1)
			else:
				print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
				sys.exit(1)
		ftp.cwd(submission_type)
		# Create submission name directory if it does not exist
		if ncbi_submission_name not in ftp.nlst():
			ftp.mkd(ncbi_submission_name)
		# CD to submission folder
		ftp.cwd(ncbi_submission_name)
		print("Submitting '" + submission_name + "'", file=sys.stdout)
		# Upload submission xml
		res = ftp.storlines("STOR " + "submission.xml", open(os.path.join(submission_files_dir, "submission.xml"), 'rb'))
		if not res.startswith('226 Transfer complete'):
			print('Submission.xml upload failed.', file=sys.stderr)
			sys.exit(1)
		# Upload raw reads
		if "SRA" in database:
			raw_read_location = os.path.join(submission_files_dir, "raw_reads_location.txt")
			if os.path.isfile(raw_read_location) is False:
				print("Error: Submission " + submission_name + " is missing raw reads file at " + raw_read_location, file=sys.stderr)
				sys.exit(1)
			else:
				# Upload SRA files
				with open(raw_read_location, "r") as file:
					for line in file:
						line = line.strip()
						if line is None or line == "":
							continue
						elif os.path.isfile(line):
							res = ftp.storbinary("STOR " + os.path.basename(line), open(line, 'rb'))
							if not res.startswith('226 Transfer complete'):
								print('SRA file upload failed. Try again.', file=sys.stderr)
								sys.exit(1)
						else:
							print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
							sys.exit(1)
		elif "GENBANK" in database:
			res = ftp.storbinary("STOR " + submission_name + ".zip", open(os.path.join(submission_files_dir, submission_name + ".zip"), 'rb'))
			if not res.startswith('226 Transfer complete'):
				print("Uploading " + os.path.join(submission_files_dir, submission_name + ".zip") + " failed.", file=sys.stderr)
				sys.exit(1)
		try:
			res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, 'rb'))
			complete = True
		except Exception as err:
			if str(err).startswith('Error:550 submit.ready: Permission denied'):
				print('The submission has been submitted and currently in pending.', file=sys.stdout)
			else:
				print(err, file=sys.stderr)
				sys.exit(1)
		if (complete == True) and (not res.startswith('226 Transfer complete')):
			print('submit.ready upload failed.', file=sys.stderr)
			sys.exit(1)
		else:
			return
	except ftplib.all_errors as e:
		print("\n" + 'Error:' + str(e), file=sys.stderr)
		sys.exit(1)

# Send table2asn file through email
def sendmail(database: str, submission_name: str, submission_dir: str, config_dict: Dict[str, Any], test: bool) -> str:
	# Create a database subfolder within the submission directory to dump all submission files
	submission_files_dir = os.path.join(submission_dir, submission_name, "submission_files", database)
	# Create submission files directory
	os.makedirs(submission_files_dir, exist_ok=True)
	submission_status = "processed-ok"
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["@email"]
		to_email = []
		cc_email = []
		if test == True:
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["@email"])
		else:
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["@email"])
		if config_dict["Description"]["Organization"]["Submitter"]["@alt_email"]:
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["@alt_email"])
		msg['From'] = from_email
		msg['To'] = ", ".join(to_email)
		if len(cc_email) != 0:
			msg['Cc'] = ", ".join(cc_email)
		with open(os.path.join(submission_dir, submission_name + ".sqn"), 'rb') as file_input:
			part = MIMEApplication(file_input.read(), Name=submission_name + ".sqn")
		part['Content-Disposition'] = "attachment; filename=" + submission_name + ".sqn"
		msg.attach(part)
		s = smtplib.SMTP('localhost')
		s.sendmail(from_email, to_email, msg.as_string())
	except Exception as e:
		print("Error: Unable to send mail automatically. If unable to email, submission can be made manually using the sqn file.", file=sys.stderr)
		print("sqn_file:" + os.path.join(submission_files_dir, submission_name + ".sqn"), file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "processed-ok-email-failure"
	return submission_status
