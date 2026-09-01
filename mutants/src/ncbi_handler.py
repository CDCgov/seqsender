
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
from typing import Any, Union, Optional

# Local imports
import src.tools as tools
import src.setup as setup
from src.settings import NCBI_FTP_HOST, TABLE2ASN_EMAIL
from typing import Annotated
from typing import Callable
from typing import ClassVar

MutantDict = Annotated[dict[str, Callable], "Mutant"] # type: ignore


def _mutmut_trampoline(orig, mutants, call_args, call_kwargs, self_arg = None): # type: ignore
	"""Forward call to original or mutated function, depending on the environment"""
	import os # type: ignore
	mutant_under_test = os.environ['MUTANT_UNDER_TEST'] # type: ignore
	if mutant_under_test == 'fail': # type: ignore
		from mutmut.__main__ import MutmutProgrammaticFailException # type: ignore
		raise MutmutProgrammaticFailException('Failed programmatically')       # type: ignore
	elif mutant_under_test == 'stats': # type: ignore
		from mutmut.__main__ import record_trampoline_hit # type: ignore
		record_trampoline_hit(orig.__module__ + '.' + orig.__name__) # type: ignore
		# (for class methods, orig is bound and thus does not need the explicit self argument)
		result = orig(*call_args, **call_kwargs) # type: ignore
		return result # type: ignore
	prefix = orig.__module__ + '.' + orig.__name__ + '__mutmut_' # type: ignore
	if not mutant_under_test.startswith(prefix): # type: ignore
		result = orig(*call_args, **call_kwargs) # type: ignore
		return result # type: ignore
	mutant_name = mutant_under_test.rpartition('.')[-1] # type: ignore
	if self_arg is not None: # type: ignore
		# call to a class method where self is not bound
		result = mutants[mutant_name](self_arg, *call_args, **call_kwargs) # type: ignore
	else:
		result = mutants[mutant_name](*call_args, **call_kwargs) # type: ignore
	return result # type: ignore

# Process NCBI Report file
def get_ncbi_report(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
	args = [database, submission_name, submission_dir, config_dict, submission_type]# type: ignore
	kwargs = {}# type: ignore
	return _mutmut_trampoline(x_get_ncbi_report__mutmut_orig, x_get_ncbi_report__mutmut_mutants, args, kwargs, None)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_orig(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_1(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
	# Check user credentials
	tools.check_credentials(config_dict=None, database="NCBI")
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Login into NCBI FTP Server
	try:
		ftp = ncbi_login(config_dict=config_dict)
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type)
		# Check if report.xml exists
		if "report.xml" in ftp.nlst():
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_2(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database=None)
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Login into NCBI FTP Server
	try:
		ftp = ncbi_login(config_dict=config_dict)
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type)
		# Check if report.xml exists
		if "report.xml" in ftp.nlst():
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_3(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
	# Check user credentials
	tools.check_credentials(database="NCBI")
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Login into NCBI FTP Server
	try:
		ftp = ncbi_login(config_dict=config_dict)
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type)
		# Check if report.xml exists
		if "report.xml" in ftp.nlst():
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_4(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, )
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Login into NCBI FTP Server
	try:
		ftp = ncbi_login(config_dict=config_dict)
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type)
		# Check if report.xml exists
		if "report.xml" in ftp.nlst():
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_5(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="XXNCBIXX")
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Login into NCBI FTP Server
	try:
		ftp = ncbi_login(config_dict=config_dict)
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type)
		# Check if report.xml exists
		if "report.xml" in ftp.nlst():
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_6(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="ncbi")
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Login into NCBI FTP Server
	try:
		ftp = ncbi_login(config_dict=config_dict)
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type)
		# Check if report.xml exists
		if "report.xml" in ftp.nlst():
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_7(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Create submission name
	ncbi_submission_name = None
	# Login into NCBI FTP Server
	try:
		ftp = ncbi_login(config_dict=config_dict)
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type)
		# Check if report.xml exists
		if "report.xml" in ftp.nlst():
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_8(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Create submission name
	ncbi_submission_name = submission_name + "_" - database
	# Login into NCBI FTP Server
	try:
		ftp = ncbi_login(config_dict=config_dict)
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type)
		# Check if report.xml exists
		if "report.xml" in ftp.nlst():
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_9(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Create submission name
	ncbi_submission_name = submission_name - "_" + database
	# Login into NCBI FTP Server
	try:
		ftp = ncbi_login(config_dict=config_dict)
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type)
		# Check if report.xml exists
		if "report.xml" in ftp.nlst():
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_10(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Create submission name
	ncbi_submission_name = submission_name + "XX_XX" + database
	# Login into NCBI FTP Server
	try:
		ftp = ncbi_login(config_dict=config_dict)
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type)
		# Check if report.xml exists
		if "report.xml" in ftp.nlst():
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_11(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Login into NCBI FTP Server
	try:
		ftp = None
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type)
		# Check if report.xml exists
		if "report.xml" in ftp.nlst():
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_12(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Login into NCBI FTP Server
	try:
		ftp = ncbi_login(config_dict=None)
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type)
		# Check if report.xml exists
		if "report.xml" in ftp.nlst():
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_13(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Login into NCBI FTP Server
	try:
		ftp = ncbi_login(config_dict=config_dict)
		ftp = None
		# Check if report.xml exists
		if "report.xml" in ftp.nlst():
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_14(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Login into NCBI FTP Server
	try:
		ftp = ncbi_login(config_dict=config_dict)
		ftp = ftp_navigate_to_folder(ftp=None, folder_name=ncbi_submission_name, submission_type=submission_type)
		# Check if report.xml exists
		if "report.xml" in ftp.nlst():
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_15(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Login into NCBI FTP Server
	try:
		ftp = ncbi_login(config_dict=config_dict)
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=None, submission_type=submission_type)
		# Check if report.xml exists
		if "report.xml" in ftp.nlst():
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_16(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Login into NCBI FTP Server
	try:
		ftp = ncbi_login(config_dict=config_dict)
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=None)
		# Check if report.xml exists
		if "report.xml" in ftp.nlst():
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_17(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Login into NCBI FTP Server
	try:
		ftp = ncbi_login(config_dict=config_dict)
		ftp = ftp_navigate_to_folder(folder_name=ncbi_submission_name, submission_type=submission_type)
		# Check if report.xml exists
		if "report.xml" in ftp.nlst():
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_18(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Login into NCBI FTP Server
	try:
		ftp = ncbi_login(config_dict=config_dict)
		ftp = ftp_navigate_to_folder(ftp=ftp, submission_type=submission_type)
		# Check if report.xml exists
		if "report.xml" in ftp.nlst():
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_19(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Login into NCBI FTP Server
	try:
		ftp = ncbi_login(config_dict=config_dict)
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, )
		# Check if report.xml exists
		if "report.xml" in ftp.nlst():
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_20(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Login into NCBI FTP Server
	try:
		ftp = ncbi_login(config_dict=config_dict)
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type)
		# Check if report.xml exists
		if "XXreport.xmlXX" in ftp.nlst():
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_21(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Login into NCBI FTP Server
	try:
		ftp = ncbi_login(config_dict=config_dict)
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type)
		# Check if report.xml exists
		if "REPORT.XML" in ftp.nlst():
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_22(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Login into NCBI FTP Server
	try:
		ftp = ncbi_login(config_dict=config_dict)
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type)
		# Check if report.xml exists
		if "report.xml" not in ftp.nlst():
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_23(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print(None)
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_24(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("XXDownloading report.xmlXX")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_25(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_26(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("DOWNLOADING REPORT.XML")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_27(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = None
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_28(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(None, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_29(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, None)
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_30(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join("report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_31(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, )
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_32(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "XXreport.xmlXX")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_33(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "REPORT.XML")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_34(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(None, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_35(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, None) as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_36(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open('wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_37(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, ) as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_38(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'XXwbXX') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_39(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'WB') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_40(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary(None, f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_41(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', None, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_42(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, None)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_43(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary(f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_44(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_45(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, )
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_46(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('XXRETR report.xmlXX', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_47(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('retr report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_48(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR REPORT.XML', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_49(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262145)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_50(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print(None)
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_51(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("XXThe report.xml has not yet been generated.XX")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_52(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("the report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_53(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("THE REPORT.XML HAS NOT YET BEEN GENERATED.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_54(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print(None, file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_55(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=None)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_56(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print(file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_57(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), )
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_58(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " - str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_59(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" - "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_60(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("XX\nXX" + "Error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_61(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "XXError: XX" + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_62(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "error: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_63(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "ERROR: " + str(e), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_64(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(None), file=sys.stderr)
		sys.exit(1)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_65(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(None)

# Process NCBI Report file
def x_get_ncbi_report__mutmut_66(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> Optional[str]:
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
			print("Downloading report.xml")
			report_file = os.path.join(submission_dir, "report.xml")
			with open(report_file, 'wb') as f:
				ftp.retrbinary('RETR report.xml', f.write, 262144)
			return report_file
		else:
			print("The report.xml has not yet been generated.")
			return None
	except ftplib.all_errors as e:
		print("\n" + "Error: " + str(e), file=sys.stderr)
		sys.exit(2)

x_get_ncbi_report__mutmut_mutants : ClassVar[MutantDict] = { # type: ignore
'x_get_ncbi_report__mutmut_1': x_get_ncbi_report__mutmut_1, 
    'x_get_ncbi_report__mutmut_2': x_get_ncbi_report__mutmut_2, 
    'x_get_ncbi_report__mutmut_3': x_get_ncbi_report__mutmut_3, 
    'x_get_ncbi_report__mutmut_4': x_get_ncbi_report__mutmut_4, 
    'x_get_ncbi_report__mutmut_5': x_get_ncbi_report__mutmut_5, 
    'x_get_ncbi_report__mutmut_6': x_get_ncbi_report__mutmut_6, 
    'x_get_ncbi_report__mutmut_7': x_get_ncbi_report__mutmut_7, 
    'x_get_ncbi_report__mutmut_8': x_get_ncbi_report__mutmut_8, 
    'x_get_ncbi_report__mutmut_9': x_get_ncbi_report__mutmut_9, 
    'x_get_ncbi_report__mutmut_10': x_get_ncbi_report__mutmut_10, 
    'x_get_ncbi_report__mutmut_11': x_get_ncbi_report__mutmut_11, 
    'x_get_ncbi_report__mutmut_12': x_get_ncbi_report__mutmut_12, 
    'x_get_ncbi_report__mutmut_13': x_get_ncbi_report__mutmut_13, 
    'x_get_ncbi_report__mutmut_14': x_get_ncbi_report__mutmut_14, 
    'x_get_ncbi_report__mutmut_15': x_get_ncbi_report__mutmut_15, 
    'x_get_ncbi_report__mutmut_16': x_get_ncbi_report__mutmut_16, 
    'x_get_ncbi_report__mutmut_17': x_get_ncbi_report__mutmut_17, 
    'x_get_ncbi_report__mutmut_18': x_get_ncbi_report__mutmut_18, 
    'x_get_ncbi_report__mutmut_19': x_get_ncbi_report__mutmut_19, 
    'x_get_ncbi_report__mutmut_20': x_get_ncbi_report__mutmut_20, 
    'x_get_ncbi_report__mutmut_21': x_get_ncbi_report__mutmut_21, 
    'x_get_ncbi_report__mutmut_22': x_get_ncbi_report__mutmut_22, 
    'x_get_ncbi_report__mutmut_23': x_get_ncbi_report__mutmut_23, 
    'x_get_ncbi_report__mutmut_24': x_get_ncbi_report__mutmut_24, 
    'x_get_ncbi_report__mutmut_25': x_get_ncbi_report__mutmut_25, 
    'x_get_ncbi_report__mutmut_26': x_get_ncbi_report__mutmut_26, 
    'x_get_ncbi_report__mutmut_27': x_get_ncbi_report__mutmut_27, 
    'x_get_ncbi_report__mutmut_28': x_get_ncbi_report__mutmut_28, 
    'x_get_ncbi_report__mutmut_29': x_get_ncbi_report__mutmut_29, 
    'x_get_ncbi_report__mutmut_30': x_get_ncbi_report__mutmut_30, 
    'x_get_ncbi_report__mutmut_31': x_get_ncbi_report__mutmut_31, 
    'x_get_ncbi_report__mutmut_32': x_get_ncbi_report__mutmut_32, 
    'x_get_ncbi_report__mutmut_33': x_get_ncbi_report__mutmut_33, 
    'x_get_ncbi_report__mutmut_34': x_get_ncbi_report__mutmut_34, 
    'x_get_ncbi_report__mutmut_35': x_get_ncbi_report__mutmut_35, 
    'x_get_ncbi_report__mutmut_36': x_get_ncbi_report__mutmut_36, 
    'x_get_ncbi_report__mutmut_37': x_get_ncbi_report__mutmut_37, 
    'x_get_ncbi_report__mutmut_38': x_get_ncbi_report__mutmut_38, 
    'x_get_ncbi_report__mutmut_39': x_get_ncbi_report__mutmut_39, 
    'x_get_ncbi_report__mutmut_40': x_get_ncbi_report__mutmut_40, 
    'x_get_ncbi_report__mutmut_41': x_get_ncbi_report__mutmut_41, 
    'x_get_ncbi_report__mutmut_42': x_get_ncbi_report__mutmut_42, 
    'x_get_ncbi_report__mutmut_43': x_get_ncbi_report__mutmut_43, 
    'x_get_ncbi_report__mutmut_44': x_get_ncbi_report__mutmut_44, 
    'x_get_ncbi_report__mutmut_45': x_get_ncbi_report__mutmut_45, 
    'x_get_ncbi_report__mutmut_46': x_get_ncbi_report__mutmut_46, 
    'x_get_ncbi_report__mutmut_47': x_get_ncbi_report__mutmut_47, 
    'x_get_ncbi_report__mutmut_48': x_get_ncbi_report__mutmut_48, 
    'x_get_ncbi_report__mutmut_49': x_get_ncbi_report__mutmut_49, 
    'x_get_ncbi_report__mutmut_50': x_get_ncbi_report__mutmut_50, 
    'x_get_ncbi_report__mutmut_51': x_get_ncbi_report__mutmut_51, 
    'x_get_ncbi_report__mutmut_52': x_get_ncbi_report__mutmut_52, 
    'x_get_ncbi_report__mutmut_53': x_get_ncbi_report__mutmut_53, 
    'x_get_ncbi_report__mutmut_54': x_get_ncbi_report__mutmut_54, 
    'x_get_ncbi_report__mutmut_55': x_get_ncbi_report__mutmut_55, 
    'x_get_ncbi_report__mutmut_56': x_get_ncbi_report__mutmut_56, 
    'x_get_ncbi_report__mutmut_57': x_get_ncbi_report__mutmut_57, 
    'x_get_ncbi_report__mutmut_58': x_get_ncbi_report__mutmut_58, 
    'x_get_ncbi_report__mutmut_59': x_get_ncbi_report__mutmut_59, 
    'x_get_ncbi_report__mutmut_60': x_get_ncbi_report__mutmut_60, 
    'x_get_ncbi_report__mutmut_61': x_get_ncbi_report__mutmut_61, 
    'x_get_ncbi_report__mutmut_62': x_get_ncbi_report__mutmut_62, 
    'x_get_ncbi_report__mutmut_63': x_get_ncbi_report__mutmut_63, 
    'x_get_ncbi_report__mutmut_64': x_get_ncbi_report__mutmut_64, 
    'x_get_ncbi_report__mutmut_65': x_get_ncbi_report__mutmut_65, 
    'x_get_ncbi_report__mutmut_66': x_get_ncbi_report__mutmut_66
}
x_get_ncbi_report__mutmut_orig.__name__ = 'x_get_ncbi_report'

# Create an empty submit.ready file if it not exists
def create_submit_ready_file(ftp, submission_dir: str):
	args = [ftp, submission_dir]# type: ignore
	kwargs = {}# type: ignore
	return _mutmut_trampoline(x_create_submit_ready_file__mutmut_orig, x_create_submit_ready_file__mutmut_mutants, args, kwargs, None)

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_orig(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_1(ftp, submission_dir: str):
	try:
		submit_ready_file = None
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_2(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(None, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_3(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, None)
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_4(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join("submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_5(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, )
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_6(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "XXsubmit.readyXX")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_7(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "SUBMIT.READY")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_8(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(None, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_9(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, None).close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_10(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open('w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_11(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, ).close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_12(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'XXw+XX').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_13(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'W+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_14(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = None
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_15(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines(None, open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_16(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", None)
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_17(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines(open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_18(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", )
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_19(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " - "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_20(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("XXSTOR XX" + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_21(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("stor " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_22(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "XXsubmit.readyXX", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_23(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "SUBMIT.READY", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_24(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(None, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_25(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, None))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_26(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open("rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_27(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, ))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_28(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "XXrbXX"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_29(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "RB"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_30(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_31(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith(None):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_32(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('XX226 Transfer completeXX'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_33(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_34(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 TRANSFER COMPLETE'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_35(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print(None, file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_36(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=None)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_37(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print(file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_38(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", )
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_39(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("XXError: submit.ready upload failed.XX", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_40(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_41(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("ERROR: SUBMIT.READY UPLOAD FAILED.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_42(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(None)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_43(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(2)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_44(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith(None):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_45(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(None).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_46(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('XXError:550 submit.ready: Permission deniedXX'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_47(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('error:550 submit.ready: permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_48(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('ERROR:550 SUBMIT.READY: PERMISSION DENIED'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_49(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print(None)
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_50(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("XXThe submission has already been made and is currently processing.XX")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_51(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("the submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_52(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("THE SUBMISSION HAS ALREADY BEEN MADE AND IS CURRENTLY PROCESSING.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_53(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(None, file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_54(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=None)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_55(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(file=sys.stderr)
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_56(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", )
			sys.exit(1)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_57(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(None)
	return ftp

# Create an empty submit.ready file if it not exists
def x_create_submit_ready_file__mutmut_58(ftp, submission_dir: str):
	try:
		submit_ready_file = os.path.join(submission_dir, "submit.ready")
		open(submit_ready_file, 'w+').close()
		res = ftp.storlines("STOR " + "submit.ready", open(submit_ready_file, "rb"))
		if not res.startswith('226 Transfer complete'):
			print("Error: submit.ready upload failed.", file=sys.stderr)
			sys.exit(1)
	except Exception as err:
		if str(err).startswith('Error:550 submit.ready: Permission denied'):
			print("The submission has already been made and is currently processing.")
		else:
			print(f"Error: Unable to upload submit.ready file. {err}", file=sys.stderr)
			sys.exit(2)
	return ftp

x_create_submit_ready_file__mutmut_mutants : ClassVar[MutantDict] = { # type: ignore
'x_create_submit_ready_file__mutmut_1': x_create_submit_ready_file__mutmut_1, 
    'x_create_submit_ready_file__mutmut_2': x_create_submit_ready_file__mutmut_2, 
    'x_create_submit_ready_file__mutmut_3': x_create_submit_ready_file__mutmut_3, 
    'x_create_submit_ready_file__mutmut_4': x_create_submit_ready_file__mutmut_4, 
    'x_create_submit_ready_file__mutmut_5': x_create_submit_ready_file__mutmut_5, 
    'x_create_submit_ready_file__mutmut_6': x_create_submit_ready_file__mutmut_6, 
    'x_create_submit_ready_file__mutmut_7': x_create_submit_ready_file__mutmut_7, 
    'x_create_submit_ready_file__mutmut_8': x_create_submit_ready_file__mutmut_8, 
    'x_create_submit_ready_file__mutmut_9': x_create_submit_ready_file__mutmut_9, 
    'x_create_submit_ready_file__mutmut_10': x_create_submit_ready_file__mutmut_10, 
    'x_create_submit_ready_file__mutmut_11': x_create_submit_ready_file__mutmut_11, 
    'x_create_submit_ready_file__mutmut_12': x_create_submit_ready_file__mutmut_12, 
    'x_create_submit_ready_file__mutmut_13': x_create_submit_ready_file__mutmut_13, 
    'x_create_submit_ready_file__mutmut_14': x_create_submit_ready_file__mutmut_14, 
    'x_create_submit_ready_file__mutmut_15': x_create_submit_ready_file__mutmut_15, 
    'x_create_submit_ready_file__mutmut_16': x_create_submit_ready_file__mutmut_16, 
    'x_create_submit_ready_file__mutmut_17': x_create_submit_ready_file__mutmut_17, 
    'x_create_submit_ready_file__mutmut_18': x_create_submit_ready_file__mutmut_18, 
    'x_create_submit_ready_file__mutmut_19': x_create_submit_ready_file__mutmut_19, 
    'x_create_submit_ready_file__mutmut_20': x_create_submit_ready_file__mutmut_20, 
    'x_create_submit_ready_file__mutmut_21': x_create_submit_ready_file__mutmut_21, 
    'x_create_submit_ready_file__mutmut_22': x_create_submit_ready_file__mutmut_22, 
    'x_create_submit_ready_file__mutmut_23': x_create_submit_ready_file__mutmut_23, 
    'x_create_submit_ready_file__mutmut_24': x_create_submit_ready_file__mutmut_24, 
    'x_create_submit_ready_file__mutmut_25': x_create_submit_ready_file__mutmut_25, 
    'x_create_submit_ready_file__mutmut_26': x_create_submit_ready_file__mutmut_26, 
    'x_create_submit_ready_file__mutmut_27': x_create_submit_ready_file__mutmut_27, 
    'x_create_submit_ready_file__mutmut_28': x_create_submit_ready_file__mutmut_28, 
    'x_create_submit_ready_file__mutmut_29': x_create_submit_ready_file__mutmut_29, 
    'x_create_submit_ready_file__mutmut_30': x_create_submit_ready_file__mutmut_30, 
    'x_create_submit_ready_file__mutmut_31': x_create_submit_ready_file__mutmut_31, 
    'x_create_submit_ready_file__mutmut_32': x_create_submit_ready_file__mutmut_32, 
    'x_create_submit_ready_file__mutmut_33': x_create_submit_ready_file__mutmut_33, 
    'x_create_submit_ready_file__mutmut_34': x_create_submit_ready_file__mutmut_34, 
    'x_create_submit_ready_file__mutmut_35': x_create_submit_ready_file__mutmut_35, 
    'x_create_submit_ready_file__mutmut_36': x_create_submit_ready_file__mutmut_36, 
    'x_create_submit_ready_file__mutmut_37': x_create_submit_ready_file__mutmut_37, 
    'x_create_submit_ready_file__mutmut_38': x_create_submit_ready_file__mutmut_38, 
    'x_create_submit_ready_file__mutmut_39': x_create_submit_ready_file__mutmut_39, 
    'x_create_submit_ready_file__mutmut_40': x_create_submit_ready_file__mutmut_40, 
    'x_create_submit_ready_file__mutmut_41': x_create_submit_ready_file__mutmut_41, 
    'x_create_submit_ready_file__mutmut_42': x_create_submit_ready_file__mutmut_42, 
    'x_create_submit_ready_file__mutmut_43': x_create_submit_ready_file__mutmut_43, 
    'x_create_submit_ready_file__mutmut_44': x_create_submit_ready_file__mutmut_44, 
    'x_create_submit_ready_file__mutmut_45': x_create_submit_ready_file__mutmut_45, 
    'x_create_submit_ready_file__mutmut_46': x_create_submit_ready_file__mutmut_46, 
    'x_create_submit_ready_file__mutmut_47': x_create_submit_ready_file__mutmut_47, 
    'x_create_submit_ready_file__mutmut_48': x_create_submit_ready_file__mutmut_48, 
    'x_create_submit_ready_file__mutmut_49': x_create_submit_ready_file__mutmut_49, 
    'x_create_submit_ready_file__mutmut_50': x_create_submit_ready_file__mutmut_50, 
    'x_create_submit_ready_file__mutmut_51': x_create_submit_ready_file__mutmut_51, 
    'x_create_submit_ready_file__mutmut_52': x_create_submit_ready_file__mutmut_52, 
    'x_create_submit_ready_file__mutmut_53': x_create_submit_ready_file__mutmut_53, 
    'x_create_submit_ready_file__mutmut_54': x_create_submit_ready_file__mutmut_54, 
    'x_create_submit_ready_file__mutmut_55': x_create_submit_ready_file__mutmut_55, 
    'x_create_submit_ready_file__mutmut_56': x_create_submit_ready_file__mutmut_56, 
    'x_create_submit_ready_file__mutmut_57': x_create_submit_ready_file__mutmut_57, 
    'x_create_submit_ready_file__mutmut_58': x_create_submit_ready_file__mutmut_58
}
x_create_submit_ready_file__mutmut_orig.__name__ = 'x_create_submit_ready_file'

def ncbi_login(config_dict: dict[str, Any], crash_on_error: bool = False):
	args = [config_dict, crash_on_error]# type: ignore
	kwargs = {}# type: ignore
	return _mutmut_trampoline(x_ncbi_login__mutmut_orig, x_ncbi_login__mutmut_mutants, args, kwargs, None)

def x_ncbi_login__mutmut_orig(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_1(config_dict: dict[str, Any], crash_on_error: bool = True):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_2(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = None
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_3(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(None)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_4(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=None, passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_5(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=None)
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_6(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_7(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], )
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_8(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["XXUsernameXX"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_9(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_10(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["USERNAME"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_11(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["XXPasswordXX"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_12(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_13(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["PASSWORD"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_14(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(None, file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_15(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=None)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_16(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_17(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", )
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_18(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(None)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_19(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(2)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_20(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print(None, file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_21(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=None)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_22(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print(file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_23(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", )
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_24(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("XXError unable to connect to FTP site. Running network test...XX", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_25(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("error unable to connect to ftp site. running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_26(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("ERROR UNABLE TO CONNECT TO FTP SITE. RUNNING NETWORK TEST...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_27(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=None)
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_28(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["XXNCBIXX"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_29(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["ncbi"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_30(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(None, file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_31(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=None)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_32(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_33(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", )
		sys.exit(1)
	return ftp

def x_ncbi_login__mutmut_34(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(None)
	return ftp

def x_ncbi_login__mutmut_35(config_dict: dict[str, Any], crash_on_error: bool = False):
	try:
		ftp = ftplib.FTP(NCBI_FTP_HOST)
		ftp.login(user=config_dict["Username"], passwd=config_dict["Password"])
	except ftplib.error_perm as err:
		print(f"Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException {err}", file=sys.stderr)
		if crash_on_error:
			sys.exit(1)
	except Exception as err:
		print("Error unable to connect to FTP site. Running network test...", file=sys.stderr)
		setup.test_internet_connection(databases=["NCBI"])
		print(f"Exception: {err}", file=sys.stderr)
		sys.exit(2)
	return ftp

x_ncbi_login__mutmut_mutants : ClassVar[MutantDict] = { # type: ignore
'x_ncbi_login__mutmut_1': x_ncbi_login__mutmut_1, 
    'x_ncbi_login__mutmut_2': x_ncbi_login__mutmut_2, 
    'x_ncbi_login__mutmut_3': x_ncbi_login__mutmut_3, 
    'x_ncbi_login__mutmut_4': x_ncbi_login__mutmut_4, 
    'x_ncbi_login__mutmut_5': x_ncbi_login__mutmut_5, 
    'x_ncbi_login__mutmut_6': x_ncbi_login__mutmut_6, 
    'x_ncbi_login__mutmut_7': x_ncbi_login__mutmut_7, 
    'x_ncbi_login__mutmut_8': x_ncbi_login__mutmut_8, 
    'x_ncbi_login__mutmut_9': x_ncbi_login__mutmut_9, 
    'x_ncbi_login__mutmut_10': x_ncbi_login__mutmut_10, 
    'x_ncbi_login__mutmut_11': x_ncbi_login__mutmut_11, 
    'x_ncbi_login__mutmut_12': x_ncbi_login__mutmut_12, 
    'x_ncbi_login__mutmut_13': x_ncbi_login__mutmut_13, 
    'x_ncbi_login__mutmut_14': x_ncbi_login__mutmut_14, 
    'x_ncbi_login__mutmut_15': x_ncbi_login__mutmut_15, 
    'x_ncbi_login__mutmut_16': x_ncbi_login__mutmut_16, 
    'x_ncbi_login__mutmut_17': x_ncbi_login__mutmut_17, 
    'x_ncbi_login__mutmut_18': x_ncbi_login__mutmut_18, 
    'x_ncbi_login__mutmut_19': x_ncbi_login__mutmut_19, 
    'x_ncbi_login__mutmut_20': x_ncbi_login__mutmut_20, 
    'x_ncbi_login__mutmut_21': x_ncbi_login__mutmut_21, 
    'x_ncbi_login__mutmut_22': x_ncbi_login__mutmut_22, 
    'x_ncbi_login__mutmut_23': x_ncbi_login__mutmut_23, 
    'x_ncbi_login__mutmut_24': x_ncbi_login__mutmut_24, 
    'x_ncbi_login__mutmut_25': x_ncbi_login__mutmut_25, 
    'x_ncbi_login__mutmut_26': x_ncbi_login__mutmut_26, 
    'x_ncbi_login__mutmut_27': x_ncbi_login__mutmut_27, 
    'x_ncbi_login__mutmut_28': x_ncbi_login__mutmut_28, 
    'x_ncbi_login__mutmut_29': x_ncbi_login__mutmut_29, 
    'x_ncbi_login__mutmut_30': x_ncbi_login__mutmut_30, 
    'x_ncbi_login__mutmut_31': x_ncbi_login__mutmut_31, 
    'x_ncbi_login__mutmut_32': x_ncbi_login__mutmut_32, 
    'x_ncbi_login__mutmut_33': x_ncbi_login__mutmut_33, 
    'x_ncbi_login__mutmut_34': x_ncbi_login__mutmut_34, 
    'x_ncbi_login__mutmut_35': x_ncbi_login__mutmut_35
}
x_ncbi_login__mutmut_orig.__name__ = 'x_ncbi_login'

def ftp_upload_file(ftp, upload_file: str, upload_name: Optional[str] = None):
	args = [ftp, upload_file, upload_name]# type: ignore
	kwargs = {}# type: ignore
	return _mutmut_trampoline(x_ftp_upload_file__mutmut_orig, x_ftp_upload_file__mutmut_mutants, args, kwargs, None)

def x_ftp_upload_file__mutmut_orig(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	res = ftp.storbinary(f"STOR {upload_name}", open(upload_file, "rb"))
	if not res.startswith('226 Transfer complete'):
		print(f"Error: Uploading {upload_file} failed.", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_1(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is not None:
		upload_name = os.path.basename(upload_file)
	res = ftp.storbinary(f"STOR {upload_name}", open(upload_file, "rb"))
	if not res.startswith('226 Transfer complete'):
		print(f"Error: Uploading {upload_file} failed.", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_2(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = None
	res = ftp.storbinary(f"STOR {upload_name}", open(upload_file, "rb"))
	if not res.startswith('226 Transfer complete'):
		print(f"Error: Uploading {upload_file} failed.", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_3(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(None)
	res = ftp.storbinary(f"STOR {upload_name}", open(upload_file, "rb"))
	if not res.startswith('226 Transfer complete'):
		print(f"Error: Uploading {upload_file} failed.", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_4(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	res = None
	if not res.startswith('226 Transfer complete'):
		print(f"Error: Uploading {upload_file} failed.", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_5(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	res = ftp.storbinary(None, open(upload_file, "rb"))
	if not res.startswith('226 Transfer complete'):
		print(f"Error: Uploading {upload_file} failed.", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_6(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	res = ftp.storbinary(f"STOR {upload_name}", None)
	if not res.startswith('226 Transfer complete'):
		print(f"Error: Uploading {upload_file} failed.", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_7(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	res = ftp.storbinary(open(upload_file, "rb"))
	if not res.startswith('226 Transfer complete'):
		print(f"Error: Uploading {upload_file} failed.", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_8(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	res = ftp.storbinary(f"STOR {upload_name}", )
	if not res.startswith('226 Transfer complete'):
		print(f"Error: Uploading {upload_file} failed.", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_9(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	res = ftp.storbinary(f"STOR {upload_name}", open(None, "rb"))
	if not res.startswith('226 Transfer complete'):
		print(f"Error: Uploading {upload_file} failed.", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_10(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	res = ftp.storbinary(f"STOR {upload_name}", open(upload_file, None))
	if not res.startswith('226 Transfer complete'):
		print(f"Error: Uploading {upload_file} failed.", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_11(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	res = ftp.storbinary(f"STOR {upload_name}", open("rb"))
	if not res.startswith('226 Transfer complete'):
		print(f"Error: Uploading {upload_file} failed.", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_12(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	res = ftp.storbinary(f"STOR {upload_name}", open(upload_file, ))
	if not res.startswith('226 Transfer complete'):
		print(f"Error: Uploading {upload_file} failed.", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_13(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	res = ftp.storbinary(f"STOR {upload_name}", open(upload_file, "XXrbXX"))
	if not res.startswith('226 Transfer complete'):
		print(f"Error: Uploading {upload_file} failed.", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_14(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	res = ftp.storbinary(f"STOR {upload_name}", open(upload_file, "RB"))
	if not res.startswith('226 Transfer complete'):
		print(f"Error: Uploading {upload_file} failed.", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_15(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	res = ftp.storbinary(f"STOR {upload_name}", open(upload_file, "rb"))
	if res.startswith('226 Transfer complete'):
		print(f"Error: Uploading {upload_file} failed.", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_16(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	res = ftp.storbinary(f"STOR {upload_name}", open(upload_file, "rb"))
	if not res.startswith(None):
		print(f"Error: Uploading {upload_file} failed.", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_17(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	res = ftp.storbinary(f"STOR {upload_name}", open(upload_file, "rb"))
	if not res.startswith('XX226 Transfer completeXX'):
		print(f"Error: Uploading {upload_file} failed.", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_18(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	res = ftp.storbinary(f"STOR {upload_name}", open(upload_file, "rb"))
	if not res.startswith('226 transfer complete'):
		print(f"Error: Uploading {upload_file} failed.", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_19(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	res = ftp.storbinary(f"STOR {upload_name}", open(upload_file, "rb"))
	if not res.startswith('226 TRANSFER COMPLETE'):
		print(f"Error: Uploading {upload_file} failed.", file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_20(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	res = ftp.storbinary(f"STOR {upload_name}", open(upload_file, "rb"))
	if not res.startswith('226 Transfer complete'):
		print(None, file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_21(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	res = ftp.storbinary(f"STOR {upload_name}", open(upload_file, "rb"))
	if not res.startswith('226 Transfer complete'):
		print(f"Error: Uploading {upload_file} failed.", file=None)
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_22(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	res = ftp.storbinary(f"STOR {upload_name}", open(upload_file, "rb"))
	if not res.startswith('226 Transfer complete'):
		print(file=sys.stderr)
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_23(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	res = ftp.storbinary(f"STOR {upload_name}", open(upload_file, "rb"))
	if not res.startswith('226 Transfer complete'):
		print(f"Error: Uploading {upload_file} failed.", )
		sys.exit(1)
	return ftp

def x_ftp_upload_file__mutmut_24(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	res = ftp.storbinary(f"STOR {upload_name}", open(upload_file, "rb"))
	if not res.startswith('226 Transfer complete'):
		print(f"Error: Uploading {upload_file} failed.", file=sys.stderr)
		sys.exit(None)
	return ftp

def x_ftp_upload_file__mutmut_25(ftp, upload_file: str, upload_name: Optional[str] = None):
	if upload_name is None:
		upload_name = os.path.basename(upload_file)
	res = ftp.storbinary(f"STOR {upload_name}", open(upload_file, "rb"))
	if not res.startswith('226 Transfer complete'):
		print(f"Error: Uploading {upload_file} failed.", file=sys.stderr)
		sys.exit(2)
	return ftp

x_ftp_upload_file__mutmut_mutants : ClassVar[MutantDict] = { # type: ignore
'x_ftp_upload_file__mutmut_1': x_ftp_upload_file__mutmut_1, 
    'x_ftp_upload_file__mutmut_2': x_ftp_upload_file__mutmut_2, 
    'x_ftp_upload_file__mutmut_3': x_ftp_upload_file__mutmut_3, 
    'x_ftp_upload_file__mutmut_4': x_ftp_upload_file__mutmut_4, 
    'x_ftp_upload_file__mutmut_5': x_ftp_upload_file__mutmut_5, 
    'x_ftp_upload_file__mutmut_6': x_ftp_upload_file__mutmut_6, 
    'x_ftp_upload_file__mutmut_7': x_ftp_upload_file__mutmut_7, 
    'x_ftp_upload_file__mutmut_8': x_ftp_upload_file__mutmut_8, 
    'x_ftp_upload_file__mutmut_9': x_ftp_upload_file__mutmut_9, 
    'x_ftp_upload_file__mutmut_10': x_ftp_upload_file__mutmut_10, 
    'x_ftp_upload_file__mutmut_11': x_ftp_upload_file__mutmut_11, 
    'x_ftp_upload_file__mutmut_12': x_ftp_upload_file__mutmut_12, 
    'x_ftp_upload_file__mutmut_13': x_ftp_upload_file__mutmut_13, 
    'x_ftp_upload_file__mutmut_14': x_ftp_upload_file__mutmut_14, 
    'x_ftp_upload_file__mutmut_15': x_ftp_upload_file__mutmut_15, 
    'x_ftp_upload_file__mutmut_16': x_ftp_upload_file__mutmut_16, 
    'x_ftp_upload_file__mutmut_17': x_ftp_upload_file__mutmut_17, 
    'x_ftp_upload_file__mutmut_18': x_ftp_upload_file__mutmut_18, 
    'x_ftp_upload_file__mutmut_19': x_ftp_upload_file__mutmut_19, 
    'x_ftp_upload_file__mutmut_20': x_ftp_upload_file__mutmut_20, 
    'x_ftp_upload_file__mutmut_21': x_ftp_upload_file__mutmut_21, 
    'x_ftp_upload_file__mutmut_22': x_ftp_upload_file__mutmut_22, 
    'x_ftp_upload_file__mutmut_23': x_ftp_upload_file__mutmut_23, 
    'x_ftp_upload_file__mutmut_24': x_ftp_upload_file__mutmut_24, 
    'x_ftp_upload_file__mutmut_25': x_ftp_upload_file__mutmut_25
}
x_ftp_upload_file__mutmut_orig.__name__ = 'x_ftp_upload_file'

def ftp_navigate_to_folder(ftp, folder_name: str, submission_type: str, make_folder=False):
	args = [ftp, folder_name, submission_type, make_folder]# type: ignore
	kwargs = {}# type: ignore
	return _mutmut_trampoline(x_ftp_navigate_to_folder__mutmut_orig, x_ftp_navigate_to_folder__mutmut_mutants, args, kwargs, None)

def x_ftp_navigate_to_folder__mutmut_orig(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_1(ftp, folder_name: str, submission_type: str, make_folder=True):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_2(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = None
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_3(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type not in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_4(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(None)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_5(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() or "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_6(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_7(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "XXsubmitXX" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_8(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "SUBMIT" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_9(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_10(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print(None, file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_11(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=None)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_12(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print(file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_13(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", )
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_14(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("XXError: Cannot find submission folder on NCBI FTP site.XX", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_15(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("error: cannot find submission folder on ncbi ftp site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_16(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("ERROR: CANNOT FIND SUBMISSION FOLDER ON NCBI FTP SITE.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_17(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(None)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_18(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(2)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_19(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd(None)
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_20(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("XXsubmitXX")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_21(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("SUBMIT")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_22(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type not in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_23(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(None)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_24(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print(None, file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_25(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=None)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_26(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print(file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_27(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", )
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_28(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("XXError: Cannot find submission folder on NCBI FTP site.XX", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_29(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("error: cannot find submission folder on ncbi ftp site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_30(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("ERROR: CANNOT FIND SUBMISSION FOLDER ON NCBI FTP SITE.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_31(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(None)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_32(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(2)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_33(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder or folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_34(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_35(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_36(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print(None, file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_37(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=None)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_38(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print(file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_39(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", )
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_40(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("XXError: Cannot find submission folder on NCBI FTP site.XX", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_41(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("error: cannot find submission folder on ncbi ftp site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_42(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("ERROR: CANNOT FIND SUBMISSION FOLDER ON NCBI FTP SITE.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_43(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(None)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_44(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(2)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_45(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder or folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_46(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_47(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(None)
	ftp.cwd(folder_name)
	return ftp

def x_ftp_navigate_to_folder__mutmut_48(ftp, folder_name: str, submission_type: str, make_folder=False):
	# Ensure correct punctuation for folders
	submission_type = submission_type.capitalize()
	# Check FTP folder structure either /submit/Production/ or /Production/
	if submission_type in ftp.nlst():
		ftp.cwd(submission_type)
	elif submission_type not in ftp.nlst() and "submit" not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	else:
		ftp.cwd("submit")
		if submission_type in ftp.nlst():
			ftp.cwd(submission_type)
		else:
			print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
			sys.exit(1)
	# Check if submission folder exists / can be created
	if not make_folder and folder_name not in ftp.nlst():
		print("Error: Cannot find submission folder on NCBI FTP site.", file=sys.stderr)
		sys.exit(1)
	elif make_folder and folder_name not in ftp.nlst():
		ftp.mkd(folder_name)
	ftp.cwd(None)
	return ftp

x_ftp_navigate_to_folder__mutmut_mutants : ClassVar[MutantDict] = { # type: ignore
'x_ftp_navigate_to_folder__mutmut_1': x_ftp_navigate_to_folder__mutmut_1, 
    'x_ftp_navigate_to_folder__mutmut_2': x_ftp_navigate_to_folder__mutmut_2, 
    'x_ftp_navigate_to_folder__mutmut_3': x_ftp_navigate_to_folder__mutmut_3, 
    'x_ftp_navigate_to_folder__mutmut_4': x_ftp_navigate_to_folder__mutmut_4, 
    'x_ftp_navigate_to_folder__mutmut_5': x_ftp_navigate_to_folder__mutmut_5, 
    'x_ftp_navigate_to_folder__mutmut_6': x_ftp_navigate_to_folder__mutmut_6, 
    'x_ftp_navigate_to_folder__mutmut_7': x_ftp_navigate_to_folder__mutmut_7, 
    'x_ftp_navigate_to_folder__mutmut_8': x_ftp_navigate_to_folder__mutmut_8, 
    'x_ftp_navigate_to_folder__mutmut_9': x_ftp_navigate_to_folder__mutmut_9, 
    'x_ftp_navigate_to_folder__mutmut_10': x_ftp_navigate_to_folder__mutmut_10, 
    'x_ftp_navigate_to_folder__mutmut_11': x_ftp_navigate_to_folder__mutmut_11, 
    'x_ftp_navigate_to_folder__mutmut_12': x_ftp_navigate_to_folder__mutmut_12, 
    'x_ftp_navigate_to_folder__mutmut_13': x_ftp_navigate_to_folder__mutmut_13, 
    'x_ftp_navigate_to_folder__mutmut_14': x_ftp_navigate_to_folder__mutmut_14, 
    'x_ftp_navigate_to_folder__mutmut_15': x_ftp_navigate_to_folder__mutmut_15, 
    'x_ftp_navigate_to_folder__mutmut_16': x_ftp_navigate_to_folder__mutmut_16, 
    'x_ftp_navigate_to_folder__mutmut_17': x_ftp_navigate_to_folder__mutmut_17, 
    'x_ftp_navigate_to_folder__mutmut_18': x_ftp_navigate_to_folder__mutmut_18, 
    'x_ftp_navigate_to_folder__mutmut_19': x_ftp_navigate_to_folder__mutmut_19, 
    'x_ftp_navigate_to_folder__mutmut_20': x_ftp_navigate_to_folder__mutmut_20, 
    'x_ftp_navigate_to_folder__mutmut_21': x_ftp_navigate_to_folder__mutmut_21, 
    'x_ftp_navigate_to_folder__mutmut_22': x_ftp_navigate_to_folder__mutmut_22, 
    'x_ftp_navigate_to_folder__mutmut_23': x_ftp_navigate_to_folder__mutmut_23, 
    'x_ftp_navigate_to_folder__mutmut_24': x_ftp_navigate_to_folder__mutmut_24, 
    'x_ftp_navigate_to_folder__mutmut_25': x_ftp_navigate_to_folder__mutmut_25, 
    'x_ftp_navigate_to_folder__mutmut_26': x_ftp_navigate_to_folder__mutmut_26, 
    'x_ftp_navigate_to_folder__mutmut_27': x_ftp_navigate_to_folder__mutmut_27, 
    'x_ftp_navigate_to_folder__mutmut_28': x_ftp_navigate_to_folder__mutmut_28, 
    'x_ftp_navigate_to_folder__mutmut_29': x_ftp_navigate_to_folder__mutmut_29, 
    'x_ftp_navigate_to_folder__mutmut_30': x_ftp_navigate_to_folder__mutmut_30, 
    'x_ftp_navigate_to_folder__mutmut_31': x_ftp_navigate_to_folder__mutmut_31, 
    'x_ftp_navigate_to_folder__mutmut_32': x_ftp_navigate_to_folder__mutmut_32, 
    'x_ftp_navigate_to_folder__mutmut_33': x_ftp_navigate_to_folder__mutmut_33, 
    'x_ftp_navigate_to_folder__mutmut_34': x_ftp_navigate_to_folder__mutmut_34, 
    'x_ftp_navigate_to_folder__mutmut_35': x_ftp_navigate_to_folder__mutmut_35, 
    'x_ftp_navigate_to_folder__mutmut_36': x_ftp_navigate_to_folder__mutmut_36, 
    'x_ftp_navigate_to_folder__mutmut_37': x_ftp_navigate_to_folder__mutmut_37, 
    'x_ftp_navigate_to_folder__mutmut_38': x_ftp_navigate_to_folder__mutmut_38, 
    'x_ftp_navigate_to_folder__mutmut_39': x_ftp_navigate_to_folder__mutmut_39, 
    'x_ftp_navigate_to_folder__mutmut_40': x_ftp_navigate_to_folder__mutmut_40, 
    'x_ftp_navigate_to_folder__mutmut_41': x_ftp_navigate_to_folder__mutmut_41, 
    'x_ftp_navigate_to_folder__mutmut_42': x_ftp_navigate_to_folder__mutmut_42, 
    'x_ftp_navigate_to_folder__mutmut_43': x_ftp_navigate_to_folder__mutmut_43, 
    'x_ftp_navigate_to_folder__mutmut_44': x_ftp_navigate_to_folder__mutmut_44, 
    'x_ftp_navigate_to_folder__mutmut_45': x_ftp_navigate_to_folder__mutmut_45, 
    'x_ftp_navigate_to_folder__mutmut_46': x_ftp_navigate_to_folder__mutmut_46, 
    'x_ftp_navigate_to_folder__mutmut_47': x_ftp_navigate_to_folder__mutmut_47, 
    'x_ftp_navigate_to_folder__mutmut_48': x_ftp_navigate_to_folder__mutmut_48
}
x_ftp_navigate_to_folder__mutmut_orig.__name__ = 'x_ftp_navigate_to_folder'

def upload_raw_reads(ftp, submission_dir: str, submission_name: str):
	args = [ftp, submission_dir, submission_name]# type: ignore
	kwargs = {}# type: ignore
	return _mutmut_trampoline(x_upload_raw_reads__mutmut_orig, x_upload_raw_reads__mutmut_mutants, args, kwargs, None)

def x_upload_raw_reads__mutmut_orig(ftp, submission_dir: str, submission_name: str):
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
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_1(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = None
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
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_2(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(None, "raw_reads_location.txt")
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
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_3(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, None)
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
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_4(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join("raw_reads_location.txt")
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
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_5(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, )
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
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_6(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, "XXraw_reads_location.txtXX")
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
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_7(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, "RAW_READS_LOCATION.TXT")
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
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_8(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, "raw_reads_location.txt")
	if os.path.isfile(None) is False:
		print(f"Error: Submission {submission_name} is missing raw reads file at {raw_reads_files}", file=sys.stderr)
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
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_9(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, "raw_reads_location.txt")
	if os.path.isfile(raw_reads_files) is not False:
		print(f"Error: Submission {submission_name} is missing raw reads file at {raw_reads_files}", file=sys.stderr)
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
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_10(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, "raw_reads_location.txt")
	if os.path.isfile(raw_reads_files) is True:
		print(f"Error: Submission {submission_name} is missing raw reads file at {raw_reads_files}", file=sys.stderr)
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
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_11(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, "raw_reads_location.txt")
	if os.path.isfile(raw_reads_files) is False:
		print(None, file=sys.stderr)
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
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_12(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, "raw_reads_location.txt")
	if os.path.isfile(raw_reads_files) is False:
		print(f"Error: Submission {submission_name} is missing raw reads file at {raw_reads_files}", file=None)
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
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_13(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, "raw_reads_location.txt")
	if os.path.isfile(raw_reads_files) is False:
		print(file=sys.stderr)
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
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_14(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, "raw_reads_location.txt")
	if os.path.isfile(raw_reads_files) is False:
		print(f"Error: Submission {submission_name} is missing raw reads file at {raw_reads_files}", )
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
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_15(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, "raw_reads_location.txt")
	if os.path.isfile(raw_reads_files) is False:
		print(f"Error: Submission {submission_name} is missing raw reads file at {raw_reads_files}", file=sys.stderr)
		sys.exit(None)
	# Upload SRA files
	with open(raw_reads_files, "r") as file:
		for line in file:
			line = line.strip()
			if line is None or line == "":
				continue
			elif os.path.isfile(line):
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_16(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, "raw_reads_location.txt")
	if os.path.isfile(raw_reads_files) is False:
		print(f"Error: Submission {submission_name} is missing raw reads file at {raw_reads_files}", file=sys.stderr)
		sys.exit(2)
	# Upload SRA files
	with open(raw_reads_files, "r") as file:
		for line in file:
			line = line.strip()
			if line is None or line == "":
				continue
			elif os.path.isfile(line):
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_17(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, "raw_reads_location.txt")
	if os.path.isfile(raw_reads_files) is False:
		print(f"Error: Submission {submission_name} is missing raw reads file at {raw_reads_files}", file=sys.stderr)
		sys.exit(1)
	# Upload SRA files
	with open(None, "r") as file:
		for line in file:
			line = line.strip()
			if line is None or line == "":
				continue
			elif os.path.isfile(line):
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_18(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, "raw_reads_location.txt")
	if os.path.isfile(raw_reads_files) is False:
		print(f"Error: Submission {submission_name} is missing raw reads file at {raw_reads_files}", file=sys.stderr)
		sys.exit(1)
	# Upload SRA files
	with open(raw_reads_files, None) as file:
		for line in file:
			line = line.strip()
			if line is None or line == "":
				continue
			elif os.path.isfile(line):
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_19(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, "raw_reads_location.txt")
	if os.path.isfile(raw_reads_files) is False:
		print(f"Error: Submission {submission_name} is missing raw reads file at {raw_reads_files}", file=sys.stderr)
		sys.exit(1)
	# Upload SRA files
	with open("r") as file:
		for line in file:
			line = line.strip()
			if line is None or line == "":
				continue
			elif os.path.isfile(line):
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_20(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, "raw_reads_location.txt")
	if os.path.isfile(raw_reads_files) is False:
		print(f"Error: Submission {submission_name} is missing raw reads file at {raw_reads_files}", file=sys.stderr)
		sys.exit(1)
	# Upload SRA files
	with open(raw_reads_files, ) as file:
		for line in file:
			line = line.strip()
			if line is None or line == "":
				continue
			elif os.path.isfile(line):
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_21(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, "raw_reads_location.txt")
	if os.path.isfile(raw_reads_files) is False:
		print(f"Error: Submission {submission_name} is missing raw reads file at {raw_reads_files}", file=sys.stderr)
		sys.exit(1)
	# Upload SRA files
	with open(raw_reads_files, "XXrXX") as file:
		for line in file:
			line = line.strip()
			if line is None or line == "":
				continue
			elif os.path.isfile(line):
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_22(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, "raw_reads_location.txt")
	if os.path.isfile(raw_reads_files) is False:
		print(f"Error: Submission {submission_name} is missing raw reads file at {raw_reads_files}", file=sys.stderr)
		sys.exit(1)
	# Upload SRA files
	with open(raw_reads_files, "R") as file:
		for line in file:
			line = line.strip()
			if line is None or line == "":
				continue
			elif os.path.isfile(line):
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_23(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, "raw_reads_location.txt")
	if os.path.isfile(raw_reads_files) is False:
		print(f"Error: Submission {submission_name} is missing raw reads file at {raw_reads_files}", file=sys.stderr)
		sys.exit(1)
	# Upload SRA files
	with open(raw_reads_files, "r") as file:
		for line in file:
			line = None
			if line is None or line == "":
				continue
			elif os.path.isfile(line):
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_24(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, "raw_reads_location.txt")
	if os.path.isfile(raw_reads_files) is False:
		print(f"Error: Submission {submission_name} is missing raw reads file at {raw_reads_files}", file=sys.stderr)
		sys.exit(1)
	# Upload SRA files
	with open(raw_reads_files, "r") as file:
		for line in file:
			line = line.strip()
			if line is None and line == "":
				continue
			elif os.path.isfile(line):
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_25(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, "raw_reads_location.txt")
	if os.path.isfile(raw_reads_files) is False:
		print(f"Error: Submission {submission_name} is missing raw reads file at {raw_reads_files}", file=sys.stderr)
		sys.exit(1)
	# Upload SRA files
	with open(raw_reads_files, "r") as file:
		for line in file:
			line = line.strip()
			if line is not None or line == "":
				continue
			elif os.path.isfile(line):
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_26(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, "raw_reads_location.txt")
	if os.path.isfile(raw_reads_files) is False:
		print(f"Error: Submission {submission_name} is missing raw reads file at {raw_reads_files}", file=sys.stderr)
		sys.exit(1)
	# Upload SRA files
	with open(raw_reads_files, "r") as file:
		for line in file:
			line = line.strip()
			if line is None or line != "":
				continue
			elif os.path.isfile(line):
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_27(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, "raw_reads_location.txt")
	if os.path.isfile(raw_reads_files) is False:
		print(f"Error: Submission {submission_name} is missing raw reads file at {raw_reads_files}", file=sys.stderr)
		sys.exit(1)
	# Upload SRA files
	with open(raw_reads_files, "r") as file:
		for line in file:
			line = line.strip()
			if line is None or line == "XXXX":
				continue
			elif os.path.isfile(line):
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_28(ftp, submission_dir: str, submission_name: str):
	raw_reads_files = os.path.join(submission_dir, "raw_reads_location.txt")
	if os.path.isfile(raw_reads_files) is False:
		print(f"Error: Submission {submission_name} is missing raw reads file at {raw_reads_files}", file=sys.stderr)
		sys.exit(1)
	# Upload SRA files
	with open(raw_reads_files, "r") as file:
		for line in file:
			line = line.strip()
			if line is None or line == "":
				break
			elif os.path.isfile(line):
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_29(ftp, submission_dir: str, submission_name: str):
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
			elif os.path.isfile(None):
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_30(ftp, submission_dir: str, submission_name: str):
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
				ftp = None
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_31(ftp, submission_dir: str, submission_name: str):
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
				ftp = ftp_upload_file(ftp=None, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_32(ftp, submission_dir: str, submission_name: str):
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
				ftp = ftp_upload_file(ftp=ftp, upload_file=None)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_33(ftp, submission_dir: str, submission_name: str):
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
				ftp = ftp_upload_file(upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_34(ftp, submission_dir: str, submission_name: str):
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
				ftp = ftp_upload_file(ftp=ftp, )
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_35(ftp, submission_dir: str, submission_name: str):
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
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print(None, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_36(ftp, submission_dir: str, submission_name: str):
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
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=None)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_37(ftp, submission_dir: str, submission_name: str):
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
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print(file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_38(ftp, submission_dir: str, submission_name: str):
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
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, )
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_39(ftp, submission_dir: str, submission_name: str):
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
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " - line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_40(ftp, submission_dir: str, submission_name: str):
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
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("XXError: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: XX" + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_41(ftp, submission_dir: str, submission_name: str):
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
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("error: uploading files to sra database failed. possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_42(ftp, submission_dir: str, submission_name: str):
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
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("ERROR: UPLOADING FILES TO SRA DATABASE FAILED. POSSIBLY FILES HAVE BEEN MOVED OR THIS IS NOT A VALID FILE: " + line, file=sys.stderr)
				sys.exit(1)
	return ftp

def x_upload_raw_reads__mutmut_43(ftp, submission_dir: str, submission_name: str):
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
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(None)
	return ftp

def x_upload_raw_reads__mutmut_44(ftp, submission_dir: str, submission_name: str):
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
				ftp = ftp_upload_file(ftp=ftp, upload_file=line)
			else:
				print("Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: " + line, file=sys.stderr)
				sys.exit(2)
	return ftp

x_upload_raw_reads__mutmut_mutants : ClassVar[MutantDict] = { # type: ignore
'x_upload_raw_reads__mutmut_1': x_upload_raw_reads__mutmut_1, 
    'x_upload_raw_reads__mutmut_2': x_upload_raw_reads__mutmut_2, 
    'x_upload_raw_reads__mutmut_3': x_upload_raw_reads__mutmut_3, 
    'x_upload_raw_reads__mutmut_4': x_upload_raw_reads__mutmut_4, 
    'x_upload_raw_reads__mutmut_5': x_upload_raw_reads__mutmut_5, 
    'x_upload_raw_reads__mutmut_6': x_upload_raw_reads__mutmut_6, 
    'x_upload_raw_reads__mutmut_7': x_upload_raw_reads__mutmut_7, 
    'x_upload_raw_reads__mutmut_8': x_upload_raw_reads__mutmut_8, 
    'x_upload_raw_reads__mutmut_9': x_upload_raw_reads__mutmut_9, 
    'x_upload_raw_reads__mutmut_10': x_upload_raw_reads__mutmut_10, 
    'x_upload_raw_reads__mutmut_11': x_upload_raw_reads__mutmut_11, 
    'x_upload_raw_reads__mutmut_12': x_upload_raw_reads__mutmut_12, 
    'x_upload_raw_reads__mutmut_13': x_upload_raw_reads__mutmut_13, 
    'x_upload_raw_reads__mutmut_14': x_upload_raw_reads__mutmut_14, 
    'x_upload_raw_reads__mutmut_15': x_upload_raw_reads__mutmut_15, 
    'x_upload_raw_reads__mutmut_16': x_upload_raw_reads__mutmut_16, 
    'x_upload_raw_reads__mutmut_17': x_upload_raw_reads__mutmut_17, 
    'x_upload_raw_reads__mutmut_18': x_upload_raw_reads__mutmut_18, 
    'x_upload_raw_reads__mutmut_19': x_upload_raw_reads__mutmut_19, 
    'x_upload_raw_reads__mutmut_20': x_upload_raw_reads__mutmut_20, 
    'x_upload_raw_reads__mutmut_21': x_upload_raw_reads__mutmut_21, 
    'x_upload_raw_reads__mutmut_22': x_upload_raw_reads__mutmut_22, 
    'x_upload_raw_reads__mutmut_23': x_upload_raw_reads__mutmut_23, 
    'x_upload_raw_reads__mutmut_24': x_upload_raw_reads__mutmut_24, 
    'x_upload_raw_reads__mutmut_25': x_upload_raw_reads__mutmut_25, 
    'x_upload_raw_reads__mutmut_26': x_upload_raw_reads__mutmut_26, 
    'x_upload_raw_reads__mutmut_27': x_upload_raw_reads__mutmut_27, 
    'x_upload_raw_reads__mutmut_28': x_upload_raw_reads__mutmut_28, 
    'x_upload_raw_reads__mutmut_29': x_upload_raw_reads__mutmut_29, 
    'x_upload_raw_reads__mutmut_30': x_upload_raw_reads__mutmut_30, 
    'x_upload_raw_reads__mutmut_31': x_upload_raw_reads__mutmut_31, 
    'x_upload_raw_reads__mutmut_32': x_upload_raw_reads__mutmut_32, 
    'x_upload_raw_reads__mutmut_33': x_upload_raw_reads__mutmut_33, 
    'x_upload_raw_reads__mutmut_34': x_upload_raw_reads__mutmut_34, 
    'x_upload_raw_reads__mutmut_35': x_upload_raw_reads__mutmut_35, 
    'x_upload_raw_reads__mutmut_36': x_upload_raw_reads__mutmut_36, 
    'x_upload_raw_reads__mutmut_37': x_upload_raw_reads__mutmut_37, 
    'x_upload_raw_reads__mutmut_38': x_upload_raw_reads__mutmut_38, 
    'x_upload_raw_reads__mutmut_39': x_upload_raw_reads__mutmut_39, 
    'x_upload_raw_reads__mutmut_40': x_upload_raw_reads__mutmut_40, 
    'x_upload_raw_reads__mutmut_41': x_upload_raw_reads__mutmut_41, 
    'x_upload_raw_reads__mutmut_42': x_upload_raw_reads__mutmut_42, 
    'x_upload_raw_reads__mutmut_43': x_upload_raw_reads__mutmut_43, 
    'x_upload_raw_reads__mutmut_44': x_upload_raw_reads__mutmut_44
}
x_upload_raw_reads__mutmut_orig.__name__ = 'x_upload_raw_reads'

# Submit to NCBI
def submit_ncbi(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	args = [database, submission_name, submission_dir, config_dict, submission_type]# type: ignore
	kwargs = {}# type: ignore
	return _mutmut_trampoline(x_submit_ncbi__mutmut_orig, x_submit_ncbi__mutmut_mutants, args, kwargs, None)

# Submit to NCBI
def x_submit_ncbi__mutmut_orig(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_1(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = None
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_2(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" - database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_3(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name - "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_4(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "XX_XX" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_5(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=None, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_6(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database=None)
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_7(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_8(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, )
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_9(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="XXNCBIXX")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_10(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="ncbi")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_11(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(None)
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_12(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(None)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_13(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(6)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_14(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = None
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_15(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(None)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_16(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(None)
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_17(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(None)
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_18(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = None
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_19(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=None, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_20(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=None, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_21(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=None, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_22(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=None)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_23(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_24(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_25(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_26(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, )
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_27(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=False)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_28(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(None)
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_29(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = None
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_30(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=None, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_31(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=None)
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_32(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_33(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, )
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_34(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(None, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_35(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, None))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_36(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join("submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_37(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, ))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_38(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "XXsubmission.xmlXX"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_39(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "SUBMISSION.XML"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_40(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "XXSRAXX" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_41(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "sra" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_42(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" not in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_43(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = None
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_44(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=None, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_45(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=None, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_46(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=None)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_47(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_48(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_49(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, )
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_50(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "XXGENBANKXX" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_51(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "genbank" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_52(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" not in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_53(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = None
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_54(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=None, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_55(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=None)
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_56(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_57(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, )
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_58(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(None, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_59(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, None))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_60(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_61(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, ))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_62(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = None
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_63(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=None, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_64(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=None)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_65(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_66(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, )
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_67(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print(None, file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_68(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=None)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_69(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print(file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_70(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), )
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_71(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' - str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_72(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" - 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_73(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("XX\nXX" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_74(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'XXError: XX' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_75(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'error: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_76(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'ERROR: ' + str(e), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_77(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(None), file=sys.stderr)
		sys.exit(1)

# Submit to NCBI
def x_submit_ncbi__mutmut_78(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(None)

# Submit to NCBI
def x_submit_ncbi__mutmut_79(database: str, submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> None:
	# Create submission name
	ncbi_submission_name = submission_name + "_" + database
	# Check user credentials
	tools.check_credentials(config_dict=config_dict, database="NCBI")
	# Submit sequences to NCBI via FTP Server
	print(f"Uploading sample files to NCBI-{database}, as a '{submission_type}' submission. If this is not intended, interrupt immediately.")
	time.sleep(5)
	try:
		# Login into NCBI FTP Server
		ftp = ncbi_login(config_dict)
		print(f"Connecting to NCBI FTP Server")
		print(f"Submission name: {ncbi_submission_name}")
		ftp = ftp_navigate_to_folder(ftp=ftp, folder_name=ncbi_submission_name, submission_type=submission_type, make_folder=True)
		print(f"Submitting '{submission_name}'")
		# Upload submission xml
		ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, "submission.xml"))
		# Upload raw reads
		if "SRA" in database:
			ftp = upload_raw_reads(ftp=ftp, submission_dir=submission_dir, submission_name=submission_name)
		# Upload zipfile
		elif "GENBANK" in database:
			ftp = ftp_upload_file(ftp=ftp, upload_file=os.path.join(submission_dir, f"{submission_name}.zip"))
		ftp = create_submit_ready_file(ftp=ftp, submission_dir=submission_dir)
	except ftplib.all_errors as e:
		print("\n" + 'Error: ' + str(e), file=sys.stderr)
		sys.exit(2)

x_submit_ncbi__mutmut_mutants : ClassVar[MutantDict] = { # type: ignore
'x_submit_ncbi__mutmut_1': x_submit_ncbi__mutmut_1, 
    'x_submit_ncbi__mutmut_2': x_submit_ncbi__mutmut_2, 
    'x_submit_ncbi__mutmut_3': x_submit_ncbi__mutmut_3, 
    'x_submit_ncbi__mutmut_4': x_submit_ncbi__mutmut_4, 
    'x_submit_ncbi__mutmut_5': x_submit_ncbi__mutmut_5, 
    'x_submit_ncbi__mutmut_6': x_submit_ncbi__mutmut_6, 
    'x_submit_ncbi__mutmut_7': x_submit_ncbi__mutmut_7, 
    'x_submit_ncbi__mutmut_8': x_submit_ncbi__mutmut_8, 
    'x_submit_ncbi__mutmut_9': x_submit_ncbi__mutmut_9, 
    'x_submit_ncbi__mutmut_10': x_submit_ncbi__mutmut_10, 
    'x_submit_ncbi__mutmut_11': x_submit_ncbi__mutmut_11, 
    'x_submit_ncbi__mutmut_12': x_submit_ncbi__mutmut_12, 
    'x_submit_ncbi__mutmut_13': x_submit_ncbi__mutmut_13, 
    'x_submit_ncbi__mutmut_14': x_submit_ncbi__mutmut_14, 
    'x_submit_ncbi__mutmut_15': x_submit_ncbi__mutmut_15, 
    'x_submit_ncbi__mutmut_16': x_submit_ncbi__mutmut_16, 
    'x_submit_ncbi__mutmut_17': x_submit_ncbi__mutmut_17, 
    'x_submit_ncbi__mutmut_18': x_submit_ncbi__mutmut_18, 
    'x_submit_ncbi__mutmut_19': x_submit_ncbi__mutmut_19, 
    'x_submit_ncbi__mutmut_20': x_submit_ncbi__mutmut_20, 
    'x_submit_ncbi__mutmut_21': x_submit_ncbi__mutmut_21, 
    'x_submit_ncbi__mutmut_22': x_submit_ncbi__mutmut_22, 
    'x_submit_ncbi__mutmut_23': x_submit_ncbi__mutmut_23, 
    'x_submit_ncbi__mutmut_24': x_submit_ncbi__mutmut_24, 
    'x_submit_ncbi__mutmut_25': x_submit_ncbi__mutmut_25, 
    'x_submit_ncbi__mutmut_26': x_submit_ncbi__mutmut_26, 
    'x_submit_ncbi__mutmut_27': x_submit_ncbi__mutmut_27, 
    'x_submit_ncbi__mutmut_28': x_submit_ncbi__mutmut_28, 
    'x_submit_ncbi__mutmut_29': x_submit_ncbi__mutmut_29, 
    'x_submit_ncbi__mutmut_30': x_submit_ncbi__mutmut_30, 
    'x_submit_ncbi__mutmut_31': x_submit_ncbi__mutmut_31, 
    'x_submit_ncbi__mutmut_32': x_submit_ncbi__mutmut_32, 
    'x_submit_ncbi__mutmut_33': x_submit_ncbi__mutmut_33, 
    'x_submit_ncbi__mutmut_34': x_submit_ncbi__mutmut_34, 
    'x_submit_ncbi__mutmut_35': x_submit_ncbi__mutmut_35, 
    'x_submit_ncbi__mutmut_36': x_submit_ncbi__mutmut_36, 
    'x_submit_ncbi__mutmut_37': x_submit_ncbi__mutmut_37, 
    'x_submit_ncbi__mutmut_38': x_submit_ncbi__mutmut_38, 
    'x_submit_ncbi__mutmut_39': x_submit_ncbi__mutmut_39, 
    'x_submit_ncbi__mutmut_40': x_submit_ncbi__mutmut_40, 
    'x_submit_ncbi__mutmut_41': x_submit_ncbi__mutmut_41, 
    'x_submit_ncbi__mutmut_42': x_submit_ncbi__mutmut_42, 
    'x_submit_ncbi__mutmut_43': x_submit_ncbi__mutmut_43, 
    'x_submit_ncbi__mutmut_44': x_submit_ncbi__mutmut_44, 
    'x_submit_ncbi__mutmut_45': x_submit_ncbi__mutmut_45, 
    'x_submit_ncbi__mutmut_46': x_submit_ncbi__mutmut_46, 
    'x_submit_ncbi__mutmut_47': x_submit_ncbi__mutmut_47, 
    'x_submit_ncbi__mutmut_48': x_submit_ncbi__mutmut_48, 
    'x_submit_ncbi__mutmut_49': x_submit_ncbi__mutmut_49, 
    'x_submit_ncbi__mutmut_50': x_submit_ncbi__mutmut_50, 
    'x_submit_ncbi__mutmut_51': x_submit_ncbi__mutmut_51, 
    'x_submit_ncbi__mutmut_52': x_submit_ncbi__mutmut_52, 
    'x_submit_ncbi__mutmut_53': x_submit_ncbi__mutmut_53, 
    'x_submit_ncbi__mutmut_54': x_submit_ncbi__mutmut_54, 
    'x_submit_ncbi__mutmut_55': x_submit_ncbi__mutmut_55, 
    'x_submit_ncbi__mutmut_56': x_submit_ncbi__mutmut_56, 
    'x_submit_ncbi__mutmut_57': x_submit_ncbi__mutmut_57, 
    'x_submit_ncbi__mutmut_58': x_submit_ncbi__mutmut_58, 
    'x_submit_ncbi__mutmut_59': x_submit_ncbi__mutmut_59, 
    'x_submit_ncbi__mutmut_60': x_submit_ncbi__mutmut_60, 
    'x_submit_ncbi__mutmut_61': x_submit_ncbi__mutmut_61, 
    'x_submit_ncbi__mutmut_62': x_submit_ncbi__mutmut_62, 
    'x_submit_ncbi__mutmut_63': x_submit_ncbi__mutmut_63, 
    'x_submit_ncbi__mutmut_64': x_submit_ncbi__mutmut_64, 
    'x_submit_ncbi__mutmut_65': x_submit_ncbi__mutmut_65, 
    'x_submit_ncbi__mutmut_66': x_submit_ncbi__mutmut_66, 
    'x_submit_ncbi__mutmut_67': x_submit_ncbi__mutmut_67, 
    'x_submit_ncbi__mutmut_68': x_submit_ncbi__mutmut_68, 
    'x_submit_ncbi__mutmut_69': x_submit_ncbi__mutmut_69, 
    'x_submit_ncbi__mutmut_70': x_submit_ncbi__mutmut_70, 
    'x_submit_ncbi__mutmut_71': x_submit_ncbi__mutmut_71, 
    'x_submit_ncbi__mutmut_72': x_submit_ncbi__mutmut_72, 
    'x_submit_ncbi__mutmut_73': x_submit_ncbi__mutmut_73, 
    'x_submit_ncbi__mutmut_74': x_submit_ncbi__mutmut_74, 
    'x_submit_ncbi__mutmut_75': x_submit_ncbi__mutmut_75, 
    'x_submit_ncbi__mutmut_76': x_submit_ncbi__mutmut_76, 
    'x_submit_ncbi__mutmut_77': x_submit_ncbi__mutmut_77, 
    'x_submit_ncbi__mutmut_78': x_submit_ncbi__mutmut_78, 
    'x_submit_ncbi__mutmut_79': x_submit_ncbi__mutmut_79
}
x_submit_ncbi__mutmut_orig.__name__ = 'x_submit_ncbi'

# Send table2asn file through email
def email_table2asn(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	args = [submission_name, submission_dir, config_dict, submission_type]# type: ignore
	kwargs = {}# type: ignore
	return _mutmut_trampoline(x_email_table2asn__mutmut_orig, x_email_table2asn__mutmut_mutants, args, kwargs, None)

# Send table2asn file through email
def x_email_table2asn__mutmut_orig(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_1(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = None
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_2(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(None, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_3(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, None)
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_4(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_5(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, )
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_6(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name - ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_7(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + "XX.sqnXX")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_8(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".SQN")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_9(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = None
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_10(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart(None)
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_11(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('XXmultipartXX')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_12(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('MULTIPART')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_13(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = None
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_14(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['XXSubjectXX'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_15(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_16(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['SUBJECT'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_17(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name - " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_18(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + "XX table2asn submissionXX"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_19(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " TABLE2ASN SUBMISSION"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_20(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = None
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_21(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["XXDescriptionXX"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_22(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_23(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["DESCRIPTION"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_24(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["XXOrganizationXX"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_25(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_26(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["ORGANIZATION"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_27(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["XXSubmitterXX"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_28(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_29(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["SUBMITTER"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_30(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["XXEmailXX"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_31(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_32(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["EMAIL"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_33(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = None
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_34(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = None
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_35(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type != "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_36(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "XXTESTXX":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_37(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "test":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_38(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(None)
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_39(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["XXDescriptionXX"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_40(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_41(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["DESCRIPTION"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_42(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["XXOrganizationXX"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_43(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_44(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["ORGANIZATION"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_45(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["XXSubmitterXX"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_46(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_47(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["SUBMITTER"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_48(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["XXEmailXX"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_49(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_50(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["EMAIL"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_51(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(None)
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_52(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['XXDescriptionXX']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_53(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_54(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['DESCRIPTION']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_55(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['XXOrganizationXX']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_56(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_57(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['ORGANIZATION']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_58(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['XXSubmitterXX']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_59(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_60(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['SUBMITTER']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_61(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['XXEmailXX']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_62(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_63(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['EMAIL']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_64(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type != "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_65(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "XXPRODUCTIONXX":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_66(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "production":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_67(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(None)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_68(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(None)
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_69(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["XXDescriptionXX"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_70(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_71(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["DESCRIPTION"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_72(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["XXOrganizationXX"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_73(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_74(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["ORGANIZATION"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_75(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["XXSubmitterXX"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_76(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_77(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["SUBMITTER"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_78(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["XXEmailXX"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_79(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_80(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["EMAIL"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_81(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(None)
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

# Send table2asn file through email
def x_email_table2asn__mutmut_82(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(None, file=sys.stderr)
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

# Send table2asn file through email
def x_email_table2asn__mutmut_83(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=None)
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

# Send table2asn file through email
def x_email_table2asn__mutmut_84(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(file=sys.stderr)
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

# Send table2asn file through email
def x_email_table2asn__mutmut_85(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", )
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

# Send table2asn file through email
def x_email_table2asn__mutmut_86(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(None)
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

# Send table2asn file through email
def x_email_table2asn__mutmut_87(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(2)
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

# Send table2asn file through email
def x_email_table2asn__mutmut_88(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(None)
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

# Send table2asn file through email
def x_email_table2asn__mutmut_89(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(6)
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

# Send table2asn file through email
def x_email_table2asn__mutmut_90(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["XXDescriptionXX"]["Organization"]["Submitter"]["Alt_Email"]:
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

# Send table2asn file through email
def x_email_table2asn__mutmut_91(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["description"]["Organization"]["Submitter"]["Alt_Email"]:
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

# Send table2asn file through email
def x_email_table2asn__mutmut_92(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["DESCRIPTION"]["Organization"]["Submitter"]["Alt_Email"]:
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

# Send table2asn file through email
def x_email_table2asn__mutmut_93(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["XXOrganizationXX"]["Submitter"]["Alt_Email"]:
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

# Send table2asn file through email
def x_email_table2asn__mutmut_94(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["organization"]["Submitter"]["Alt_Email"]:
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

# Send table2asn file through email
def x_email_table2asn__mutmut_95(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["ORGANIZATION"]["Submitter"]["Alt_Email"]:
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

# Send table2asn file through email
def x_email_table2asn__mutmut_96(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["XXSubmitterXX"]["Alt_Email"]:
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

# Send table2asn file through email
def x_email_table2asn__mutmut_97(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["submitter"]["Alt_Email"]:
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

# Send table2asn file through email
def x_email_table2asn__mutmut_98(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["SUBMITTER"]["Alt_Email"]:
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

# Send table2asn file through email
def x_email_table2asn__mutmut_99(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["XXAlt_EmailXX"]:
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

# Send table2asn file through email
def x_email_table2asn__mutmut_100(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["alt_email"]:
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

# Send table2asn file through email
def x_email_table2asn__mutmut_101(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["ALT_EMAIL"]:
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

# Send table2asn file through email
def x_email_table2asn__mutmut_102(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(None)
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

# Send table2asn file through email
def x_email_table2asn__mutmut_103(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["XXDescriptionXX"]["Organization"]["Submitter"]["Alt_Email"])
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

# Send table2asn file through email
def x_email_table2asn__mutmut_104(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["description"]["Organization"]["Submitter"]["Alt_Email"])
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

# Send table2asn file through email
def x_email_table2asn__mutmut_105(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["DESCRIPTION"]["Organization"]["Submitter"]["Alt_Email"])
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

# Send table2asn file through email
def x_email_table2asn__mutmut_106(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["XXOrganizationXX"]["Submitter"]["Alt_Email"])
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

# Send table2asn file through email
def x_email_table2asn__mutmut_107(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["organization"]["Submitter"]["Alt_Email"])
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

# Send table2asn file through email
def x_email_table2asn__mutmut_108(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["ORGANIZATION"]["Submitter"]["Alt_Email"])
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

# Send table2asn file through email
def x_email_table2asn__mutmut_109(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["XXSubmitterXX"]["Alt_Email"])
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

# Send table2asn file through email
def x_email_table2asn__mutmut_110(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["submitter"]["Alt_Email"])
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

# Send table2asn file through email
def x_email_table2asn__mutmut_111(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["SUBMITTER"]["Alt_Email"])
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

# Send table2asn file through email
def x_email_table2asn__mutmut_112(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["XXAlt_EmailXX"])
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

# Send table2asn file through email
def x_email_table2asn__mutmut_113(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["alt_email"])
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

# Send table2asn file through email
def x_email_table2asn__mutmut_114(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["ALT_EMAIL"])
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

# Send table2asn file through email
def x_email_table2asn__mutmut_115(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"])
		msg['From'] = None
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

# Send table2asn file through email
def x_email_table2asn__mutmut_116(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"])
		msg['XXFromXX'] = from_email
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

# Send table2asn file through email
def x_email_table2asn__mutmut_117(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"])
		msg['from'] = from_email
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

# Send table2asn file through email
def x_email_table2asn__mutmut_118(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"])
		msg['FROM'] = from_email
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

# Send table2asn file through email
def x_email_table2asn__mutmut_119(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"])
		msg['From'] = from_email
		msg['To'] = None
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

# Send table2asn file through email
def x_email_table2asn__mutmut_120(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"])
		msg['From'] = from_email
		msg['XXToXX'] = ", ".join(to_email)
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

# Send table2asn file through email
def x_email_table2asn__mutmut_121(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"])
		msg['From'] = from_email
		msg['to'] = ", ".join(to_email)
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

# Send table2asn file through email
def x_email_table2asn__mutmut_122(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"])
		msg['From'] = from_email
		msg['TO'] = ", ".join(to_email)
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

# Send table2asn file through email
def x_email_table2asn__mutmut_123(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"])
		msg['From'] = from_email
		msg['To'] = ", ".join(None)
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

# Send table2asn file through email
def x_email_table2asn__mutmut_124(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"])
		msg['From'] = from_email
		msg['To'] = "XX, XX".join(to_email)
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

# Send table2asn file through email
def x_email_table2asn__mutmut_125(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"])
		msg['From'] = from_email
		msg['To'] = ", ".join(to_email)
		if len(cc_email) == 0:
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

# Send table2asn file through email
def x_email_table2asn__mutmut_126(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"])
		msg['From'] = from_email
		msg['To'] = ", ".join(to_email)
		if len(cc_email) != 1:
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

# Send table2asn file through email
def x_email_table2asn__mutmut_127(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"])
		msg['From'] = from_email
		msg['To'] = ", ".join(to_email)
		if len(cc_email) != 0:
			msg['Cc'] = None
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

# Send table2asn file through email
def x_email_table2asn__mutmut_128(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"])
		msg['From'] = from_email
		msg['To'] = ", ".join(to_email)
		if len(cc_email) != 0:
			msg['XXCcXX'] = ", ".join(cc_email)
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

# Send table2asn file through email
def x_email_table2asn__mutmut_129(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"])
		msg['From'] = from_email
		msg['To'] = ", ".join(to_email)
		if len(cc_email) != 0:
			msg['cc'] = ", ".join(cc_email)
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

# Send table2asn file through email
def x_email_table2asn__mutmut_130(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"])
		msg['From'] = from_email
		msg['To'] = ", ".join(to_email)
		if len(cc_email) != 0:
			msg['CC'] = ", ".join(cc_email)
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

# Send table2asn file through email
def x_email_table2asn__mutmut_131(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"])
		msg['From'] = from_email
		msg['To'] = ", ".join(to_email)
		if len(cc_email) != 0:
			msg['Cc'] = ", ".join(None)
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

# Send table2asn file through email
def x_email_table2asn__mutmut_132(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
		else:
			print(f"Error: Submission type '{submission_type}' is not a valid option.", file=sys.stderr)
			sys.exit(1)
		time.sleep(5)
		if config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"]:
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Alt_Email"])
		msg['From'] = from_email
		msg['To'] = ", ".join(to_email)
		if len(cc_email) != 0:
			msg['Cc'] = "XX, XX".join(cc_email)
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

# Send table2asn file through email
def x_email_table2asn__mutmut_133(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		with open(None, 'rb') as file_input:
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

# Send table2asn file through email
def x_email_table2asn__mutmut_134(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		with open(sqn_file, None) as file_input:
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

# Send table2asn file through email
def x_email_table2asn__mutmut_135(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		with open('rb') as file_input:
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

# Send table2asn file through email
def x_email_table2asn__mutmut_136(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		with open(sqn_file, ) as file_input:
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

# Send table2asn file through email
def x_email_table2asn__mutmut_137(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		with open(sqn_file, 'XXrbXX') as file_input:
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

# Send table2asn file through email
def x_email_table2asn__mutmut_138(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		with open(sqn_file, 'RB') as file_input:
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

# Send table2asn file through email
def x_email_table2asn__mutmut_139(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
			part = None
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

# Send table2asn file through email
def x_email_table2asn__mutmut_140(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
			part = MIMEApplication(None, Name=submission_name + ".sqn")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_141(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
			part = MIMEApplication(file_input.read(), Name=None)
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

# Send table2asn file through email
def x_email_table2asn__mutmut_142(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
			part = MIMEApplication(Name=submission_name + ".sqn")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_143(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
			part = MIMEApplication(file_input.read(), )
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

# Send table2asn file through email
def x_email_table2asn__mutmut_144(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
			part = MIMEApplication(file_input.read(), Name=submission_name - ".sqn")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_145(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
			part = MIMEApplication(file_input.read(), Name=submission_name + "XX.sqnXX")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_146(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
			part = MIMEApplication(file_input.read(), Name=submission_name + ".SQN")
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

# Send table2asn file through email
def x_email_table2asn__mutmut_147(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		part['Content-Disposition'] = None
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

# Send table2asn file through email
def x_email_table2asn__mutmut_148(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		part['XXContent-DispositionXX'] = "attachment; filename=" + submission_name + ".sqn"
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

# Send table2asn file through email
def x_email_table2asn__mutmut_149(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		part['content-disposition'] = "attachment; filename=" + submission_name + ".sqn"
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

# Send table2asn file through email
def x_email_table2asn__mutmut_150(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		part['CONTENT-DISPOSITION'] = "attachment; filename=" + submission_name + ".sqn"
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

# Send table2asn file through email
def x_email_table2asn__mutmut_151(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		part['Content-Disposition'] = "attachment; filename=" + submission_name - ".sqn"
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

# Send table2asn file through email
def x_email_table2asn__mutmut_152(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		part['Content-Disposition'] = "attachment; filename=" - submission_name + ".sqn"
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

# Send table2asn file through email
def x_email_table2asn__mutmut_153(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		part['Content-Disposition'] = "XXattachment; filename=XX" + submission_name + ".sqn"
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

# Send table2asn file through email
def x_email_table2asn__mutmut_154(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		part['Content-Disposition'] = "ATTACHMENT; FILENAME=" + submission_name + ".sqn"
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

# Send table2asn file through email
def x_email_table2asn__mutmut_155(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		part['Content-Disposition'] = "attachment; filename=" + submission_name + "XX.sqnXX"
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

# Send table2asn file through email
def x_email_table2asn__mutmut_156(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		part['Content-Disposition'] = "attachment; filename=" + submission_name + ".SQN"
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

# Send table2asn file through email
def x_email_table2asn__mutmut_157(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		msg.attach(None)
		s = smtplib.SMTP('localhost')
		s.sendmail(from_email, to_email, msg.as_string())
		submission_status = "PROCESSED"
	except Exception as e:
		print("Error: Unable to send mail automatically. If unable to email, submission can be made manually using the sqn file.", file=sys.stderr)
		print(f"sqn_file:{sqn_file}", file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_158(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		s = None
		s.sendmail(from_email, to_email, msg.as_string())
		submission_status = "PROCESSED"
	except Exception as e:
		print("Error: Unable to send mail automatically. If unable to email, submission can be made manually using the sqn file.", file=sys.stderr)
		print(f"sqn_file:{sqn_file}", file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_159(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		s = smtplib.SMTP(None)
		s.sendmail(from_email, to_email, msg.as_string())
		submission_status = "PROCESSED"
	except Exception as e:
		print("Error: Unable to send mail automatically. If unable to email, submission can be made manually using the sqn file.", file=sys.stderr)
		print(f"sqn_file:{sqn_file}", file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_160(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		s = smtplib.SMTP('XXlocalhostXX')
		s.sendmail(from_email, to_email, msg.as_string())
		submission_status = "PROCESSED"
	except Exception as e:
		print("Error: Unable to send mail automatically. If unable to email, submission can be made manually using the sqn file.", file=sys.stderr)
		print(f"sqn_file:{sqn_file}", file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_161(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		s = smtplib.SMTP('LOCALHOST')
		s.sendmail(from_email, to_email, msg.as_string())
		submission_status = "PROCESSED"
	except Exception as e:
		print("Error: Unable to send mail automatically. If unable to email, submission can be made manually using the sqn file.", file=sys.stderr)
		print(f"sqn_file:{sqn_file}", file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_162(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		s.sendmail(None, to_email, msg.as_string())
		submission_status = "PROCESSED"
	except Exception as e:
		print("Error: Unable to send mail automatically. If unable to email, submission can be made manually using the sqn file.", file=sys.stderr)
		print(f"sqn_file:{sqn_file}", file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_163(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		s.sendmail(from_email, None, msg.as_string())
		submission_status = "PROCESSED"
	except Exception as e:
		print("Error: Unable to send mail automatically. If unable to email, submission can be made manually using the sqn file.", file=sys.stderr)
		print(f"sqn_file:{sqn_file}", file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_164(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		s.sendmail(from_email, to_email, None)
		submission_status = "PROCESSED"
	except Exception as e:
		print("Error: Unable to send mail automatically. If unable to email, submission can be made manually using the sqn file.", file=sys.stderr)
		print(f"sqn_file:{sqn_file}", file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_165(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		s.sendmail(to_email, msg.as_string())
		submission_status = "PROCESSED"
	except Exception as e:
		print("Error: Unable to send mail automatically. If unable to email, submission can be made manually using the sqn file.", file=sys.stderr)
		print(f"sqn_file:{sqn_file}", file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_166(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		s.sendmail(from_email, msg.as_string())
		submission_status = "PROCESSED"
	except Exception as e:
		print("Error: Unable to send mail automatically. If unable to email, submission can be made manually using the sqn file.", file=sys.stderr)
		print(f"sqn_file:{sqn_file}", file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_167(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		s.sendmail(from_email, to_email, )
		submission_status = "PROCESSED"
	except Exception as e:
		print("Error: Unable to send mail automatically. If unable to email, submission can be made manually using the sqn file.", file=sys.stderr)
		print(f"sqn_file:{sqn_file}", file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_168(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		submission_status = None
	except Exception as e:
		print("Error: Unable to send mail automatically. If unable to email, submission can be made manually using the sqn file.", file=sys.stderr)
		print(f"sqn_file:{sqn_file}", file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_169(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		submission_status = "XXPROCESSEDXX"
	except Exception as e:
		print("Error: Unable to send mail automatically. If unable to email, submission can be made manually using the sqn file.", file=sys.stderr)
		print(f"sqn_file:{sqn_file}", file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_170(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		submission_status = "processed"
	except Exception as e:
		print("Error: Unable to send mail automatically. If unable to email, submission can be made manually using the sqn file.", file=sys.stderr)
		print(f"sqn_file:{sqn_file}", file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_171(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		print(None, file=sys.stderr)
		print(f"sqn_file:{sqn_file}", file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_172(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		print("Error: Unable to send mail automatically. If unable to email, submission can be made manually using the sqn file.", file=None)
		print(f"sqn_file:{sqn_file}", file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_173(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		print(file=sys.stderr)
		print(f"sqn_file:{sqn_file}", file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_174(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		print("Error: Unable to send mail automatically. If unable to email, submission can be made manually using the sqn file.", )
		print(f"sqn_file:{sqn_file}", file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_175(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		print("XXError: Unable to send mail automatically. If unable to email, submission can be made manually using the sqn file.XX", file=sys.stderr)
		print(f"sqn_file:{sqn_file}", file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_176(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		print("error: unable to send mail automatically. if unable to email, submission can be made manually using the sqn file.", file=sys.stderr)
		print(f"sqn_file:{sqn_file}", file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_177(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		print("ERROR: UNABLE TO SEND MAIL AUTOMATICALLY. IF UNABLE TO EMAIL, SUBMISSION CAN BE MADE MANUALLY USING THE SQN FILE.", file=sys.stderr)
		print(f"sqn_file:{sqn_file}", file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_178(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		print(None, file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_179(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		print(f"sqn_file:{sqn_file}", file=None)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_180(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		print(file=sys.stderr)
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_181(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		print(f"sqn_file:{sqn_file}", )
		print(e, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_182(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		print(None, file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_183(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		print(e, file=None)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_184(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		print(file=sys.stderr)
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_185(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		print(e, )
		submission_status = "ERROR"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_186(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		submission_status = None
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_187(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		submission_status = "XXERRORXX"
	return submission_status

# Send table2asn file through email
def x_email_table2asn__mutmut_188(submission_name: str, submission_dir: str, config_dict: dict[str, Any], submission_type: str) -> str:
	sqn_file = os.path.join(submission_dir, submission_name + ".sqn")
	try:
		msg = MIMEMultipart('multipart')
		msg['Subject'] = submission_name + " table2asn submission"
		from_email = config_dict["Description"]["Organization"]["Submitter"]["Email"]
		to_email = []
		cc_email = []
		if submission_type == "TEST":
			to_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to submitter '{config_dict['Description']['Organization']['Submitter']['Email']}' as a 'TEST' submission. If this is not intended, interrupt immediately.")
		elif submission_type == "PRODUCTION":
			to_email.append(TABLE2ASN_EMAIL)
			cc_email.append(config_dict["Description"]["Organization"]["Submitter"]["Email"])
			print(f"Emailing table2asn sqn file to NCBI-GENBANK '{TABLE2ASN_EMAIL}', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.")
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
		submission_status = "error"
	return submission_status

x_email_table2asn__mutmut_mutants : ClassVar[MutantDict] = { # type: ignore
'x_email_table2asn__mutmut_1': x_email_table2asn__mutmut_1, 
    'x_email_table2asn__mutmut_2': x_email_table2asn__mutmut_2, 
    'x_email_table2asn__mutmut_3': x_email_table2asn__mutmut_3, 
    'x_email_table2asn__mutmut_4': x_email_table2asn__mutmut_4, 
    'x_email_table2asn__mutmut_5': x_email_table2asn__mutmut_5, 
    'x_email_table2asn__mutmut_6': x_email_table2asn__mutmut_6, 
    'x_email_table2asn__mutmut_7': x_email_table2asn__mutmut_7, 
    'x_email_table2asn__mutmut_8': x_email_table2asn__mutmut_8, 
    'x_email_table2asn__mutmut_9': x_email_table2asn__mutmut_9, 
    'x_email_table2asn__mutmut_10': x_email_table2asn__mutmut_10, 
    'x_email_table2asn__mutmut_11': x_email_table2asn__mutmut_11, 
    'x_email_table2asn__mutmut_12': x_email_table2asn__mutmut_12, 
    'x_email_table2asn__mutmut_13': x_email_table2asn__mutmut_13, 
    'x_email_table2asn__mutmut_14': x_email_table2asn__mutmut_14, 
    'x_email_table2asn__mutmut_15': x_email_table2asn__mutmut_15, 
    'x_email_table2asn__mutmut_16': x_email_table2asn__mutmut_16, 
    'x_email_table2asn__mutmut_17': x_email_table2asn__mutmut_17, 
    'x_email_table2asn__mutmut_18': x_email_table2asn__mutmut_18, 
    'x_email_table2asn__mutmut_19': x_email_table2asn__mutmut_19, 
    'x_email_table2asn__mutmut_20': x_email_table2asn__mutmut_20, 
    'x_email_table2asn__mutmut_21': x_email_table2asn__mutmut_21, 
    'x_email_table2asn__mutmut_22': x_email_table2asn__mutmut_22, 
    'x_email_table2asn__mutmut_23': x_email_table2asn__mutmut_23, 
    'x_email_table2asn__mutmut_24': x_email_table2asn__mutmut_24, 
    'x_email_table2asn__mutmut_25': x_email_table2asn__mutmut_25, 
    'x_email_table2asn__mutmut_26': x_email_table2asn__mutmut_26, 
    'x_email_table2asn__mutmut_27': x_email_table2asn__mutmut_27, 
    'x_email_table2asn__mutmut_28': x_email_table2asn__mutmut_28, 
    'x_email_table2asn__mutmut_29': x_email_table2asn__mutmut_29, 
    'x_email_table2asn__mutmut_30': x_email_table2asn__mutmut_30, 
    'x_email_table2asn__mutmut_31': x_email_table2asn__mutmut_31, 
    'x_email_table2asn__mutmut_32': x_email_table2asn__mutmut_32, 
    'x_email_table2asn__mutmut_33': x_email_table2asn__mutmut_33, 
    'x_email_table2asn__mutmut_34': x_email_table2asn__mutmut_34, 
    'x_email_table2asn__mutmut_35': x_email_table2asn__mutmut_35, 
    'x_email_table2asn__mutmut_36': x_email_table2asn__mutmut_36, 
    'x_email_table2asn__mutmut_37': x_email_table2asn__mutmut_37, 
    'x_email_table2asn__mutmut_38': x_email_table2asn__mutmut_38, 
    'x_email_table2asn__mutmut_39': x_email_table2asn__mutmut_39, 
    'x_email_table2asn__mutmut_40': x_email_table2asn__mutmut_40, 
    'x_email_table2asn__mutmut_41': x_email_table2asn__mutmut_41, 
    'x_email_table2asn__mutmut_42': x_email_table2asn__mutmut_42, 
    'x_email_table2asn__mutmut_43': x_email_table2asn__mutmut_43, 
    'x_email_table2asn__mutmut_44': x_email_table2asn__mutmut_44, 
    'x_email_table2asn__mutmut_45': x_email_table2asn__mutmut_45, 
    'x_email_table2asn__mutmut_46': x_email_table2asn__mutmut_46, 
    'x_email_table2asn__mutmut_47': x_email_table2asn__mutmut_47, 
    'x_email_table2asn__mutmut_48': x_email_table2asn__mutmut_48, 
    'x_email_table2asn__mutmut_49': x_email_table2asn__mutmut_49, 
    'x_email_table2asn__mutmut_50': x_email_table2asn__mutmut_50, 
    'x_email_table2asn__mutmut_51': x_email_table2asn__mutmut_51, 
    'x_email_table2asn__mutmut_52': x_email_table2asn__mutmut_52, 
    'x_email_table2asn__mutmut_53': x_email_table2asn__mutmut_53, 
    'x_email_table2asn__mutmut_54': x_email_table2asn__mutmut_54, 
    'x_email_table2asn__mutmut_55': x_email_table2asn__mutmut_55, 
    'x_email_table2asn__mutmut_56': x_email_table2asn__mutmut_56, 
    'x_email_table2asn__mutmut_57': x_email_table2asn__mutmut_57, 
    'x_email_table2asn__mutmut_58': x_email_table2asn__mutmut_58, 
    'x_email_table2asn__mutmut_59': x_email_table2asn__mutmut_59, 
    'x_email_table2asn__mutmut_60': x_email_table2asn__mutmut_60, 
    'x_email_table2asn__mutmut_61': x_email_table2asn__mutmut_61, 
    'x_email_table2asn__mutmut_62': x_email_table2asn__mutmut_62, 
    'x_email_table2asn__mutmut_63': x_email_table2asn__mutmut_63, 
    'x_email_table2asn__mutmut_64': x_email_table2asn__mutmut_64, 
    'x_email_table2asn__mutmut_65': x_email_table2asn__mutmut_65, 
    'x_email_table2asn__mutmut_66': x_email_table2asn__mutmut_66, 
    'x_email_table2asn__mutmut_67': x_email_table2asn__mutmut_67, 
    'x_email_table2asn__mutmut_68': x_email_table2asn__mutmut_68, 
    'x_email_table2asn__mutmut_69': x_email_table2asn__mutmut_69, 
    'x_email_table2asn__mutmut_70': x_email_table2asn__mutmut_70, 
    'x_email_table2asn__mutmut_71': x_email_table2asn__mutmut_71, 
    'x_email_table2asn__mutmut_72': x_email_table2asn__mutmut_72, 
    'x_email_table2asn__mutmut_73': x_email_table2asn__mutmut_73, 
    'x_email_table2asn__mutmut_74': x_email_table2asn__mutmut_74, 
    'x_email_table2asn__mutmut_75': x_email_table2asn__mutmut_75, 
    'x_email_table2asn__mutmut_76': x_email_table2asn__mutmut_76, 
    'x_email_table2asn__mutmut_77': x_email_table2asn__mutmut_77, 
    'x_email_table2asn__mutmut_78': x_email_table2asn__mutmut_78, 
    'x_email_table2asn__mutmut_79': x_email_table2asn__mutmut_79, 
    'x_email_table2asn__mutmut_80': x_email_table2asn__mutmut_80, 
    'x_email_table2asn__mutmut_81': x_email_table2asn__mutmut_81, 
    'x_email_table2asn__mutmut_82': x_email_table2asn__mutmut_82, 
    'x_email_table2asn__mutmut_83': x_email_table2asn__mutmut_83, 
    'x_email_table2asn__mutmut_84': x_email_table2asn__mutmut_84, 
    'x_email_table2asn__mutmut_85': x_email_table2asn__mutmut_85, 
    'x_email_table2asn__mutmut_86': x_email_table2asn__mutmut_86, 
    'x_email_table2asn__mutmut_87': x_email_table2asn__mutmut_87, 
    'x_email_table2asn__mutmut_88': x_email_table2asn__mutmut_88, 
    'x_email_table2asn__mutmut_89': x_email_table2asn__mutmut_89, 
    'x_email_table2asn__mutmut_90': x_email_table2asn__mutmut_90, 
    'x_email_table2asn__mutmut_91': x_email_table2asn__mutmut_91, 
    'x_email_table2asn__mutmut_92': x_email_table2asn__mutmut_92, 
    'x_email_table2asn__mutmut_93': x_email_table2asn__mutmut_93, 
    'x_email_table2asn__mutmut_94': x_email_table2asn__mutmut_94, 
    'x_email_table2asn__mutmut_95': x_email_table2asn__mutmut_95, 
    'x_email_table2asn__mutmut_96': x_email_table2asn__mutmut_96, 
    'x_email_table2asn__mutmut_97': x_email_table2asn__mutmut_97, 
    'x_email_table2asn__mutmut_98': x_email_table2asn__mutmut_98, 
    'x_email_table2asn__mutmut_99': x_email_table2asn__mutmut_99, 
    'x_email_table2asn__mutmut_100': x_email_table2asn__mutmut_100, 
    'x_email_table2asn__mutmut_101': x_email_table2asn__mutmut_101, 
    'x_email_table2asn__mutmut_102': x_email_table2asn__mutmut_102, 
    'x_email_table2asn__mutmut_103': x_email_table2asn__mutmut_103, 
    'x_email_table2asn__mutmut_104': x_email_table2asn__mutmut_104, 
    'x_email_table2asn__mutmut_105': x_email_table2asn__mutmut_105, 
    'x_email_table2asn__mutmut_106': x_email_table2asn__mutmut_106, 
    'x_email_table2asn__mutmut_107': x_email_table2asn__mutmut_107, 
    'x_email_table2asn__mutmut_108': x_email_table2asn__mutmut_108, 
    'x_email_table2asn__mutmut_109': x_email_table2asn__mutmut_109, 
    'x_email_table2asn__mutmut_110': x_email_table2asn__mutmut_110, 
    'x_email_table2asn__mutmut_111': x_email_table2asn__mutmut_111, 
    'x_email_table2asn__mutmut_112': x_email_table2asn__mutmut_112, 
    'x_email_table2asn__mutmut_113': x_email_table2asn__mutmut_113, 
    'x_email_table2asn__mutmut_114': x_email_table2asn__mutmut_114, 
    'x_email_table2asn__mutmut_115': x_email_table2asn__mutmut_115, 
    'x_email_table2asn__mutmut_116': x_email_table2asn__mutmut_116, 
    'x_email_table2asn__mutmut_117': x_email_table2asn__mutmut_117, 
    'x_email_table2asn__mutmut_118': x_email_table2asn__mutmut_118, 
    'x_email_table2asn__mutmut_119': x_email_table2asn__mutmut_119, 
    'x_email_table2asn__mutmut_120': x_email_table2asn__mutmut_120, 
    'x_email_table2asn__mutmut_121': x_email_table2asn__mutmut_121, 
    'x_email_table2asn__mutmut_122': x_email_table2asn__mutmut_122, 
    'x_email_table2asn__mutmut_123': x_email_table2asn__mutmut_123, 
    'x_email_table2asn__mutmut_124': x_email_table2asn__mutmut_124, 
    'x_email_table2asn__mutmut_125': x_email_table2asn__mutmut_125, 
    'x_email_table2asn__mutmut_126': x_email_table2asn__mutmut_126, 
    'x_email_table2asn__mutmut_127': x_email_table2asn__mutmut_127, 
    'x_email_table2asn__mutmut_128': x_email_table2asn__mutmut_128, 
    'x_email_table2asn__mutmut_129': x_email_table2asn__mutmut_129, 
    'x_email_table2asn__mutmut_130': x_email_table2asn__mutmut_130, 
    'x_email_table2asn__mutmut_131': x_email_table2asn__mutmut_131, 
    'x_email_table2asn__mutmut_132': x_email_table2asn__mutmut_132, 
    'x_email_table2asn__mutmut_133': x_email_table2asn__mutmut_133, 
    'x_email_table2asn__mutmut_134': x_email_table2asn__mutmut_134, 
    'x_email_table2asn__mutmut_135': x_email_table2asn__mutmut_135, 
    'x_email_table2asn__mutmut_136': x_email_table2asn__mutmut_136, 
    'x_email_table2asn__mutmut_137': x_email_table2asn__mutmut_137, 
    'x_email_table2asn__mutmut_138': x_email_table2asn__mutmut_138, 
    'x_email_table2asn__mutmut_139': x_email_table2asn__mutmut_139, 
    'x_email_table2asn__mutmut_140': x_email_table2asn__mutmut_140, 
    'x_email_table2asn__mutmut_141': x_email_table2asn__mutmut_141, 
    'x_email_table2asn__mutmut_142': x_email_table2asn__mutmut_142, 
    'x_email_table2asn__mutmut_143': x_email_table2asn__mutmut_143, 
    'x_email_table2asn__mutmut_144': x_email_table2asn__mutmut_144, 
    'x_email_table2asn__mutmut_145': x_email_table2asn__mutmut_145, 
    'x_email_table2asn__mutmut_146': x_email_table2asn__mutmut_146, 
    'x_email_table2asn__mutmut_147': x_email_table2asn__mutmut_147, 
    'x_email_table2asn__mutmut_148': x_email_table2asn__mutmut_148, 
    'x_email_table2asn__mutmut_149': x_email_table2asn__mutmut_149, 
    'x_email_table2asn__mutmut_150': x_email_table2asn__mutmut_150, 
    'x_email_table2asn__mutmut_151': x_email_table2asn__mutmut_151, 
    'x_email_table2asn__mutmut_152': x_email_table2asn__mutmut_152, 
    'x_email_table2asn__mutmut_153': x_email_table2asn__mutmut_153, 
    'x_email_table2asn__mutmut_154': x_email_table2asn__mutmut_154, 
    'x_email_table2asn__mutmut_155': x_email_table2asn__mutmut_155, 
    'x_email_table2asn__mutmut_156': x_email_table2asn__mutmut_156, 
    'x_email_table2asn__mutmut_157': x_email_table2asn__mutmut_157, 
    'x_email_table2asn__mutmut_158': x_email_table2asn__mutmut_158, 
    'x_email_table2asn__mutmut_159': x_email_table2asn__mutmut_159, 
    'x_email_table2asn__mutmut_160': x_email_table2asn__mutmut_160, 
    'x_email_table2asn__mutmut_161': x_email_table2asn__mutmut_161, 
    'x_email_table2asn__mutmut_162': x_email_table2asn__mutmut_162, 
    'x_email_table2asn__mutmut_163': x_email_table2asn__mutmut_163, 
    'x_email_table2asn__mutmut_164': x_email_table2asn__mutmut_164, 
    'x_email_table2asn__mutmut_165': x_email_table2asn__mutmut_165, 
    'x_email_table2asn__mutmut_166': x_email_table2asn__mutmut_166, 
    'x_email_table2asn__mutmut_167': x_email_table2asn__mutmut_167, 
    'x_email_table2asn__mutmut_168': x_email_table2asn__mutmut_168, 
    'x_email_table2asn__mutmut_169': x_email_table2asn__mutmut_169, 
    'x_email_table2asn__mutmut_170': x_email_table2asn__mutmut_170, 
    'x_email_table2asn__mutmut_171': x_email_table2asn__mutmut_171, 
    'x_email_table2asn__mutmut_172': x_email_table2asn__mutmut_172, 
    'x_email_table2asn__mutmut_173': x_email_table2asn__mutmut_173, 
    'x_email_table2asn__mutmut_174': x_email_table2asn__mutmut_174, 
    'x_email_table2asn__mutmut_175': x_email_table2asn__mutmut_175, 
    'x_email_table2asn__mutmut_176': x_email_table2asn__mutmut_176, 
    'x_email_table2asn__mutmut_177': x_email_table2asn__mutmut_177, 
    'x_email_table2asn__mutmut_178': x_email_table2asn__mutmut_178, 
    'x_email_table2asn__mutmut_179': x_email_table2asn__mutmut_179, 
    'x_email_table2asn__mutmut_180': x_email_table2asn__mutmut_180, 
    'x_email_table2asn__mutmut_181': x_email_table2asn__mutmut_181, 
    'x_email_table2asn__mutmut_182': x_email_table2asn__mutmut_182, 
    'x_email_table2asn__mutmut_183': x_email_table2asn__mutmut_183, 
    'x_email_table2asn__mutmut_184': x_email_table2asn__mutmut_184, 
    'x_email_table2asn__mutmut_185': x_email_table2asn__mutmut_185, 
    'x_email_table2asn__mutmut_186': x_email_table2asn__mutmut_186, 
    'x_email_table2asn__mutmut_187': x_email_table2asn__mutmut_187, 
    'x_email_table2asn__mutmut_188': x_email_table2asn__mutmut_188
}
x_email_table2asn__mutmut_orig.__name__ = 'x_email_table2asn'

def standardize_submission_status(submission_status: str) -> str:
	args = [submission_status]# type: ignore
	kwargs = {}# type: ignore
	return _mutmut_trampoline(x_standardize_submission_status__mutmut_orig, x_standardize_submission_status__mutmut_mutants, args, kwargs, None)

def x_standardize_submission_status__mutmut_orig(submission_status: str) -> str:
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

def x_standardize_submission_status__mutmut_1(submission_status: str) -> str:
	submission_status = None
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

def x_standardize_submission_status__mutmut_2(submission_status: str) -> str:
	submission_status = submission_status.strip().upper()
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

def x_standardize_submission_status__mutmut_3(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "XXsubmittedXX" in submission_status:
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

def x_standardize_submission_status__mutmut_4(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "SUBMITTED" in submission_status:
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

def x_standardize_submission_status__mutmut_5(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "submitted" not in submission_status:
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

def x_standardize_submission_status__mutmut_6(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "submitted" in submission_status:
		return "XXSUBMITTEDXX"
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

def x_standardize_submission_status__mutmut_7(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "submitted" in submission_status:
		return "submitted"
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

def x_standardize_submission_status__mutmut_8(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "submitted" in submission_status:
		return "SUBMITTED"
	elif "XXcreatedXX" in submission_status:
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

def x_standardize_submission_status__mutmut_9(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "submitted" in submission_status:
		return "SUBMITTED"
	elif "CREATED" in submission_status:
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

def x_standardize_submission_status__mutmut_10(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "submitted" in submission_status:
		return "SUBMITTED"
	elif "created" not in submission_status:
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

def x_standardize_submission_status__mutmut_11(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "submitted" in submission_status:
		return "SUBMITTED"
	elif "created" in submission_status:
		return "XXCREATEDXX"
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

def x_standardize_submission_status__mutmut_12(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "submitted" in submission_status:
		return "SUBMITTED"
	elif "created" in submission_status:
		return "created"
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

def x_standardize_submission_status__mutmut_13(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "submitted" in submission_status:
		return "SUBMITTED"
	elif "created" in submission_status:
		return "CREATED"
	elif "XXqueuedXX" in submission_status:
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

def x_standardize_submission_status__mutmut_14(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "submitted" in submission_status:
		return "SUBMITTED"
	elif "created" in submission_status:
		return "CREATED"
	elif "QUEUED" in submission_status:
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

def x_standardize_submission_status__mutmut_15(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "submitted" in submission_status:
		return "SUBMITTED"
	elif "created" in submission_status:
		return "CREATED"
	elif "queued" not in submission_status:
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

def x_standardize_submission_status__mutmut_16(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "submitted" in submission_status:
		return "SUBMITTED"
	elif "created" in submission_status:
		return "CREATED"
	elif "queued" in submission_status:
		return "XXQUEUEDXX"
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

def x_standardize_submission_status__mutmut_17(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "submitted" in submission_status:
		return "SUBMITTED"
	elif "created" in submission_status:
		return "CREATED"
	elif "queued" in submission_status:
		return "queued"
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

def x_standardize_submission_status__mutmut_18(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "submitted" in submission_status:
		return "SUBMITTED"
	elif "created" in submission_status:
		return "CREATED"
	elif "queued" in submission_status:
		return "QUEUED"
	elif "XXprocessingXX" in submission_status:
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

def x_standardize_submission_status__mutmut_19(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "submitted" in submission_status:
		return "SUBMITTED"
	elif "created" in submission_status:
		return "CREATED"
	elif "queued" in submission_status:
		return "QUEUED"
	elif "PROCESSING" in submission_status:
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

def x_standardize_submission_status__mutmut_20(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "submitted" in submission_status:
		return "SUBMITTED"
	elif "created" in submission_status:
		return "CREATED"
	elif "queued" in submission_status:
		return "QUEUED"
	elif "processing" not in submission_status:
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

def x_standardize_submission_status__mutmut_21(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "submitted" in submission_status:
		return "SUBMITTED"
	elif "created" in submission_status:
		return "CREATED"
	elif "queued" in submission_status:
		return "QUEUED"
	elif "processing" in submission_status:
		return "XXPROCESSINGXX"
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

def x_standardize_submission_status__mutmut_22(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "submitted" in submission_status:
		return "SUBMITTED"
	elif "created" in submission_status:
		return "CREATED"
	elif "queued" in submission_status:
		return "QUEUED"
	elif "processing" in submission_status:
		return "processing"
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

def x_standardize_submission_status__mutmut_23(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "submitted" in submission_status:
		return "SUBMITTED"
	elif "created" in submission_status:
		return "CREATED"
	elif "queued" in submission_status:
		return "QUEUED"
	elif "processing" in submission_status:
		return "PROCESSING"
	elif "XXfailedXX" in submission_status:
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

def x_standardize_submission_status__mutmut_24(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "submitted" in submission_status:
		return "SUBMITTED"
	elif "created" in submission_status:
		return "CREATED"
	elif "queued" in submission_status:
		return "QUEUED"
	elif "processing" in submission_status:
		return "PROCESSING"
	elif "FAILED" in submission_status:
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

def x_standardize_submission_status__mutmut_25(submission_status: str) -> str:
	submission_status = submission_status.strip().lower()
	if "submitted" in submission_status:
		return "SUBMITTED"
	elif "created" in submission_status:
		return "CREATED"
	elif "queued" in submission_status:
		return "QUEUED"
	elif "processing" in submission_status:
		return "PROCESSING"
	elif "failed" not in submission_status:
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

def x_standardize_submission_status__mutmut_26(submission_status: str) -> str:
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
		return "XXFAILEDXX"
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

def x_standardize_submission_status__mutmut_27(submission_status: str) -> str:
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
		return "failed"
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

def x_standardize_submission_status__mutmut_28(submission_status: str) -> str:
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
	elif "XXprocessed-okXX" in submission_status:
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

def x_standardize_submission_status__mutmut_29(submission_status: str) -> str:
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
	elif "PROCESSED-OK" in submission_status:
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

def x_standardize_submission_status__mutmut_30(submission_status: str) -> str:
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
	elif "processed-ok" not in submission_status:
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

def x_standardize_submission_status__mutmut_31(submission_status: str) -> str:
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
		return "XXPROCESSEDXX"
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

def x_standardize_submission_status__mutmut_32(submission_status: str) -> str:
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
		return "processed"
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

def x_standardize_submission_status__mutmut_33(submission_status: str) -> str:
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
	elif "XXprocessed-errorXX" in submission_status:
		return "ERROR"
	elif "deleted" in submission_status:
		return "DELETED"
	elif "waiting" in submission_status:
		return "WAITING"
	elif "retried" in submission_status:
		return "RETRIED"
	else:
		return "ERROR"

def x_standardize_submission_status__mutmut_34(submission_status: str) -> str:
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
	elif "PROCESSED-ERROR" in submission_status:
		return "ERROR"
	elif "deleted" in submission_status:
		return "DELETED"
	elif "waiting" in submission_status:
		return "WAITING"
	elif "retried" in submission_status:
		return "RETRIED"
	else:
		return "ERROR"

def x_standardize_submission_status__mutmut_35(submission_status: str) -> str:
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
	elif "processed-error" not in submission_status:
		return "ERROR"
	elif "deleted" in submission_status:
		return "DELETED"
	elif "waiting" in submission_status:
		return "WAITING"
	elif "retried" in submission_status:
		return "RETRIED"
	else:
		return "ERROR"

def x_standardize_submission_status__mutmut_36(submission_status: str) -> str:
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
		return "XXERRORXX"
	elif "deleted" in submission_status:
		return "DELETED"
	elif "waiting" in submission_status:
		return "WAITING"
	elif "retried" in submission_status:
		return "RETRIED"
	else:
		return "ERROR"

def x_standardize_submission_status__mutmut_37(submission_status: str) -> str:
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
		return "error"
	elif "deleted" in submission_status:
		return "DELETED"
	elif "waiting" in submission_status:
		return "WAITING"
	elif "retried" in submission_status:
		return "RETRIED"
	else:
		return "ERROR"

def x_standardize_submission_status__mutmut_38(submission_status: str) -> str:
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
	elif "XXdeletedXX" in submission_status:
		return "DELETED"
	elif "waiting" in submission_status:
		return "WAITING"
	elif "retried" in submission_status:
		return "RETRIED"
	else:
		return "ERROR"

def x_standardize_submission_status__mutmut_39(submission_status: str) -> str:
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
	elif "DELETED" in submission_status:
		return "DELETED"
	elif "waiting" in submission_status:
		return "WAITING"
	elif "retried" in submission_status:
		return "RETRIED"
	else:
		return "ERROR"

def x_standardize_submission_status__mutmut_40(submission_status: str) -> str:
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
	elif "deleted" not in submission_status:
		return "DELETED"
	elif "waiting" in submission_status:
		return "WAITING"
	elif "retried" in submission_status:
		return "RETRIED"
	else:
		return "ERROR"

def x_standardize_submission_status__mutmut_41(submission_status: str) -> str:
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
		return "XXDELETEDXX"
	elif "waiting" in submission_status:
		return "WAITING"
	elif "retried" in submission_status:
		return "RETRIED"
	else:
		return "ERROR"

def x_standardize_submission_status__mutmut_42(submission_status: str) -> str:
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
		return "deleted"
	elif "waiting" in submission_status:
		return "WAITING"
	elif "retried" in submission_status:
		return "RETRIED"
	else:
		return "ERROR"

def x_standardize_submission_status__mutmut_43(submission_status: str) -> str:
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
	elif "XXwaitingXX" in submission_status:
		return "WAITING"
	elif "retried" in submission_status:
		return "RETRIED"
	else:
		return "ERROR"

def x_standardize_submission_status__mutmut_44(submission_status: str) -> str:
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
	elif "WAITING" in submission_status:
		return "WAITING"
	elif "retried" in submission_status:
		return "RETRIED"
	else:
		return "ERROR"

def x_standardize_submission_status__mutmut_45(submission_status: str) -> str:
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
	elif "waiting" not in submission_status:
		return "WAITING"
	elif "retried" in submission_status:
		return "RETRIED"
	else:
		return "ERROR"

def x_standardize_submission_status__mutmut_46(submission_status: str) -> str:
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
		return "XXWAITINGXX"
	elif "retried" in submission_status:
		return "RETRIED"
	else:
		return "ERROR"

def x_standardize_submission_status__mutmut_47(submission_status: str) -> str:
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
		return "waiting"
	elif "retried" in submission_status:
		return "RETRIED"
	else:
		return "ERROR"

def x_standardize_submission_status__mutmut_48(submission_status: str) -> str:
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
	elif "XXretriedXX" in submission_status:
		return "RETRIED"
	else:
		return "ERROR"

def x_standardize_submission_status__mutmut_49(submission_status: str) -> str:
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
	elif "RETRIED" in submission_status:
		return "RETRIED"
	else:
		return "ERROR"

def x_standardize_submission_status__mutmut_50(submission_status: str) -> str:
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
	elif "retried" not in submission_status:
		return "RETRIED"
	else:
		return "ERROR"

def x_standardize_submission_status__mutmut_51(submission_status: str) -> str:
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
		return "XXRETRIEDXX"
	else:
		return "ERROR"

def x_standardize_submission_status__mutmut_52(submission_status: str) -> str:
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
		return "retried"
	else:
		return "ERROR"

def x_standardize_submission_status__mutmut_53(submission_status: str) -> str:
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
		return "XXERRORXX"

def x_standardize_submission_status__mutmut_54(submission_status: str) -> str:
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
		return "error"

x_standardize_submission_status__mutmut_mutants : ClassVar[MutantDict] = { # type: ignore
'x_standardize_submission_status__mutmut_1': x_standardize_submission_status__mutmut_1, 
    'x_standardize_submission_status__mutmut_2': x_standardize_submission_status__mutmut_2, 
    'x_standardize_submission_status__mutmut_3': x_standardize_submission_status__mutmut_3, 
    'x_standardize_submission_status__mutmut_4': x_standardize_submission_status__mutmut_4, 
    'x_standardize_submission_status__mutmut_5': x_standardize_submission_status__mutmut_5, 
    'x_standardize_submission_status__mutmut_6': x_standardize_submission_status__mutmut_6, 
    'x_standardize_submission_status__mutmut_7': x_standardize_submission_status__mutmut_7, 
    'x_standardize_submission_status__mutmut_8': x_standardize_submission_status__mutmut_8, 
    'x_standardize_submission_status__mutmut_9': x_standardize_submission_status__mutmut_9, 
    'x_standardize_submission_status__mutmut_10': x_standardize_submission_status__mutmut_10, 
    'x_standardize_submission_status__mutmut_11': x_standardize_submission_status__mutmut_11, 
    'x_standardize_submission_status__mutmut_12': x_standardize_submission_status__mutmut_12, 
    'x_standardize_submission_status__mutmut_13': x_standardize_submission_status__mutmut_13, 
    'x_standardize_submission_status__mutmut_14': x_standardize_submission_status__mutmut_14, 
    'x_standardize_submission_status__mutmut_15': x_standardize_submission_status__mutmut_15, 
    'x_standardize_submission_status__mutmut_16': x_standardize_submission_status__mutmut_16, 
    'x_standardize_submission_status__mutmut_17': x_standardize_submission_status__mutmut_17, 
    'x_standardize_submission_status__mutmut_18': x_standardize_submission_status__mutmut_18, 
    'x_standardize_submission_status__mutmut_19': x_standardize_submission_status__mutmut_19, 
    'x_standardize_submission_status__mutmut_20': x_standardize_submission_status__mutmut_20, 
    'x_standardize_submission_status__mutmut_21': x_standardize_submission_status__mutmut_21, 
    'x_standardize_submission_status__mutmut_22': x_standardize_submission_status__mutmut_22, 
    'x_standardize_submission_status__mutmut_23': x_standardize_submission_status__mutmut_23, 
    'x_standardize_submission_status__mutmut_24': x_standardize_submission_status__mutmut_24, 
    'x_standardize_submission_status__mutmut_25': x_standardize_submission_status__mutmut_25, 
    'x_standardize_submission_status__mutmut_26': x_standardize_submission_status__mutmut_26, 
    'x_standardize_submission_status__mutmut_27': x_standardize_submission_status__mutmut_27, 
    'x_standardize_submission_status__mutmut_28': x_standardize_submission_status__mutmut_28, 
    'x_standardize_submission_status__mutmut_29': x_standardize_submission_status__mutmut_29, 
    'x_standardize_submission_status__mutmut_30': x_standardize_submission_status__mutmut_30, 
    'x_standardize_submission_status__mutmut_31': x_standardize_submission_status__mutmut_31, 
    'x_standardize_submission_status__mutmut_32': x_standardize_submission_status__mutmut_32, 
    'x_standardize_submission_status__mutmut_33': x_standardize_submission_status__mutmut_33, 
    'x_standardize_submission_status__mutmut_34': x_standardize_submission_status__mutmut_34, 
    'x_standardize_submission_status__mutmut_35': x_standardize_submission_status__mutmut_35, 
    'x_standardize_submission_status__mutmut_36': x_standardize_submission_status__mutmut_36, 
    'x_standardize_submission_status__mutmut_37': x_standardize_submission_status__mutmut_37, 
    'x_standardize_submission_status__mutmut_38': x_standardize_submission_status__mutmut_38, 
    'x_standardize_submission_status__mutmut_39': x_standardize_submission_status__mutmut_39, 
    'x_standardize_submission_status__mutmut_40': x_standardize_submission_status__mutmut_40, 
    'x_standardize_submission_status__mutmut_41': x_standardize_submission_status__mutmut_41, 
    'x_standardize_submission_status__mutmut_42': x_standardize_submission_status__mutmut_42, 
    'x_standardize_submission_status__mutmut_43': x_standardize_submission_status__mutmut_43, 
    'x_standardize_submission_status__mutmut_44': x_standardize_submission_status__mutmut_44, 
    'x_standardize_submission_status__mutmut_45': x_standardize_submission_status__mutmut_45, 
    'x_standardize_submission_status__mutmut_46': x_standardize_submission_status__mutmut_46, 
    'x_standardize_submission_status__mutmut_47': x_standardize_submission_status__mutmut_47, 
    'x_standardize_submission_status__mutmut_48': x_standardize_submission_status__mutmut_48, 
    'x_standardize_submission_status__mutmut_49': x_standardize_submission_status__mutmut_49, 
    'x_standardize_submission_status__mutmut_50': x_standardize_submission_status__mutmut_50, 
    'x_standardize_submission_status__mutmut_51': x_standardize_submission_status__mutmut_51, 
    'x_standardize_submission_status__mutmut_52': x_standardize_submission_status__mutmut_52, 
    'x_standardize_submission_status__mutmut_53': x_standardize_submission_status__mutmut_53, 
    'x_standardize_submission_status__mutmut_54': x_standardize_submission_status__mutmut_54
}
x_standardize_submission_status__mutmut_orig.__name__ = 'x_standardize_submission_status'

def process_report_header(report_file: str) -> tuple[dict[str, Any], str, str]:
	args = [report_file]# type: ignore
	kwargs = {}# type: ignore
	return _mutmut_trampoline(x_process_report_header__mutmut_orig, x_process_report_header__mutmut_mutants, args, kwargs, None)

def x_process_report_header__mutmut_orig(report_file: str) -> tuple[dict[str, Any], str, str]:
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

def x_process_report_header__mutmut_1(report_file: str) -> tuple[dict[str, Any], str, str]:
	# Read in report.xml
	tree = None
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

def x_process_report_header__mutmut_2(report_file: str) -> tuple[dict[str, Any], str, str]:
	# Read in report.xml
	tree = ET.parse(None)
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

def x_process_report_header__mutmut_3(report_file: str) -> tuple[dict[str, Any], str, str]:
	# Read in report.xml
	tree = ET.parse(report_file)
	root = None
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

def x_process_report_header__mutmut_4(report_file: str) -> tuple[dict[str, Any], str, str]:
	# Read in report.xml
	tree = ET.parse(report_file)
	root = tree.getroot()
	xmlstr = None
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

def x_process_report_header__mutmut_5(report_file: str) -> tuple[dict[str, Any], str, str]:
	# Read in report.xml
	tree = ET.parse(report_file)
	root = tree.getroot()
	xmlstr = ET.tostring(None, encoding='utf8', method='xml')
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

def x_process_report_header__mutmut_6(report_file: str) -> tuple[dict[str, Any], str, str]:
	# Read in report.xml
	tree = ET.parse(report_file)
	root = tree.getroot()
	xmlstr = ET.tostring(root, encoding=None, method='xml')
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

def x_process_report_header__mutmut_7(report_file: str) -> tuple[dict[str, Any], str, str]:
	# Read in report.xml
	tree = ET.parse(report_file)
	root = tree.getroot()
	xmlstr = ET.tostring(root, encoding='utf8', method=None)
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

def x_process_report_header__mutmut_8(report_file: str) -> tuple[dict[str, Any], str, str]:
	# Read in report.xml
	tree = ET.parse(report_file)
	root = tree.getroot()
	xmlstr = ET.tostring(encoding='utf8', method='xml')
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

def x_process_report_header__mutmut_9(report_file: str) -> tuple[dict[str, Any], str, str]:
	# Read in report.xml
	tree = ET.parse(report_file)
	root = tree.getroot()
	xmlstr = ET.tostring(root, method='xml')
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

def x_process_report_header__mutmut_10(report_file: str) -> tuple[dict[str, Any], str, str]:
	# Read in report.xml
	tree = ET.parse(report_file)
	root = tree.getroot()
	xmlstr = ET.tostring(root, encoding='utf8', )
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

def x_process_report_header__mutmut_11(report_file: str) -> tuple[dict[str, Any], str, str]:
	# Read in report.xml
	tree = ET.parse(report_file)
	root = tree.getroot()
	xmlstr = ET.tostring(root, encoding='XXutf8XX', method='xml')
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

def x_process_report_header__mutmut_12(report_file: str) -> tuple[dict[str, Any], str, str]:
	# Read in report.xml
	tree = ET.parse(report_file)
	root = tree.getroot()
	xmlstr = ET.tostring(root, encoding='UTF8', method='xml')
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

def x_process_report_header__mutmut_13(report_file: str) -> tuple[dict[str, Any], str, str]:
	# Read in report.xml
	tree = ET.parse(report_file)
	root = tree.getroot()
	xmlstr = ET.tostring(root, encoding='utf8', method='XXxmlXX')
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

def x_process_report_header__mutmut_14(report_file: str) -> tuple[dict[str, Any], str, str]:
	# Read in report.xml
	tree = ET.parse(report_file)
	root = tree.getroot()
	xmlstr = ET.tostring(root, encoding='utf8', method='XML')
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

def x_process_report_header__mutmut_15(report_file: str) -> tuple[dict[str, Any], str, str]:
	# Read in report.xml
	tree = ET.parse(report_file)
	root = tree.getroot()
	xmlstr = ET.tostring(root, encoding='utf8', method='xml')
	# Convert xml to dictionary
	report_dict = None
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

def x_process_report_header__mutmut_16(report_file: str) -> tuple[dict[str, Any], str, str]:
	# Read in report.xml
	tree = ET.parse(report_file)
	root = tree.getroot()
	xmlstr = ET.tostring(root, encoding='utf8', method='xml')
	# Convert xml to dictionary
	report_dict = xmltodict.parse(None)
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

def x_process_report_header__mutmut_17(report_file: str) -> tuple[dict[str, Any], str, str]:
	# Read in report.xml
	tree = ET.parse(report_file)
	root = tree.getroot()
	xmlstr = ET.tostring(root, encoding='utf8', method='xml')
	# Convert xml to dictionary
	report_dict = xmltodict.parse(xmlstr)
	# Get submission status
	try:
		# Get submission status and id from report.xml
		submission_status = None
	except KeyError:
		submission_status = "SUBMITTED"
	submission_status = standardize_submission_status(submission_status=submission_status)
	# Get submission id
	try:
		submission_id = report_dict["SubmissionStatus"]["@submission_id"]
	except KeyError:
		submission_id = "PENDING"
	return report_dict, submission_status, submission_id

def x_process_report_header__mutmut_18(report_file: str) -> tuple[dict[str, Any], str, str]:
	# Read in report.xml
	tree = ET.parse(report_file)
	root = tree.getroot()
	xmlstr = ET.tostring(root, encoding='utf8', method='xml')
	# Convert xml to dictionary
	report_dict = xmltodict.parse(xmlstr)
	# Get submission status
	try:
		# Get submission status and id from report.xml
		submission_status = report_dict["XXSubmissionStatusXX"]["@status"]
	except KeyError:
		submission_status = "SUBMITTED"
	submission_status = standardize_submission_status(submission_status=submission_status)
	# Get submission id
	try:
		submission_id = report_dict["SubmissionStatus"]["@submission_id"]
	except KeyError:
		submission_id = "PENDING"
	return report_dict, submission_status, submission_id

def x_process_report_header__mutmut_19(report_file: str) -> tuple[dict[str, Any], str, str]:
	# Read in report.xml
	tree = ET.parse(report_file)
	root = tree.getroot()
	xmlstr = ET.tostring(root, encoding='utf8', method='xml')
	# Convert xml to dictionary
	report_dict = xmltodict.parse(xmlstr)
	# Get submission status
	try:
		# Get submission status and id from report.xml
		submission_status = report_dict["submissionstatus"]["@status"]
	except KeyError:
		submission_status = "SUBMITTED"
	submission_status = standardize_submission_status(submission_status=submission_status)
	# Get submission id
	try:
		submission_id = report_dict["SubmissionStatus"]["@submission_id"]
	except KeyError:
		submission_id = "PENDING"
	return report_dict, submission_status, submission_id

def x_process_report_header__mutmut_20(report_file: str) -> tuple[dict[str, Any], str, str]:
	# Read in report.xml
	tree = ET.parse(report_file)
	root = tree.getroot()
	xmlstr = ET.tostring(root, encoding='utf8', method='xml')
	# Convert xml to dictionary
	report_dict = xmltodict.parse(xmlstr)
	# Get submission status
	try:
		# Get submission status and id from report.xml
		submission_status = report_dict["SUBMISSIONSTATUS"]["@status"]
	except KeyError:
		submission_status = "SUBMITTED"
	submission_status = standardize_submission_status(submission_status=submission_status)
	# Get submission id
	try:
		submission_id = report_dict["SubmissionStatus"]["@submission_id"]
	except KeyError:
		submission_id = "PENDING"
	return report_dict, submission_status, submission_id

def x_process_report_header__mutmut_21(report_file: str) -> tuple[dict[str, Any], str, str]:
	# Read in report.xml
	tree = ET.parse(report_file)
	root = tree.getroot()
	xmlstr = ET.tostring(root, encoding='utf8', method='xml')
	# Convert xml to dictionary
	report_dict = xmltodict.parse(xmlstr)
	# Get submission status
	try:
		# Get submission status and id from report.xml
		submission_status = report_dict["SubmissionStatus"]["XX@statusXX"]
	except KeyError:
		submission_status = "SUBMITTED"
	submission_status = standardize_submission_status(submission_status=submission_status)
	# Get submission id
	try:
		submission_id = report_dict["SubmissionStatus"]["@submission_id"]
	except KeyError:
		submission_id = "PENDING"
	return report_dict, submission_status, submission_id

def x_process_report_header__mutmut_22(report_file: str) -> tuple[dict[str, Any], str, str]:
	# Read in report.xml
	tree = ET.parse(report_file)
	root = tree.getroot()
	xmlstr = ET.tostring(root, encoding='utf8', method='xml')
	# Convert xml to dictionary
	report_dict = xmltodict.parse(xmlstr)
	# Get submission status
	try:
		# Get submission status and id from report.xml
		submission_status = report_dict["SubmissionStatus"]["@STATUS"]
	except KeyError:
		submission_status = "SUBMITTED"
	submission_status = standardize_submission_status(submission_status=submission_status)
	# Get submission id
	try:
		submission_id = report_dict["SubmissionStatus"]["@submission_id"]
	except KeyError:
		submission_id = "PENDING"
	return report_dict, submission_status, submission_id

def x_process_report_header__mutmut_23(report_file: str) -> tuple[dict[str, Any], str, str]:
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
		submission_status = None
	submission_status = standardize_submission_status(submission_status=submission_status)
	# Get submission id
	try:
		submission_id = report_dict["SubmissionStatus"]["@submission_id"]
	except KeyError:
		submission_id = "PENDING"
	return report_dict, submission_status, submission_id

def x_process_report_header__mutmut_24(report_file: str) -> tuple[dict[str, Any], str, str]:
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
		submission_status = "XXSUBMITTEDXX"
	submission_status = standardize_submission_status(submission_status=submission_status)
	# Get submission id
	try:
		submission_id = report_dict["SubmissionStatus"]["@submission_id"]
	except KeyError:
		submission_id = "PENDING"
	return report_dict, submission_status, submission_id

def x_process_report_header__mutmut_25(report_file: str) -> tuple[dict[str, Any], str, str]:
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
		submission_status = "submitted"
	submission_status = standardize_submission_status(submission_status=submission_status)
	# Get submission id
	try:
		submission_id = report_dict["SubmissionStatus"]["@submission_id"]
	except KeyError:
		submission_id = "PENDING"
	return report_dict, submission_status, submission_id

def x_process_report_header__mutmut_26(report_file: str) -> tuple[dict[str, Any], str, str]:
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
	submission_status = None
	# Get submission id
	try:
		submission_id = report_dict["SubmissionStatus"]["@submission_id"]
	except KeyError:
		submission_id = "PENDING"
	return report_dict, submission_status, submission_id

def x_process_report_header__mutmut_27(report_file: str) -> tuple[dict[str, Any], str, str]:
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
	submission_status = standardize_submission_status(submission_status=None)
	# Get submission id
	try:
		submission_id = report_dict["SubmissionStatus"]["@submission_id"]
	except KeyError:
		submission_id = "PENDING"
	return report_dict, submission_status, submission_id

def x_process_report_header__mutmut_28(report_file: str) -> tuple[dict[str, Any], str, str]:
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
		submission_id = None
	except KeyError:
		submission_id = "PENDING"
	return report_dict, submission_status, submission_id

def x_process_report_header__mutmut_29(report_file: str) -> tuple[dict[str, Any], str, str]:
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
		submission_id = report_dict["XXSubmissionStatusXX"]["@submission_id"]
	except KeyError:
		submission_id = "PENDING"
	return report_dict, submission_status, submission_id

def x_process_report_header__mutmut_30(report_file: str) -> tuple[dict[str, Any], str, str]:
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
		submission_id = report_dict["submissionstatus"]["@submission_id"]
	except KeyError:
		submission_id = "PENDING"
	return report_dict, submission_status, submission_id

def x_process_report_header__mutmut_31(report_file: str) -> tuple[dict[str, Any], str, str]:
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
		submission_id = report_dict["SUBMISSIONSTATUS"]["@submission_id"]
	except KeyError:
		submission_id = "PENDING"
	return report_dict, submission_status, submission_id

def x_process_report_header__mutmut_32(report_file: str) -> tuple[dict[str, Any], str, str]:
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
		submission_id = report_dict["SubmissionStatus"]["XX@submission_idXX"]
	except KeyError:
		submission_id = "PENDING"
	return report_dict, submission_status, submission_id

def x_process_report_header__mutmut_33(report_file: str) -> tuple[dict[str, Any], str, str]:
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
		submission_id = report_dict["SubmissionStatus"]["@SUBMISSION_ID"]
	except KeyError:
		submission_id = "PENDING"
	return report_dict, submission_status, submission_id

def x_process_report_header__mutmut_34(report_file: str) -> tuple[dict[str, Any], str, str]:
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
		submission_id = None
	return report_dict, submission_status, submission_id

def x_process_report_header__mutmut_35(report_file: str) -> tuple[dict[str, Any], str, str]:
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
		submission_id = "XXPENDINGXX"
	return report_dict, submission_status, submission_id

def x_process_report_header__mutmut_36(report_file: str) -> tuple[dict[str, Any], str, str]:
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
		submission_id = "pending"
	return report_dict, submission_status, submission_id

x_process_report_header__mutmut_mutants : ClassVar[MutantDict] = { # type: ignore
'x_process_report_header__mutmut_1': x_process_report_header__mutmut_1, 
    'x_process_report_header__mutmut_2': x_process_report_header__mutmut_2, 
    'x_process_report_header__mutmut_3': x_process_report_header__mutmut_3, 
    'x_process_report_header__mutmut_4': x_process_report_header__mutmut_4, 
    'x_process_report_header__mutmut_5': x_process_report_header__mutmut_5, 
    'x_process_report_header__mutmut_6': x_process_report_header__mutmut_6, 
    'x_process_report_header__mutmut_7': x_process_report_header__mutmut_7, 
    'x_process_report_header__mutmut_8': x_process_report_header__mutmut_8, 
    'x_process_report_header__mutmut_9': x_process_report_header__mutmut_9, 
    'x_process_report_header__mutmut_10': x_process_report_header__mutmut_10, 
    'x_process_report_header__mutmut_11': x_process_report_header__mutmut_11, 
    'x_process_report_header__mutmut_12': x_process_report_header__mutmut_12, 
    'x_process_report_header__mutmut_13': x_process_report_header__mutmut_13, 
    'x_process_report_header__mutmut_14': x_process_report_header__mutmut_14, 
    'x_process_report_header__mutmut_15': x_process_report_header__mutmut_15, 
    'x_process_report_header__mutmut_16': x_process_report_header__mutmut_16, 
    'x_process_report_header__mutmut_17': x_process_report_header__mutmut_17, 
    'x_process_report_header__mutmut_18': x_process_report_header__mutmut_18, 
    'x_process_report_header__mutmut_19': x_process_report_header__mutmut_19, 
    'x_process_report_header__mutmut_20': x_process_report_header__mutmut_20, 
    'x_process_report_header__mutmut_21': x_process_report_header__mutmut_21, 
    'x_process_report_header__mutmut_22': x_process_report_header__mutmut_22, 
    'x_process_report_header__mutmut_23': x_process_report_header__mutmut_23, 
    'x_process_report_header__mutmut_24': x_process_report_header__mutmut_24, 
    'x_process_report_header__mutmut_25': x_process_report_header__mutmut_25, 
    'x_process_report_header__mutmut_26': x_process_report_header__mutmut_26, 
    'x_process_report_header__mutmut_27': x_process_report_header__mutmut_27, 
    'x_process_report_header__mutmut_28': x_process_report_header__mutmut_28, 
    'x_process_report_header__mutmut_29': x_process_report_header__mutmut_29, 
    'x_process_report_header__mutmut_30': x_process_report_header__mutmut_30, 
    'x_process_report_header__mutmut_31': x_process_report_header__mutmut_31, 
    'x_process_report_header__mutmut_32': x_process_report_header__mutmut_32, 
    'x_process_report_header__mutmut_33': x_process_report_header__mutmut_33, 
    'x_process_report_header__mutmut_34': x_process_report_header__mutmut_34, 
    'x_process_report_header__mutmut_35': x_process_report_header__mutmut_35, 
    'x_process_report_header__mutmut_36': x_process_report_header__mutmut_36
}
x_process_report_header__mutmut_orig.__name__ = 'x_process_report_header'
