#!/usr/bin/env python3

###########################    Description    ##################################
# Global values for SeqSender functions
################################################################################

import os
from typing import List, Dict

##### SeqSender settings #####
# Script directory
PROG_DIR: str = os.path.dirname(os.path.abspath(__file__))

# SeqSender version
VERSION: str = "1.2.3 (Beta)"

# Organism options with unique submission options
ORGANISM_CHOICES: List[str] = ["FLU", "COV", "POX", "ARBO", "OTHER"]

# Database submisison options
DATABASE_CHOICES: List[str] = ["BIOSAMPLE", "SRA", "GENBANK", "GISAID"]

# metadata prefix for each database
SAMPLE_NAME_DATABASE_PREFIX: Dict[str, str] = {"BIOSAMPLE":"bs-", "SRA":"sra-", "GENBANK":"gb-", "GISAID":"gs-"}

# Submission status report columns
BIOSAMPLE_SUBMISSION_STATUS_COLUMNS: List[str] = ["biosample_status", "biosample_accession", "biosample_message"]
SRA_SUBMISSION_STATUS_COLUMNS: List[str] = ["sra_status", "sra_accession", "sra_message"]
GENBANK_SUBMISSION_STATUS_COLUMNS: List[str] = ["genbank_status", "genbank_accession", "genbank_message"]
GISAID_SUBMISSION_STATUS_COLUMNS: List[str] = ["gisaid_accession_epi_isl_id", "gisaid_accession_epi_id", "gisaid_message"]

# Upload log columns
SUBMISSION_LOG_COLUMNS: List[str] = ["Submission_Name", "Organism", "Database", "Submission_Type", "Submission_Date", "Submission_ID", "Submission_Status", "Submission_Directory", "Config_File", "Update_Date"]

# Shiny schema options, exclusion list
SCHEMA_EXCLUSIONS = ["config.seqsender.upload_log_schema","config_file.ncbi_schema","config_file.ncbi_gisaid_schema", "config_file.gisaid_schema"]

##### NCBI settings #####
# FTP website to submit samples to
NCBI_FTP_HOST: str = "sftp-private.ncbi.nlm.nih.gov"
NCBI_FTP_PORT: int = 22
NCBI_FTP_PUBLIC_KEY: bytes = b"""AAAAB3NzaC1yc2EAAAADAQABAAABAQCjEuJc/1pAQxGxwuk8s/1u7tDq2PiOqutAXNjPOl0BkOiJJnUa6gdRx5WcVd6HApPrw7RUvUXraOmc82Zlzeno6DkzTjLk3ztuRs1hNF3KKx044AxDC1Jew6pK96At3kAvStfamNbY8URB9tj0SF3rFQxLrWo8WfhFr9AU4SuY5S0A2a/njucmTim3vDRPbBFO7Dlz2dxzeKPBdUDB/NvmdPWHxkviXkqUby/QzjglAQTts5lNEMH8qaRcU/dnWYRKvY17RHGu9/gNf8Z9DvBg8TbU7Ek3MT5A9F3nw1ZcZe7M8L7JwIKgQUl1gzOBHFfC248YES9REQffLYPWOCy9"""

# URL structure to download NCBI output files
NCBI_API_URL: str = "https://submit.ncbi.nlm.nih.gov/api/2.0/files/FILE_ID/?format=attachment"

# Table2asn email to submit samples to
TABLE2ASN_EMAIL:str  = "gb-admin@ncbi.nlm.nih.gov"

# BioSample metadata regex
BIOSAMPLE_REGEX = "^bs-|^bioproject$|^organism$|^collection_date$"

# SRA metadata regex
SRA_REGEX = "^sra-|^bioproject$|bs-sample_name|^organism$|^collection_date$"

# Genbank metadata regex
GENBANK_REGEX = "^gb-sample_name$"

# GenBank source file metadata regex
GENBANK_REGEX_SRC = "^gb-sample_name$|^src-|^bioproject$|^organism$|^collection_date$"

# GenBank comment file metadata regex
GENBANK_REGEX_CMT = "^gb-sample_name$|^cmt-|^organism$|^collection_date$"

##### GISAID settings #####
# GISAID metadata regex
GISAID_REGEX = "^gs-|^collection_date$|^authors$"
