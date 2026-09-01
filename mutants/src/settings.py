#!/usr/bin/env python3

###########################    Description    ##################################
# Global values for SeqSender functions
################################################################################

import os

##### SeqSender settings #####
# Script directory
PROG_DIR: str = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# SeqSender version
VERSION: str = "1.5.0 (Beta)"

# Organism options with unique submission options
ORGANISM_CHOICES: list[str] = ["FLU", "COV", "OTHER"]

# Database submisison options
DATABASE_CHOICES: list[str] = ["BIOSAMPLE", "SRA", "GENBANK"]

# metadata prefix for each database
SAMPLE_NAME_DATABASE_PREFIX: dict[str, str] = {"BIOSAMPLE":"bs-", "SRA":"sra-", "GENBANK":"gb-"}

# Submission status report columns
BIOSAMPLE_SUBMISSION_STATUS_COLUMNS: list[str] = ["biosample_status", "biosample_accession", "biosample_message"]
SRA_SUBMISSION_STATUS_COLUMNS: list[str] = ["sra_status", "sra_accession", "sra_message"]
GENBANK_SUBMISSION_STATUS_COLUMNS: list[str] = ["genbank_status", "genbank_accession", "genbank_message"]

# Upload log columns
SUBMISSION_LOG_COLUMNS: list[str] = ["Submission_Name", "Organism", "Database", "Submission_Type", "Submission_Date", "Submission_ID", "Submission_Status", "Submission_Directory", "Config_File", "Update_Date"]

# Shiny schema options, exclusion list
SCHEMA_EXCLUSIONS = ["config.seqsender.upload_log_schema","config_file.ncbi_schema"]

##### NCBI settings #####
# FTP website to submit samples to
NCBI_FTP_HOST: str = "ftp-private.ncbi.nlm.nih.gov"

# URL structure to download NCBI output files
NCBI_API_URL: str = "https://submit.ncbi.nlm.nih.gov/api/2.0/files/FILE_ID/?format=attachment"

# Table2asn email to submit samples to
TABLE2ASN_EMAIL:str  = "gb-admin@ncbi.nlm.nih.gov"

# GenBank FTP options based on organism
GENBANK_FTP_ORGANISMS: list[str] = ["FLU", "COV"]

# BioSample metadata regex
BIOSAMPLE_REGEX = "^bs-|^bioproject$|^organism$|^collection_date$"

# SRA metadata regex
SRA_REGEX = "^sra-|^bioproject$|bs-sample_name|^organism$|^collection_date$"

# Genbank metadata regex
GENBANK_REGEX = "^gb-sample_name$"

# GenBank source file metadata regex
GENBANK_REGEX_SRC = "^gb-sample_name$|^src-|^bioproject$|^organism$|^collection_date$"

# GenBank comment file metadata regex
GENBANK_REGEX_CMT = "^gb-sample_name$|^cmt-"

# Deprecated GenBank columns not allowed to be used
GENBANK_DEPRECATED_COLUMNS: list[str] = ["src-Authority", "src-Biotype", "src-Biovar", "src-Chemovar", "src-Forma", "src-Forma_specialis", "src-Identified_by", "src-Pathovar", "src-Pop_variant", "src-Serogroup", "src-Subclone", "src-Subtype", "src-Substrain", "src-Type"]
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
