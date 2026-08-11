#!/usr/bin/env python3

###########################	Description	##################################
# Functions to handle validating file locations and loading files
################################################################################

import sys
import os
from typing import Dict, Any
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import Optional
import shutil
import yaml
import time

from src.settings import SAMPLE_NAME_DATABASE_PREFIX, PROG_DIR
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

def copy_file(source: str, destination: str):
	args = [source, destination]# type: ignore
	kwargs = {}# type: ignore
	return _mutmut_trampoline(x_copy_file__mutmut_orig, x_copy_file__mutmut_mutants, args, kwargs, None)

def x_copy_file__mutmut_orig(source: str, destination: str):
	shutil.copy(source, destination)

def x_copy_file__mutmut_1(source: str, destination: str):
	shutil.copy(None, destination)

def x_copy_file__mutmut_2(source: str, destination: str):
	shutil.copy(source, None)

def x_copy_file__mutmut_3(source: str, destination: str):
	shutil.copy(destination)

def x_copy_file__mutmut_4(source: str, destination: str):
	shutil.copy(source, )

x_copy_file__mutmut_mutants : ClassVar[MutantDict] = { # type: ignore
'x_copy_file__mutmut_1': x_copy_file__mutmut_1, 
    'x_copy_file__mutmut_2': x_copy_file__mutmut_2, 
    'x_copy_file__mutmut_3': x_copy_file__mutmut_3, 
    'x_copy_file__mutmut_4': x_copy_file__mutmut_4
}
x_copy_file__mutmut_orig.__name__ = 'x_copy_file'

# Validate file exists or error out
def validate_file(file_type: str, file_path: str):
	args = [file_type, file_path]# type: ignore
	kwargs = {}# type: ignore
	return _mutmut_trampoline(x_validate_file__mutmut_orig, x_validate_file__mutmut_mutants, args, kwargs, None)

# Validate file exists or error out
def x_validate_file__mutmut_orig(file_type: str, file_path: str):
	if not os.path.isfile(file_path):
		print(f"Error: Input {file_type.replace('_',' ')} does not exist at: {file_path}", file=sys.stderr)
		sys.exit(1)

# Validate file exists or error out
def x_validate_file__mutmut_1(file_type: str, file_path: str):
	if os.path.isfile(file_path):
		print(f"Error: Input {file_type.replace('_',' ')} does not exist at: {file_path}", file=sys.stderr)
		sys.exit(1)

# Validate file exists or error out
def x_validate_file__mutmut_2(file_type: str, file_path: str):
	if not os.path.isfile(None):
		print(f"Error: Input {file_type.replace('_',' ')} does not exist at: {file_path}", file=sys.stderr)
		sys.exit(1)

# Validate file exists or error out
def x_validate_file__mutmut_3(file_type: str, file_path: str):
	if not os.path.isfile(file_path):
		print(None, file=sys.stderr)
		sys.exit(1)

# Validate file exists or error out
def x_validate_file__mutmut_4(file_type: str, file_path: str):
	if not os.path.isfile(file_path):
		print(f"Error: Input {file_type.replace('_',' ')} does not exist at: {file_path}", file=None)
		sys.exit(1)

# Validate file exists or error out
def x_validate_file__mutmut_5(file_type: str, file_path: str):
	if not os.path.isfile(file_path):
		print(file=sys.stderr)
		sys.exit(1)

# Validate file exists or error out
def x_validate_file__mutmut_6(file_type: str, file_path: str):
	if not os.path.isfile(file_path):
		print(f"Error: Input {file_type.replace('_',' ')} does not exist at: {file_path}", )
		sys.exit(1)

# Validate file exists or error out
def x_validate_file__mutmut_7(file_type: str, file_path: str):
	if not os.path.isfile(file_path):
		print(f"Error: Input {file_type.replace(None,' ')} does not exist at: {file_path}", file=sys.stderr)
		sys.exit(1)

# Validate file exists or error out
def x_validate_file__mutmut_8(file_type: str, file_path: str):
	if not os.path.isfile(file_path):
		print(f"Error: Input {file_type.replace('_',None)} does not exist at: {file_path}", file=sys.stderr)
		sys.exit(1)

# Validate file exists or error out
def x_validate_file__mutmut_9(file_type: str, file_path: str):
	if not os.path.isfile(file_path):
		print(f"Error: Input {file_type.replace(' ')} does not exist at: {file_path}", file=sys.stderr)
		sys.exit(1)

# Validate file exists or error out
def x_validate_file__mutmut_10(file_type: str, file_path: str):
	if not os.path.isfile(file_path):
		print(f"Error: Input {file_type.replace('_',)} does not exist at: {file_path}", file=sys.stderr)
		sys.exit(1)

# Validate file exists or error out
def x_validate_file__mutmut_11(file_type: str, file_path: str):
	if not os.path.isfile(file_path):
		print(f"Error: Input {file_type.replace('XX_XX',' ')} does not exist at: {file_path}", file=sys.stderr)
		sys.exit(1)

# Validate file exists or error out
def x_validate_file__mutmut_12(file_type: str, file_path: str):
	if not os.path.isfile(file_path):
		print(f"Error: Input {file_type.replace('_','XX XX')} does not exist at: {file_path}", file=sys.stderr)
		sys.exit(1)

# Validate file exists or error out
def x_validate_file__mutmut_13(file_type: str, file_path: str):
	if not os.path.isfile(file_path):
		print(f"Error: Input {file_type.replace('_',' ')} does not exist at: {file_path}", file=sys.stderr)
		sys.exit(None)

# Validate file exists or error out
def x_validate_file__mutmut_14(file_type: str, file_path: str):
	if not os.path.isfile(file_path):
		print(f"Error: Input {file_type.replace('_',' ')} does not exist at: {file_path}", file=sys.stderr)
		sys.exit(2)

x_validate_file__mutmut_mutants : ClassVar[MutantDict] = { # type: ignore
'x_validate_file__mutmut_1': x_validate_file__mutmut_1, 
    'x_validate_file__mutmut_2': x_validate_file__mutmut_2, 
    'x_validate_file__mutmut_3': x_validate_file__mutmut_3, 
    'x_validate_file__mutmut_4': x_validate_file__mutmut_4, 
    'x_validate_file__mutmut_5': x_validate_file__mutmut_5, 
    'x_validate_file__mutmut_6': x_validate_file__mutmut_6, 
    'x_validate_file__mutmut_7': x_validate_file__mutmut_7, 
    'x_validate_file__mutmut_8': x_validate_file__mutmut_8, 
    'x_validate_file__mutmut_9': x_validate_file__mutmut_9, 
    'x_validate_file__mutmut_10': x_validate_file__mutmut_10, 
    'x_validate_file__mutmut_11': x_validate_file__mutmut_11, 
    'x_validate_file__mutmut_12': x_validate_file__mutmut_12, 
    'x_validate_file__mutmut_13': x_validate_file__mutmut_13, 
    'x_validate_file__mutmut_14': x_validate_file__mutmut_14
}
x_validate_file__mutmut_orig.__name__ = 'x_validate_file'

def validate_directory(name: str, path: str):
	args = [name, path]# type: ignore
	kwargs = {}# type: ignore
	return _mutmut_trampoline(x_validate_directory__mutmut_orig, x_validate_directory__mutmut_mutants, args, kwargs, None)

def x_validate_directory__mutmut_orig(name: str, path: str):
	if not os.path.exists(path):
		print(f"There is no {name} at: {path}", file=sys.stderr)
		sys.exit(1)

def x_validate_directory__mutmut_1(name: str, path: str):
	if os.path.exists(path):
		print(f"There is no {name} at: {path}", file=sys.stderr)
		sys.exit(1)

def x_validate_directory__mutmut_2(name: str, path: str):
	if not os.path.exists(None):
		print(f"There is no {name} at: {path}", file=sys.stderr)
		sys.exit(1)

def x_validate_directory__mutmut_3(name: str, path: str):
	if not os.path.exists(path):
		print(None, file=sys.stderr)
		sys.exit(1)

def x_validate_directory__mutmut_4(name: str, path: str):
	if not os.path.exists(path):
		print(f"There is no {name} at: {path}", file=None)
		sys.exit(1)

def x_validate_directory__mutmut_5(name: str, path: str):
	if not os.path.exists(path):
		print(file=sys.stderr)
		sys.exit(1)

def x_validate_directory__mutmut_6(name: str, path: str):
	if not os.path.exists(path):
		print(f"There is no {name} at: {path}", )
		sys.exit(1)

def x_validate_directory__mutmut_7(name: str, path: str):
	if not os.path.exists(path):
		print(f"There is no {name} at: {path}", file=sys.stderr)
		sys.exit(None)

def x_validate_directory__mutmut_8(name: str, path: str):
	if not os.path.exists(path):
		print(f"There is no {name} at: {path}", file=sys.stderr)
		sys.exit(2)

x_validate_directory__mutmut_mutants : ClassVar[MutantDict] = { # type: ignore
'x_validate_directory__mutmut_1': x_validate_directory__mutmut_1, 
    'x_validate_directory__mutmut_2': x_validate_directory__mutmut_2, 
    'x_validate_directory__mutmut_3': x_validate_directory__mutmut_3, 
    'x_validate_directory__mutmut_4': x_validate_directory__mutmut_4, 
    'x_validate_directory__mutmut_5': x_validate_directory__mutmut_5, 
    'x_validate_directory__mutmut_6': x_validate_directory__mutmut_6, 
    'x_validate_directory__mutmut_7': x_validate_directory__mutmut_7, 
    'x_validate_directory__mutmut_8': x_validate_directory__mutmut_8
}
x_validate_directory__mutmut_orig.__name__ = 'x_validate_directory'

# Validate gisaid cli exists or error out
def validate_gisaid_installer(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	args = [submission_dir, organism, config_dict]# type: ignore
	kwargs = {}# type: ignore
	return _mutmut_trampoline(x_validate_gisaid_installer__mutmut_orig, x_validate_gisaid_installer__mutmut_mutants, args, kwargs, None)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_orig(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_1(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = None
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_2(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(None, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_3(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, None, organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_4(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", None)
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_5(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join("gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_6(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_7(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", )
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_8(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "XXgisaid_cliXX", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_9(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "GISAID_CLI", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_10(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower() - "CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_11(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.upper()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_12(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"XXCLIXX")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_13(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"cli")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_14(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = None
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_15(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(None, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_16(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, None, organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_17(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", None)
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_18(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join("gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_19(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_20(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", )
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_21(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "XXgisaid_cliXX", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_22(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "GISAID_CLI", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_23(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower() - "CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_24(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.upper()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_25(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"XXCLIXX")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_26(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"cli")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_27(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = None
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_28(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(None, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_29(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, None, organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_30(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", None, organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_31(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", None)
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_32(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join("gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_33(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_34(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_35(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", )
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_36(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "XXgisaid_cliXX", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_37(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "GISAID_CLI", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_38(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower() - "CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_39(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.upper()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_40(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"XXCLIXX", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_41(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"cli", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_42(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower() - "CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_43(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.upper()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_44(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"XXCLIXX")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_45(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"cli")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_46(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = None
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_47(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(None, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_48(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, None, organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_49(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", None, organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_50(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", None)
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_51(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join("gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_52(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_53(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_54(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", )
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_55(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "XXgisaid_cliXX", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_56(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "GISAID_CLI", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_57(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower() - "CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_58(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.upper()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_59(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"XXCLIXX", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_60(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"cli", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_61(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower() - "CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_62(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.upper()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_63(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"XXCLIXX")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_64(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"cli")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_65(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" or os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_66(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None or config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_67(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict or config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_68(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "XXCLI_PathXX" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_69(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "cli_path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_70(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_PATH" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_71(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" not in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_72(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["XXCLI_PathXX"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_73(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["cli_path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_74(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_PATH"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_75(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_76(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["XXCLI_PathXX"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_77(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["cli_path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_78(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_PATH"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_79(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() == "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_80(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "XXXX" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_81(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(None):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_82(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["XXCLI_PathXX"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_83(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["cli_path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_84(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_PATH"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_85(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["XXCLI_PathXX"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_86(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["cli_path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_87(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_PATH"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_88(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(None):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_89(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(None):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_90(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(None):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_91(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(None):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_92(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None or config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_93(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict or config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_94(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "XXCLI_PathXX" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_95(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "cli_path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_96(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_PATH" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_97(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" not in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_98(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["XXCLI_PathXX"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_99(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["cli_path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_100(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_PATH"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_101(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_102(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["XXCLI_PathXX"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_103(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["cli_path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_104(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_PATH"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_105(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() == "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_106(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "XXXX":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_107(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = None
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_108(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["XXCLI_PathXX"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_109(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["cli_path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_110(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_PATH"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_111(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(None, file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_112(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=None)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_113(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_114(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", )
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_115(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(None, file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_116(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=None)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_117(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_118(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", )
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_119(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(None, file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_120(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=None)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_121(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_122(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", )
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_123(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(None, file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_124(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=None)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_125(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(file=sys.stderr)
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_126(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", )
		sys.exit(1)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_127(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(None)

# Validate gisaid cli exists or error out
def x_validate_gisaid_installer__mutmut_128(submission_dir: str, organism: str, config_dict: dict[str, Any]) -> str:
	# /<submission_dir>/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_one = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI")
	# /seqsender/gisaid_cli/<organism>_CLI
	gisaid_cli_path_option_two = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI")
	# /<submission_dir>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_three = os.path.join(submission_dir, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# /seqsender>/gisaid_cli/<organism>_CLI/<organism>_CLI
	gisaid_cli_path_option_four = os.path.join(PROG_DIR, "gisaid_cli", organism.lower()+"CLI", organism.lower()+"CLI")
	# gisaid cli path provided by config file
	if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "" and os.path.isfile(config_dict["CLI_Path"].strip()):
		return config_dict["CLI_Path"].strip()
	elif os.path.isfile(gisaid_cli_path_option_one):
		return gisaid_cli_path_option_one
	elif os.path.isfile(gisaid_cli_path_option_two):
		return gisaid_cli_path_option_two
	elif os.path.isfile(gisaid_cli_path_option_three):
		return gisaid_cli_path_option_three
	elif os.path.isfile(gisaid_cli_path_option_four):
		return gisaid_cli_path_option_four
	else:
		if "CLI_Path" in config_dict and config_dict["CLI_Path"] is not None and config_dict["CLI_Path"].strip() != "":
			cli_path_error = config_dict["CLI_Path"]
			print(f"Error: There is not a GISAID CLI for {organism} provided via config file at: '{cli_path_error}'", file=sys.stderr)
		print(f"Error: There is not a GISAID CLI for {organism} located at: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		print(f"Download the GISAID CLI for {organism} from \"https://gisaid.org/\".", file=sys.stderr)
		print(f"Extract the zip file and place the CLI binary at either: '{gisaid_cli_path_option_one}' or '{gisaid_cli_path_option_two}'", file=sys.stderr)
		sys.exit(2)

x_validate_gisaid_installer__mutmut_mutants : ClassVar[MutantDict] = { # type: ignore
'x_validate_gisaid_installer__mutmut_1': x_validate_gisaid_installer__mutmut_1, 
    'x_validate_gisaid_installer__mutmut_2': x_validate_gisaid_installer__mutmut_2, 
    'x_validate_gisaid_installer__mutmut_3': x_validate_gisaid_installer__mutmut_3, 
    'x_validate_gisaid_installer__mutmut_4': x_validate_gisaid_installer__mutmut_4, 
    'x_validate_gisaid_installer__mutmut_5': x_validate_gisaid_installer__mutmut_5, 
    'x_validate_gisaid_installer__mutmut_6': x_validate_gisaid_installer__mutmut_6, 
    'x_validate_gisaid_installer__mutmut_7': x_validate_gisaid_installer__mutmut_7, 
    'x_validate_gisaid_installer__mutmut_8': x_validate_gisaid_installer__mutmut_8, 
    'x_validate_gisaid_installer__mutmut_9': x_validate_gisaid_installer__mutmut_9, 
    'x_validate_gisaid_installer__mutmut_10': x_validate_gisaid_installer__mutmut_10, 
    'x_validate_gisaid_installer__mutmut_11': x_validate_gisaid_installer__mutmut_11, 
    'x_validate_gisaid_installer__mutmut_12': x_validate_gisaid_installer__mutmut_12, 
    'x_validate_gisaid_installer__mutmut_13': x_validate_gisaid_installer__mutmut_13, 
    'x_validate_gisaid_installer__mutmut_14': x_validate_gisaid_installer__mutmut_14, 
    'x_validate_gisaid_installer__mutmut_15': x_validate_gisaid_installer__mutmut_15, 
    'x_validate_gisaid_installer__mutmut_16': x_validate_gisaid_installer__mutmut_16, 
    'x_validate_gisaid_installer__mutmut_17': x_validate_gisaid_installer__mutmut_17, 
    'x_validate_gisaid_installer__mutmut_18': x_validate_gisaid_installer__mutmut_18, 
    'x_validate_gisaid_installer__mutmut_19': x_validate_gisaid_installer__mutmut_19, 
    'x_validate_gisaid_installer__mutmut_20': x_validate_gisaid_installer__mutmut_20, 
    'x_validate_gisaid_installer__mutmut_21': x_validate_gisaid_installer__mutmut_21, 
    'x_validate_gisaid_installer__mutmut_22': x_validate_gisaid_installer__mutmut_22, 
    'x_validate_gisaid_installer__mutmut_23': x_validate_gisaid_installer__mutmut_23, 
    'x_validate_gisaid_installer__mutmut_24': x_validate_gisaid_installer__mutmut_24, 
    'x_validate_gisaid_installer__mutmut_25': x_validate_gisaid_installer__mutmut_25, 
    'x_validate_gisaid_installer__mutmut_26': x_validate_gisaid_installer__mutmut_26, 
    'x_validate_gisaid_installer__mutmut_27': x_validate_gisaid_installer__mutmut_27, 
    'x_validate_gisaid_installer__mutmut_28': x_validate_gisaid_installer__mutmut_28, 
    'x_validate_gisaid_installer__mutmut_29': x_validate_gisaid_installer__mutmut_29, 
    'x_validate_gisaid_installer__mutmut_30': x_validate_gisaid_installer__mutmut_30, 
    'x_validate_gisaid_installer__mutmut_31': x_validate_gisaid_installer__mutmut_31, 
    'x_validate_gisaid_installer__mutmut_32': x_validate_gisaid_installer__mutmut_32, 
    'x_validate_gisaid_installer__mutmut_33': x_validate_gisaid_installer__mutmut_33, 
    'x_validate_gisaid_installer__mutmut_34': x_validate_gisaid_installer__mutmut_34, 
    'x_validate_gisaid_installer__mutmut_35': x_validate_gisaid_installer__mutmut_35, 
    'x_validate_gisaid_installer__mutmut_36': x_validate_gisaid_installer__mutmut_36, 
    'x_validate_gisaid_installer__mutmut_37': x_validate_gisaid_installer__mutmut_37, 
    'x_validate_gisaid_installer__mutmut_38': x_validate_gisaid_installer__mutmut_38, 
    'x_validate_gisaid_installer__mutmut_39': x_validate_gisaid_installer__mutmut_39, 
    'x_validate_gisaid_installer__mutmut_40': x_validate_gisaid_installer__mutmut_40, 
    'x_validate_gisaid_installer__mutmut_41': x_validate_gisaid_installer__mutmut_41, 
    'x_validate_gisaid_installer__mutmut_42': x_validate_gisaid_installer__mutmut_42, 
    'x_validate_gisaid_installer__mutmut_43': x_validate_gisaid_installer__mutmut_43, 
    'x_validate_gisaid_installer__mutmut_44': x_validate_gisaid_installer__mutmut_44, 
    'x_validate_gisaid_installer__mutmut_45': x_validate_gisaid_installer__mutmut_45, 
    'x_validate_gisaid_installer__mutmut_46': x_validate_gisaid_installer__mutmut_46, 
    'x_validate_gisaid_installer__mutmut_47': x_validate_gisaid_installer__mutmut_47, 
    'x_validate_gisaid_installer__mutmut_48': x_validate_gisaid_installer__mutmut_48, 
    'x_validate_gisaid_installer__mutmut_49': x_validate_gisaid_installer__mutmut_49, 
    'x_validate_gisaid_installer__mutmut_50': x_validate_gisaid_installer__mutmut_50, 
    'x_validate_gisaid_installer__mutmut_51': x_validate_gisaid_installer__mutmut_51, 
    'x_validate_gisaid_installer__mutmut_52': x_validate_gisaid_installer__mutmut_52, 
    'x_validate_gisaid_installer__mutmut_53': x_validate_gisaid_installer__mutmut_53, 
    'x_validate_gisaid_installer__mutmut_54': x_validate_gisaid_installer__mutmut_54, 
    'x_validate_gisaid_installer__mutmut_55': x_validate_gisaid_installer__mutmut_55, 
    'x_validate_gisaid_installer__mutmut_56': x_validate_gisaid_installer__mutmut_56, 
    'x_validate_gisaid_installer__mutmut_57': x_validate_gisaid_installer__mutmut_57, 
    'x_validate_gisaid_installer__mutmut_58': x_validate_gisaid_installer__mutmut_58, 
    'x_validate_gisaid_installer__mutmut_59': x_validate_gisaid_installer__mutmut_59, 
    'x_validate_gisaid_installer__mutmut_60': x_validate_gisaid_installer__mutmut_60, 
    'x_validate_gisaid_installer__mutmut_61': x_validate_gisaid_installer__mutmut_61, 
    'x_validate_gisaid_installer__mutmut_62': x_validate_gisaid_installer__mutmut_62, 
    'x_validate_gisaid_installer__mutmut_63': x_validate_gisaid_installer__mutmut_63, 
    'x_validate_gisaid_installer__mutmut_64': x_validate_gisaid_installer__mutmut_64, 
    'x_validate_gisaid_installer__mutmut_65': x_validate_gisaid_installer__mutmut_65, 
    'x_validate_gisaid_installer__mutmut_66': x_validate_gisaid_installer__mutmut_66, 
    'x_validate_gisaid_installer__mutmut_67': x_validate_gisaid_installer__mutmut_67, 
    'x_validate_gisaid_installer__mutmut_68': x_validate_gisaid_installer__mutmut_68, 
    'x_validate_gisaid_installer__mutmut_69': x_validate_gisaid_installer__mutmut_69, 
    'x_validate_gisaid_installer__mutmut_70': x_validate_gisaid_installer__mutmut_70, 
    'x_validate_gisaid_installer__mutmut_71': x_validate_gisaid_installer__mutmut_71, 
    'x_validate_gisaid_installer__mutmut_72': x_validate_gisaid_installer__mutmut_72, 
    'x_validate_gisaid_installer__mutmut_73': x_validate_gisaid_installer__mutmut_73, 
    'x_validate_gisaid_installer__mutmut_74': x_validate_gisaid_installer__mutmut_74, 
    'x_validate_gisaid_installer__mutmut_75': x_validate_gisaid_installer__mutmut_75, 
    'x_validate_gisaid_installer__mutmut_76': x_validate_gisaid_installer__mutmut_76, 
    'x_validate_gisaid_installer__mutmut_77': x_validate_gisaid_installer__mutmut_77, 
    'x_validate_gisaid_installer__mutmut_78': x_validate_gisaid_installer__mutmut_78, 
    'x_validate_gisaid_installer__mutmut_79': x_validate_gisaid_installer__mutmut_79, 
    'x_validate_gisaid_installer__mutmut_80': x_validate_gisaid_installer__mutmut_80, 
    'x_validate_gisaid_installer__mutmut_81': x_validate_gisaid_installer__mutmut_81, 
    'x_validate_gisaid_installer__mutmut_82': x_validate_gisaid_installer__mutmut_82, 
    'x_validate_gisaid_installer__mutmut_83': x_validate_gisaid_installer__mutmut_83, 
    'x_validate_gisaid_installer__mutmut_84': x_validate_gisaid_installer__mutmut_84, 
    'x_validate_gisaid_installer__mutmut_85': x_validate_gisaid_installer__mutmut_85, 
    'x_validate_gisaid_installer__mutmut_86': x_validate_gisaid_installer__mutmut_86, 
    'x_validate_gisaid_installer__mutmut_87': x_validate_gisaid_installer__mutmut_87, 
    'x_validate_gisaid_installer__mutmut_88': x_validate_gisaid_installer__mutmut_88, 
    'x_validate_gisaid_installer__mutmut_89': x_validate_gisaid_installer__mutmut_89, 
    'x_validate_gisaid_installer__mutmut_90': x_validate_gisaid_installer__mutmut_90, 
    'x_validate_gisaid_installer__mutmut_91': x_validate_gisaid_installer__mutmut_91, 
    'x_validate_gisaid_installer__mutmut_92': x_validate_gisaid_installer__mutmut_92, 
    'x_validate_gisaid_installer__mutmut_93': x_validate_gisaid_installer__mutmut_93, 
    'x_validate_gisaid_installer__mutmut_94': x_validate_gisaid_installer__mutmut_94, 
    'x_validate_gisaid_installer__mutmut_95': x_validate_gisaid_installer__mutmut_95, 
    'x_validate_gisaid_installer__mutmut_96': x_validate_gisaid_installer__mutmut_96, 
    'x_validate_gisaid_installer__mutmut_97': x_validate_gisaid_installer__mutmut_97, 
    'x_validate_gisaid_installer__mutmut_98': x_validate_gisaid_installer__mutmut_98, 
    'x_validate_gisaid_installer__mutmut_99': x_validate_gisaid_installer__mutmut_99, 
    'x_validate_gisaid_installer__mutmut_100': x_validate_gisaid_installer__mutmut_100, 
    'x_validate_gisaid_installer__mutmut_101': x_validate_gisaid_installer__mutmut_101, 
    'x_validate_gisaid_installer__mutmut_102': x_validate_gisaid_installer__mutmut_102, 
    'x_validate_gisaid_installer__mutmut_103': x_validate_gisaid_installer__mutmut_103, 
    'x_validate_gisaid_installer__mutmut_104': x_validate_gisaid_installer__mutmut_104, 
    'x_validate_gisaid_installer__mutmut_105': x_validate_gisaid_installer__mutmut_105, 
    'x_validate_gisaid_installer__mutmut_106': x_validate_gisaid_installer__mutmut_106, 
    'x_validate_gisaid_installer__mutmut_107': x_validate_gisaid_installer__mutmut_107, 
    'x_validate_gisaid_installer__mutmut_108': x_validate_gisaid_installer__mutmut_108, 
    'x_validate_gisaid_installer__mutmut_109': x_validate_gisaid_installer__mutmut_109, 
    'x_validate_gisaid_installer__mutmut_110': x_validate_gisaid_installer__mutmut_110, 
    'x_validate_gisaid_installer__mutmut_111': x_validate_gisaid_installer__mutmut_111, 
    'x_validate_gisaid_installer__mutmut_112': x_validate_gisaid_installer__mutmut_112, 
    'x_validate_gisaid_installer__mutmut_113': x_validate_gisaid_installer__mutmut_113, 
    'x_validate_gisaid_installer__mutmut_114': x_validate_gisaid_installer__mutmut_114, 
    'x_validate_gisaid_installer__mutmut_115': x_validate_gisaid_installer__mutmut_115, 
    'x_validate_gisaid_installer__mutmut_116': x_validate_gisaid_installer__mutmut_116, 
    'x_validate_gisaid_installer__mutmut_117': x_validate_gisaid_installer__mutmut_117, 
    'x_validate_gisaid_installer__mutmut_118': x_validate_gisaid_installer__mutmut_118, 
    'x_validate_gisaid_installer__mutmut_119': x_validate_gisaid_installer__mutmut_119, 
    'x_validate_gisaid_installer__mutmut_120': x_validate_gisaid_installer__mutmut_120, 
    'x_validate_gisaid_installer__mutmut_121': x_validate_gisaid_installer__mutmut_121, 
    'x_validate_gisaid_installer__mutmut_122': x_validate_gisaid_installer__mutmut_122, 
    'x_validate_gisaid_installer__mutmut_123': x_validate_gisaid_installer__mutmut_123, 
    'x_validate_gisaid_installer__mutmut_124': x_validate_gisaid_installer__mutmut_124, 
    'x_validate_gisaid_installer__mutmut_125': x_validate_gisaid_installer__mutmut_125, 
    'x_validate_gisaid_installer__mutmut_126': x_validate_gisaid_installer__mutmut_126, 
    'x_validate_gisaid_installer__mutmut_127': x_validate_gisaid_installer__mutmut_127, 
    'x_validate_gisaid_installer__mutmut_128': x_validate_gisaid_installer__mutmut_128
}
x_validate_gisaid_installer__mutmut_orig.__name__ = 'x_validate_gisaid_installer'

# Create directory and don't error out if it already exists
def create_directory(path: str):
	args = [path]# type: ignore
	kwargs = {}# type: ignore
	return _mutmut_trampoline(x_create_directory__mutmut_orig, x_create_directory__mutmut_mutants, args, kwargs, None)

# Create directory and don't error out if it already exists
def x_create_directory__mutmut_orig(path: str):
	os.makedirs(path, exist_ok = True)

# Create directory and don't error out if it already exists
def x_create_directory__mutmut_1(path: str):
	os.makedirs(None, exist_ok = True)

# Create directory and don't error out if it already exists
def x_create_directory__mutmut_2(path: str):
	os.makedirs(path, exist_ok = None)

# Create directory and don't error out if it already exists
def x_create_directory__mutmut_3(path: str):
	os.makedirs(exist_ok = True)

# Create directory and don't error out if it already exists
def x_create_directory__mutmut_4(path: str):
	os.makedirs(path, )

# Create directory and don't error out if it already exists
def x_create_directory__mutmut_5(path: str):
	os.makedirs(path, exist_ok = False)

x_create_directory__mutmut_mutants : ClassVar[MutantDict] = { # type: ignore
'x_create_directory__mutmut_1': x_create_directory__mutmut_1, 
    'x_create_directory__mutmut_2': x_create_directory__mutmut_2, 
    'x_create_directory__mutmut_3': x_create_directory__mutmut_3, 
    'x_create_directory__mutmut_4': x_create_directory__mutmut_4, 
    'x_create_directory__mutmut_5': x_create_directory__mutmut_5
}
x_create_directory__mutmut_orig.__name__ = 'x_create_directory'

# Load yaml file or error out
def load_yaml(yaml_type: str, yaml_path: str):
	args = [yaml_type, yaml_path]# type: ignore
	kwargs = {}# type: ignore
	return _mutmut_trampoline(x_load_yaml__mutmut_orig, x_load_yaml__mutmut_mutants, args, kwargs, None)

# Load yaml file or error out
def x_load_yaml__mutmut_orig(yaml_type: str, yaml_path: str):
	with open(yaml_path, "r") as file:
		try:
			config_dict = yaml.load(file, Loader = yaml.FullLoader)
		except:
			print(f"Error: {yaml_type} is incorrect. File must be a valid yaml format.", file=sys.stderr)
			sys.exit(1)
	return config_dict

# Load yaml file or error out
def x_load_yaml__mutmut_1(yaml_type: str, yaml_path: str):
	with open(None, "r") as file:
		try:
			config_dict = yaml.load(file, Loader = yaml.FullLoader)
		except:
			print(f"Error: {yaml_type} is incorrect. File must be a valid yaml format.", file=sys.stderr)
			sys.exit(1)
	return config_dict

# Load yaml file or error out
def x_load_yaml__mutmut_2(yaml_type: str, yaml_path: str):
	with open(yaml_path, None) as file:
		try:
			config_dict = yaml.load(file, Loader = yaml.FullLoader)
		except:
			print(f"Error: {yaml_type} is incorrect. File must be a valid yaml format.", file=sys.stderr)
			sys.exit(1)
	return config_dict

# Load yaml file or error out
def x_load_yaml__mutmut_3(yaml_type: str, yaml_path: str):
	with open("r") as file:
		try:
			config_dict = yaml.load(file, Loader = yaml.FullLoader)
		except:
			print(f"Error: {yaml_type} is incorrect. File must be a valid yaml format.", file=sys.stderr)
			sys.exit(1)
	return config_dict

# Load yaml file or error out
def x_load_yaml__mutmut_4(yaml_type: str, yaml_path: str):
	with open(yaml_path, ) as file:
		try:
			config_dict = yaml.load(file, Loader = yaml.FullLoader)
		except:
			print(f"Error: {yaml_type} is incorrect. File must be a valid yaml format.", file=sys.stderr)
			sys.exit(1)
	return config_dict

# Load yaml file or error out
def x_load_yaml__mutmut_5(yaml_type: str, yaml_path: str):
	with open(yaml_path, "XXrXX") as file:
		try:
			config_dict = yaml.load(file, Loader = yaml.FullLoader)
		except:
			print(f"Error: {yaml_type} is incorrect. File must be a valid yaml format.", file=sys.stderr)
			sys.exit(1)
	return config_dict

# Load yaml file or error out
def x_load_yaml__mutmut_6(yaml_type: str, yaml_path: str):
	with open(yaml_path, "R") as file:
		try:
			config_dict = yaml.load(file, Loader = yaml.FullLoader)
		except:
			print(f"Error: {yaml_type} is incorrect. File must be a valid yaml format.", file=sys.stderr)
			sys.exit(1)
	return config_dict

# Load yaml file or error out
def x_load_yaml__mutmut_7(yaml_type: str, yaml_path: str):
	with open(yaml_path, "r") as file:
		try:
			config_dict = None
		except:
			print(f"Error: {yaml_type} is incorrect. File must be a valid yaml format.", file=sys.stderr)
			sys.exit(1)
	return config_dict

# Load yaml file or error out
def x_load_yaml__mutmut_8(yaml_type: str, yaml_path: str):
	with open(yaml_path, "r") as file:
		try:
			config_dict = yaml.load(None, Loader = yaml.FullLoader)
		except:
			print(f"Error: {yaml_type} is incorrect. File must be a valid yaml format.", file=sys.stderr)
			sys.exit(1)
	return config_dict

# Load yaml file or error out
def x_load_yaml__mutmut_9(yaml_type: str, yaml_path: str):
	with open(yaml_path, "r") as file:
		try:
			config_dict = yaml.load(file, Loader = None)
		except:
			print(f"Error: {yaml_type} is incorrect. File must be a valid yaml format.", file=sys.stderr)
			sys.exit(1)
	return config_dict

# Load yaml file or error out
def x_load_yaml__mutmut_10(yaml_type: str, yaml_path: str):
	with open(yaml_path, "r") as file:
		try:
			config_dict = yaml.load(Loader = yaml.FullLoader)
		except:
			print(f"Error: {yaml_type} is incorrect. File must be a valid yaml format.", file=sys.stderr)
			sys.exit(1)
	return config_dict

# Load yaml file or error out
def x_load_yaml__mutmut_11(yaml_type: str, yaml_path: str):
	with open(yaml_path, "r") as file:
		try:
			config_dict = yaml.load(file, )
		except:
			print(f"Error: {yaml_type} is incorrect. File must be a valid yaml format.", file=sys.stderr)
			sys.exit(1)
	return config_dict

# Load yaml file or error out
def x_load_yaml__mutmut_12(yaml_type: str, yaml_path: str):
	with open(yaml_path, "r") as file:
		try:
			config_dict = yaml.load(file, Loader = yaml.FullLoader)
		except:
			print(None, file=sys.stderr)
			sys.exit(1)
	return config_dict

# Load yaml file or error out
def x_load_yaml__mutmut_13(yaml_type: str, yaml_path: str):
	with open(yaml_path, "r") as file:
		try:
			config_dict = yaml.load(file, Loader = yaml.FullLoader)
		except:
			print(f"Error: {yaml_type} is incorrect. File must be a valid yaml format.", file=None)
			sys.exit(1)
	return config_dict

# Load yaml file or error out
def x_load_yaml__mutmut_14(yaml_type: str, yaml_path: str):
	with open(yaml_path, "r") as file:
		try:
			config_dict = yaml.load(file, Loader = yaml.FullLoader)
		except:
			print(file=sys.stderr)
			sys.exit(1)
	return config_dict

# Load yaml file or error out
def x_load_yaml__mutmut_15(yaml_type: str, yaml_path: str):
	with open(yaml_path, "r") as file:
		try:
			config_dict = yaml.load(file, Loader = yaml.FullLoader)
		except:
			print(f"Error: {yaml_type} is incorrect. File must be a valid yaml format.", )
			sys.exit(1)
	return config_dict

# Load yaml file or error out
def x_load_yaml__mutmut_16(yaml_type: str, yaml_path: str):
	with open(yaml_path, "r") as file:
		try:
			config_dict = yaml.load(file, Loader = yaml.FullLoader)
		except:
			print(f"Error: {yaml_type} is incorrect. File must be a valid yaml format.", file=sys.stderr)
			sys.exit(None)
	return config_dict

# Load yaml file or error out
def x_load_yaml__mutmut_17(yaml_type: str, yaml_path: str):
	with open(yaml_path, "r") as file:
		try:
			config_dict = yaml.load(file, Loader = yaml.FullLoader)
		except:
			print(f"Error: {yaml_type} is incorrect. File must be a valid yaml format.", file=sys.stderr)
			sys.exit(2)
	return config_dict

x_load_yaml__mutmut_mutants : ClassVar[MutantDict] = { # type: ignore
'x_load_yaml__mutmut_1': x_load_yaml__mutmut_1, 
    'x_load_yaml__mutmut_2': x_load_yaml__mutmut_2, 
    'x_load_yaml__mutmut_3': x_load_yaml__mutmut_3, 
    'x_load_yaml__mutmut_4': x_load_yaml__mutmut_4, 
    'x_load_yaml__mutmut_5': x_load_yaml__mutmut_5, 
    'x_load_yaml__mutmut_6': x_load_yaml__mutmut_6, 
    'x_load_yaml__mutmut_7': x_load_yaml__mutmut_7, 
    'x_load_yaml__mutmut_8': x_load_yaml__mutmut_8, 
    'x_load_yaml__mutmut_9': x_load_yaml__mutmut_9, 
    'x_load_yaml__mutmut_10': x_load_yaml__mutmut_10, 
    'x_load_yaml__mutmut_11': x_load_yaml__mutmut_11, 
    'x_load_yaml__mutmut_12': x_load_yaml__mutmut_12, 
    'x_load_yaml__mutmut_13': x_load_yaml__mutmut_13, 
    'x_load_yaml__mutmut_14': x_load_yaml__mutmut_14, 
    'x_load_yaml__mutmut_15': x_load_yaml__mutmut_15, 
    'x_load_yaml__mutmut_16': x_load_yaml__mutmut_16, 
    'x_load_yaml__mutmut_17': x_load_yaml__mutmut_17
}
x_load_yaml__mutmut_orig.__name__ = 'x_load_yaml'

# Is a entire pandas row made of just whitespace, empty strings, or None
def is_row_empty(row: pd.Series) -> bool:
	args = [row]# type: ignore
	kwargs = {}# type: ignore
	return _mutmut_trampoline(x_is_row_empty__mutmut_orig, x_is_row_empty__mutmut_mutants, args, kwargs, None)

# Is a entire pandas row made of just whitespace, empty strings, or None
def x_is_row_empty__mutmut_orig(row: pd.Series) -> bool:
	return all(cell is None or (isinstance(cell, str) and cell.strip() == "") for cell in row)

# Is a entire pandas row made of just whitespace, empty strings, or None
def x_is_row_empty__mutmut_1(row: pd.Series) -> bool:
	return all(None)

# Is a entire pandas row made of just whitespace, empty strings, or None
def x_is_row_empty__mutmut_2(row: pd.Series) -> bool:
	return all(cell is None and (isinstance(cell, str) and cell.strip() == "") for cell in row)

# Is a entire pandas row made of just whitespace, empty strings, or None
def x_is_row_empty__mutmut_3(row: pd.Series) -> bool:
	return all(cell is not None or (isinstance(cell, str) and cell.strip() == "") for cell in row)

# Is a entire pandas row made of just whitespace, empty strings, or None
def x_is_row_empty__mutmut_4(row: pd.Series) -> bool:
	return all(cell is None or (isinstance(cell, str) or cell.strip() == "") for cell in row)

# Is a entire pandas row made of just whitespace, empty strings, or None
def x_is_row_empty__mutmut_5(row: pd.Series) -> bool:
	return all(cell is None or (isinstance(cell, str) and cell.strip() != "") for cell in row)

# Is a entire pandas row made of just whitespace, empty strings, or None
def x_is_row_empty__mutmut_6(row: pd.Series) -> bool:
	return all(cell is None or (isinstance(cell, str) and cell.strip() == "XXXX") for cell in row)

x_is_row_empty__mutmut_mutants : ClassVar[MutantDict] = { # type: ignore
'x_is_row_empty__mutmut_1': x_is_row_empty__mutmut_1, 
    'x_is_row_empty__mutmut_2': x_is_row_empty__mutmut_2, 
    'x_is_row_empty__mutmut_3': x_is_row_empty__mutmut_3, 
    'x_is_row_empty__mutmut_4': x_is_row_empty__mutmut_4, 
    'x_is_row_empty__mutmut_5': x_is_row_empty__mutmut_5, 
    'x_is_row_empty__mutmut_6': x_is_row_empty__mutmut_6
}
x_is_row_empty__mutmut_orig.__name__ = 'x_is_row_empty'

# Load csv into pandas df and clean then return
def load_csv(file_path: str, sep: str = ",") -> pd.DataFrame:
	args = [file_path, sep]# type: ignore
	kwargs = {}# type: ignore
	return _mutmut_trampoline(x_load_csv__mutmut_orig, x_load_csv__mutmut_mutants, args, kwargs, None)

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_orig(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_1(file_path: str, sep: str = "XX,XX") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_2(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = None
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_3(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(None, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_4(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = None, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_5(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = None, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_6(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = None, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_7(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = None, encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_8(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = None, index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_9(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = None, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_10(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = None)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_11(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_12(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_13(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_14(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_15(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_16(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_17(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_18(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, )
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_19(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 1, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_20(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "XXpythonXX", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_21(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "PYTHON", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_22(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "XXutf-8XX", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_23(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "UTF-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_24(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = True, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_25(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = True)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_26(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = None # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_27(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(None, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_28(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = None) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_29(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_30(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, ) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_31(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_32(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) or is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_33(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(None) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_34(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 2) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_35(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = None
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_36(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = None)
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_37(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "XXallXX")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_38(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "ALL")
	# Remove extra spaces from column names that could cause issues
	df.columns = df.columns.str.strip()
	return df

# Load csv into pandas df and clean then return
def x_load_csv__mutmut_39(file_path: str, sep: str = ",") -> pd.DataFrame:
	df = pd.read_csv(file_path, header = 0, dtype = str, sep = sep, engine = "python", encoding = "utf-8", index_col = False, na_filter = False)
	# Replace empty strings and strings entirely made of whitespaces with nan
	df = df.apply(lambda row: None if isinstance(row, pd.Series) and is_row_empty(row) else row, axis = 1) # type: ignore
	# Drop rows if entirely empty
	df = df.dropna(how = "all")
	# Remove extra spaces from column names that could cause issues
	df.columns = None
	return df

x_load_csv__mutmut_mutants : ClassVar[MutantDict] = { # type: ignore
'x_load_csv__mutmut_1': x_load_csv__mutmut_1, 
    'x_load_csv__mutmut_2': x_load_csv__mutmut_2, 
    'x_load_csv__mutmut_3': x_load_csv__mutmut_3, 
    'x_load_csv__mutmut_4': x_load_csv__mutmut_4, 
    'x_load_csv__mutmut_5': x_load_csv__mutmut_5, 
    'x_load_csv__mutmut_6': x_load_csv__mutmut_6, 
    'x_load_csv__mutmut_7': x_load_csv__mutmut_7, 
    'x_load_csv__mutmut_8': x_load_csv__mutmut_8, 
    'x_load_csv__mutmut_9': x_load_csv__mutmut_9, 
    'x_load_csv__mutmut_10': x_load_csv__mutmut_10, 
    'x_load_csv__mutmut_11': x_load_csv__mutmut_11, 
    'x_load_csv__mutmut_12': x_load_csv__mutmut_12, 
    'x_load_csv__mutmut_13': x_load_csv__mutmut_13, 
    'x_load_csv__mutmut_14': x_load_csv__mutmut_14, 
    'x_load_csv__mutmut_15': x_load_csv__mutmut_15, 
    'x_load_csv__mutmut_16': x_load_csv__mutmut_16, 
    'x_load_csv__mutmut_17': x_load_csv__mutmut_17, 
    'x_load_csv__mutmut_18': x_load_csv__mutmut_18, 
    'x_load_csv__mutmut_19': x_load_csv__mutmut_19, 
    'x_load_csv__mutmut_20': x_load_csv__mutmut_20, 
    'x_load_csv__mutmut_21': x_load_csv__mutmut_21, 
    'x_load_csv__mutmut_22': x_load_csv__mutmut_22, 
    'x_load_csv__mutmut_23': x_load_csv__mutmut_23, 
    'x_load_csv__mutmut_24': x_load_csv__mutmut_24, 
    'x_load_csv__mutmut_25': x_load_csv__mutmut_25, 
    'x_load_csv__mutmut_26': x_load_csv__mutmut_26, 
    'x_load_csv__mutmut_27': x_load_csv__mutmut_27, 
    'x_load_csv__mutmut_28': x_load_csv__mutmut_28, 
    'x_load_csv__mutmut_29': x_load_csv__mutmut_29, 
    'x_load_csv__mutmut_30': x_load_csv__mutmut_30, 
    'x_load_csv__mutmut_31': x_load_csv__mutmut_31, 
    'x_load_csv__mutmut_32': x_load_csv__mutmut_32, 
    'x_load_csv__mutmut_33': x_load_csv__mutmut_33, 
    'x_load_csv__mutmut_34': x_load_csv__mutmut_34, 
    'x_load_csv__mutmut_35': x_load_csv__mutmut_35, 
    'x_load_csv__mutmut_36': x_load_csv__mutmut_36, 
    'x_load_csv__mutmut_37': x_load_csv__mutmut_37, 
    'x_load_csv__mutmut_38': x_load_csv__mutmut_38, 
    'x_load_csv__mutmut_39': x_load_csv__mutmut_39
}
x_load_csv__mutmut_orig.__name__ = 'x_load_csv'

# Load fasta file into pandas df and clean then return
def load_fasta_file(fasta_file: str) -> pd.DataFrame:
	args = [fasta_file]# type: ignore
	kwargs = {}# type: ignore
	return _mutmut_trampoline(x_load_fasta_file__mutmut_orig, x_load_fasta_file__mutmut_mutants, args, kwargs, None)

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_orig(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "r") as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_1(fasta_file: str) -> pd.DataFrame:
	fasta_dict = None
	with open(fasta_file, "r") as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_2(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(None, "r") as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_3(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, None) as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_4(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open("r") as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_5(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, ) as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_6(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "XXrXX") as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_7(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "R") as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_8(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "r") as fsa:
		records = None
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_9(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "r") as fsa:
		records = SeqIO.parse(None, "fasta")
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_10(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "r") as fsa:
		records = SeqIO.parse(fsa, None)
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_11(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "r") as fsa:
		records = SeqIO.parse("fasta")
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_12(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "r") as fsa:
		records = SeqIO.parse(fsa, )
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_13(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "r") as fsa:
		records = SeqIO.parse(fsa, "XXfastaXX")
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_14(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "r") as fsa:
		records = SeqIO.parse(fsa, "FASTA")
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_15(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "r") as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			fasta_dict.append(None)
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_16(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "r") as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			fasta_dict.append({"XXfasta_name_origXX":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_17(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "r") as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			fasta_dict.append({"FASTA_NAME_ORIG":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_18(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "r") as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "XXfasta_sequence_origXX":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_19(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "r") as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "FASTA_SEQUENCE_ORIG":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_20(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "r") as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "XXfasta_description_origXX":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_21(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "r") as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "FASTA_DESCRIPTION_ORIG":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_22(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "r") as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = None
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_23(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "r") as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(None)
	fasta_df = fasta_df.dropna(how = "all")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_24(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "r") as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = None
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_25(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "r") as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = None)
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_26(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "r") as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "XXallXX")
	return fasta_df

# Load fasta file into pandas df and clean then return
def x_load_fasta_file__mutmut_27(fasta_file: str) -> pd.DataFrame:
	fasta_dict = []
	with open(fasta_file, "r") as fsa:
		records = SeqIO.parse(fsa, "fasta")
		for record in records:
			fasta_dict.append({"fasta_name_orig":record.id, "fasta_sequence_orig":record.seq, "fasta_description_orig":record.description})
	fasta_df = pd.DataFrame(fasta_dict)
	fasta_df = fasta_df.dropna(how = "ALL")
	return fasta_df

x_load_fasta_file__mutmut_mutants : ClassVar[MutantDict] = { # type: ignore
'x_load_fasta_file__mutmut_1': x_load_fasta_file__mutmut_1, 
    'x_load_fasta_file__mutmut_2': x_load_fasta_file__mutmut_2, 
    'x_load_fasta_file__mutmut_3': x_load_fasta_file__mutmut_3, 
    'x_load_fasta_file__mutmut_4': x_load_fasta_file__mutmut_4, 
    'x_load_fasta_file__mutmut_5': x_load_fasta_file__mutmut_5, 
    'x_load_fasta_file__mutmut_6': x_load_fasta_file__mutmut_6, 
    'x_load_fasta_file__mutmut_7': x_load_fasta_file__mutmut_7, 
    'x_load_fasta_file__mutmut_8': x_load_fasta_file__mutmut_8, 
    'x_load_fasta_file__mutmut_9': x_load_fasta_file__mutmut_9, 
    'x_load_fasta_file__mutmut_10': x_load_fasta_file__mutmut_10, 
    'x_load_fasta_file__mutmut_11': x_load_fasta_file__mutmut_11, 
    'x_load_fasta_file__mutmut_12': x_load_fasta_file__mutmut_12, 
    'x_load_fasta_file__mutmut_13': x_load_fasta_file__mutmut_13, 
    'x_load_fasta_file__mutmut_14': x_load_fasta_file__mutmut_14, 
    'x_load_fasta_file__mutmut_15': x_load_fasta_file__mutmut_15, 
    'x_load_fasta_file__mutmut_16': x_load_fasta_file__mutmut_16, 
    'x_load_fasta_file__mutmut_17': x_load_fasta_file__mutmut_17, 
    'x_load_fasta_file__mutmut_18': x_load_fasta_file__mutmut_18, 
    'x_load_fasta_file__mutmut_19': x_load_fasta_file__mutmut_19, 
    'x_load_fasta_file__mutmut_20': x_load_fasta_file__mutmut_20, 
    'x_load_fasta_file__mutmut_21': x_load_fasta_file__mutmut_21, 
    'x_load_fasta_file__mutmut_22': x_load_fasta_file__mutmut_22, 
    'x_load_fasta_file__mutmut_23': x_load_fasta_file__mutmut_23, 
    'x_load_fasta_file__mutmut_24': x_load_fasta_file__mutmut_24, 
    'x_load_fasta_file__mutmut_25': x_load_fasta_file__mutmut_25, 
    'x_load_fasta_file__mutmut_26': x_load_fasta_file__mutmut_26, 
    'x_load_fasta_file__mutmut_27': x_load_fasta_file__mutmut_27
}
x_load_fasta_file__mutmut_orig.__name__ = 'x_load_fasta_file'

# Save submission xml
def save_xml(submission_xml: bytes, submission_dir: str) -> None:
	args = [submission_xml, submission_dir]# type: ignore
	kwargs = {}# type: ignore
	return _mutmut_trampoline(x_save_xml__mutmut_orig, x_save_xml__mutmut_mutants, args, kwargs, None)

# Save submission xml
def x_save_xml__mutmut_orig(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_1(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(None, "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_2(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), None) as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_3(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open("wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_4(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), ) as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_5(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(None, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_6(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, None), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_7(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join("submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_8(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, ), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_9(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "XXsubmission.xmlXX"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_10(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "SUBMISSION.XML"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_11(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "XXwbXX") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_12(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "WB") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_13(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(None)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_14(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(None, file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_15(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=None)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_16(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_17(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", )
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_18(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(None, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_19(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=None)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_20(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_21(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, )
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_22(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(None)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_23(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(2)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_24(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(None, file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_25(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=None)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_26(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_27(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", )
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_28(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(None, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_29(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=None)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_30(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_31(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, )
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_32(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(None)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_33(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(2)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_34(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_35(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(None):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_36(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(None, "submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_37(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, None)):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_38(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join("submission.xml")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_39(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, )):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_40(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "XXsubmission.xmlXX")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_41(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "SUBMISSION.XML")):
		time.sleep(10)

# Save submission xml
def x_save_xml__mutmut_42(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(None)

# Save submission xml
def x_save_xml__mutmut_43(submission_xml: bytes, submission_dir: str) -> None:
	# Save string as submission.xml
	try:
		with open(os.path.join(submission_dir, "submission.xml"), "wb") as file:
			file.write(submission_xml)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	# Waiting for the xml file to write
	while not os.path.exists(os.path.join(submission_dir, "submission.xml")):
		time.sleep(11)

x_save_xml__mutmut_mutants : ClassVar[MutantDict] = { # type: ignore
'x_save_xml__mutmut_1': x_save_xml__mutmut_1, 
    'x_save_xml__mutmut_2': x_save_xml__mutmut_2, 
    'x_save_xml__mutmut_3': x_save_xml__mutmut_3, 
    'x_save_xml__mutmut_4': x_save_xml__mutmut_4, 
    'x_save_xml__mutmut_5': x_save_xml__mutmut_5, 
    'x_save_xml__mutmut_6': x_save_xml__mutmut_6, 
    'x_save_xml__mutmut_7': x_save_xml__mutmut_7, 
    'x_save_xml__mutmut_8': x_save_xml__mutmut_8, 
    'x_save_xml__mutmut_9': x_save_xml__mutmut_9, 
    'x_save_xml__mutmut_10': x_save_xml__mutmut_10, 
    'x_save_xml__mutmut_11': x_save_xml__mutmut_11, 
    'x_save_xml__mutmut_12': x_save_xml__mutmut_12, 
    'x_save_xml__mutmut_13': x_save_xml__mutmut_13, 
    'x_save_xml__mutmut_14': x_save_xml__mutmut_14, 
    'x_save_xml__mutmut_15': x_save_xml__mutmut_15, 
    'x_save_xml__mutmut_16': x_save_xml__mutmut_16, 
    'x_save_xml__mutmut_17': x_save_xml__mutmut_17, 
    'x_save_xml__mutmut_18': x_save_xml__mutmut_18, 
    'x_save_xml__mutmut_19': x_save_xml__mutmut_19, 
    'x_save_xml__mutmut_20': x_save_xml__mutmut_20, 
    'x_save_xml__mutmut_21': x_save_xml__mutmut_21, 
    'x_save_xml__mutmut_22': x_save_xml__mutmut_22, 
    'x_save_xml__mutmut_23': x_save_xml__mutmut_23, 
    'x_save_xml__mutmut_24': x_save_xml__mutmut_24, 
    'x_save_xml__mutmut_25': x_save_xml__mutmut_25, 
    'x_save_xml__mutmut_26': x_save_xml__mutmut_26, 
    'x_save_xml__mutmut_27': x_save_xml__mutmut_27, 
    'x_save_xml__mutmut_28': x_save_xml__mutmut_28, 
    'x_save_xml__mutmut_29': x_save_xml__mutmut_29, 
    'x_save_xml__mutmut_30': x_save_xml__mutmut_30, 
    'x_save_xml__mutmut_31': x_save_xml__mutmut_31, 
    'x_save_xml__mutmut_32': x_save_xml__mutmut_32, 
    'x_save_xml__mutmut_33': x_save_xml__mutmut_33, 
    'x_save_xml__mutmut_34': x_save_xml__mutmut_34, 
    'x_save_xml__mutmut_35': x_save_xml__mutmut_35, 
    'x_save_xml__mutmut_36': x_save_xml__mutmut_36, 
    'x_save_xml__mutmut_37': x_save_xml__mutmut_37, 
    'x_save_xml__mutmut_38': x_save_xml__mutmut_38, 
    'x_save_xml__mutmut_39': x_save_xml__mutmut_39, 
    'x_save_xml__mutmut_40': x_save_xml__mutmut_40, 
    'x_save_xml__mutmut_41': x_save_xml__mutmut_41, 
    'x_save_xml__mutmut_42': x_save_xml__mutmut_42, 
    'x_save_xml__mutmut_43': x_save_xml__mutmut_43
}
x_save_xml__mutmut_orig.__name__ = 'x_save_xml'

# Save pandas df to csv
def save_csv(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	args = [df, file_path, file_name, sep]# type: ignore
	kwargs = {}# type: ignore
	return _mutmut_trampoline(x_save_csv__mutmut_orig, x_save_csv__mutmut_mutants, args, kwargs, None)

# Save pandas df to csv
def x_save_csv__mutmut_orig(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_1(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = "XX,XX") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_2(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = None
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_3(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(None, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_4(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, None)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_5(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_6(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, )
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_7(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(None, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_8(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = None, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_9(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = None, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_10(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = None)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_11(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_12(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_13(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_14(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, )
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_15(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = False, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_16(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = True, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_17(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(None, file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_18(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=None)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_19(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_20(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", )
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_21(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(None, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_22(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=None)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_23(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_24(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, )
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_25(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(None)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_26(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(2)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_27(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(None, file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_28(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=None)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_29(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_30(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", )
		print(e, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_31(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(None, file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_32(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=None)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_33(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(file=sys.stderr)
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_34(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, )
		sys.exit(1)

# Save pandas df to csv
def x_save_csv__mutmut_35(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(None)

# Save pandas df to csv
def x_save_csv__mutmut_36(df: pd.DataFrame, file_path: str, file_name: Optional[str] = None, sep: str = ",") -> None:
	if file_name:
		file_path = os.path.join(file_path, file_name)
	try:
		df.to_csv(file_path, header = True, index = False, sep = sep)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save '{file_name}' to path: {file_path}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(2)

x_save_csv__mutmut_mutants : ClassVar[MutantDict] = { # type: ignore
'x_save_csv__mutmut_1': x_save_csv__mutmut_1, 
    'x_save_csv__mutmut_2': x_save_csv__mutmut_2, 
    'x_save_csv__mutmut_3': x_save_csv__mutmut_3, 
    'x_save_csv__mutmut_4': x_save_csv__mutmut_4, 
    'x_save_csv__mutmut_5': x_save_csv__mutmut_5, 
    'x_save_csv__mutmut_6': x_save_csv__mutmut_6, 
    'x_save_csv__mutmut_7': x_save_csv__mutmut_7, 
    'x_save_csv__mutmut_8': x_save_csv__mutmut_8, 
    'x_save_csv__mutmut_9': x_save_csv__mutmut_9, 
    'x_save_csv__mutmut_10': x_save_csv__mutmut_10, 
    'x_save_csv__mutmut_11': x_save_csv__mutmut_11, 
    'x_save_csv__mutmut_12': x_save_csv__mutmut_12, 
    'x_save_csv__mutmut_13': x_save_csv__mutmut_13, 
    'x_save_csv__mutmut_14': x_save_csv__mutmut_14, 
    'x_save_csv__mutmut_15': x_save_csv__mutmut_15, 
    'x_save_csv__mutmut_16': x_save_csv__mutmut_16, 
    'x_save_csv__mutmut_17': x_save_csv__mutmut_17, 
    'x_save_csv__mutmut_18': x_save_csv__mutmut_18, 
    'x_save_csv__mutmut_19': x_save_csv__mutmut_19, 
    'x_save_csv__mutmut_20': x_save_csv__mutmut_20, 
    'x_save_csv__mutmut_21': x_save_csv__mutmut_21, 
    'x_save_csv__mutmut_22': x_save_csv__mutmut_22, 
    'x_save_csv__mutmut_23': x_save_csv__mutmut_23, 
    'x_save_csv__mutmut_24': x_save_csv__mutmut_24, 
    'x_save_csv__mutmut_25': x_save_csv__mutmut_25, 
    'x_save_csv__mutmut_26': x_save_csv__mutmut_26, 
    'x_save_csv__mutmut_27': x_save_csv__mutmut_27, 
    'x_save_csv__mutmut_28': x_save_csv__mutmut_28, 
    'x_save_csv__mutmut_29': x_save_csv__mutmut_29, 
    'x_save_csv__mutmut_30': x_save_csv__mutmut_30, 
    'x_save_csv__mutmut_31': x_save_csv__mutmut_31, 
    'x_save_csv__mutmut_32': x_save_csv__mutmut_32, 
    'x_save_csv__mutmut_33': x_save_csv__mutmut_33, 
    'x_save_csv__mutmut_34': x_save_csv__mutmut_34, 
    'x_save_csv__mutmut_35': x_save_csv__mutmut_35, 
    'x_save_csv__mutmut_36': x_save_csv__mutmut_36
}
x_save_csv__mutmut_orig.__name__ = 'x_save_csv'

# Create fasta file based on database
def create_fasta(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	args = [database, metadata, submission_dir, config_dict]# type: ignore
	kwargs = {}# type: ignore
	return _mutmut_trampoline(x_create_fasta__mutmut_orig, x_create_fasta__mutmut_mutants, args, kwargs, None)

# Create fasta file based on database
def x_create_fasta__mutmut_orig(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_1(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = None
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_2(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = None
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_3(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] - "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_4(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "XXsample_nameXX"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_5(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "SAMPLE_NAME"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_6(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) or row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_7(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata or pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_8(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True or "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_9(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict or config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_10(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "XXAdd_Definition_Line_AccessionsXX" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_11(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "add_definition_line_accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_12(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "ADD_DEFINITION_LINE_ACCESSIONS" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_13(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" not in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_14(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["XXAdd_Definition_Line_AccessionsXX"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_15(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["add_definition_line_accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_16(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["ADD_DEFINITION_LINE_ACCESSIONS"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_17(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] != True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_18(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == False and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_19(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "XXbioprojectXX" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_20(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "BIOPROJECT" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_21(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" not in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_22(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(None) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_23(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["XXbioprojectXX"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_24(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["BIOPROJECT"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_25(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["XXbioprojectXX"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_26(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["BIOPROJECT"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_27(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() == "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_28(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "XXXX":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_29(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = None
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_30(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['XXbioprojectXX']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_31(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['BIOPROJECT']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_32(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = None
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_33(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = "XXXX"
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_34(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) or row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_35(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata or pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_36(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database or "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_37(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "XXGENBANKXX" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_38(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "genbank" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_39(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" not in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_40(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "XXgb-fasta_definition_line_modifiersXX" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_41(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "GB-FASTA_DEFINITION_LINE_MODIFIERS" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_42(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" not in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_43(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(None) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_44(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["XXgb-fasta_definition_line_modifiersXX"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_45(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["GB-FASTA_DEFINITION_LINE_MODIFIERS"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_46(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["XXgb-fasta_definition_line_modifiersXX"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_47(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["GB-FASTA_DEFINITION_LINE_MODIFIERS"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_48(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() == "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_49(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "XXXX":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_50(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(None)
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_51(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(None, id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_52(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =None, description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_53(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = None))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_54(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_55(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_56(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), ))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_57(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["XXfasta_sequence_origXX"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_58(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["FASTA_SEQUENCE_ORIG"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_59(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " - row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_60(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession - " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_61(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() - bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_62(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + "XX XX" + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_63(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["XXgb-fasta_definition_line_modifiersXX"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_64(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["GB-FASTA_DEFINITION_LINE_MODIFIERS"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_65(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = "XXXX"))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_66(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(None)
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_67(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(None, id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_68(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = None, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_69(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = None))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_70(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_71(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_72(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, ))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_73(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["XXfasta_sequence_origXX"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_74(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["FASTA_SEQUENCE_ORIG"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_75(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() - bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_76(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = "XXXX"))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_77(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(None, "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_78(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), None) as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_79(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open("w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_80(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), ) as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_81(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(None, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_82(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, None), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_83(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join("sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_84(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, ), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_85(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "XXsequence.fsaXX"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_86(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "SEQUENCE.FSA"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_87(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "XXw+XX") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_88(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "W+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_89(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(None, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_90(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, None, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_91(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, None)
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_92(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_93(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_94(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, )
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_95(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "XXfastaXX")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_96(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "FASTA")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_97(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(None, file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_98(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=None)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_99(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_100(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", )
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_101(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(None, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_102(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=None)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_103(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_104(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, )
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_105(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(None)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_106(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(2)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_107(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(None, file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_108(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=None)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_109(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_110(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", )
		print(e, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_111(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(None, file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_112(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=None)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_113(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(file=sys.stderr)
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_114(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, )
		sys.exit(1)

# Create fasta file based on database
def x_create_fasta__mutmut_115(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(None)

# Create fasta file based on database
def x_create_fasta__mutmut_116(database: str, metadata: pd.DataFrame, submission_dir: str, config_dict: dict[str, Any]) -> None:
	records = []
	for index, row in metadata.iterrows():
		column_name = SAMPLE_NAME_DATABASE_PREFIX[database] + "sample_name"
		if "Add_Definition_Line_Accessions" in config_dict and config_dict["Add_Definition_Line_Accessions"] == True and "bioproject" in metadata and pd.notnull(row["bioproject"]) and row["bioproject"].strip() != "":
			bioproject_accession = f" [BioProject={row['bioproject']}]"
		else:
			bioproject_accession = ""
		if "GENBANK" in database and "gb-fasta_definition_line_modifiers" in metadata and pd.notnull(row["gb-fasta_definition_line_modifiers"]) and row["gb-fasta_definition_line_modifiers"].strip() != "":
			records.append(SeqRecord(row["fasta_sequence_orig"], id =(row[column_name].strip() + bioproject_accession + " " + row["gb-fasta_definition_line_modifiers"].strip()), description = ""))
		else:
			records.append(SeqRecord(row["fasta_sequence_orig"], id = row[column_name].strip() + bioproject_accession, description = ""))
	try:
		with open(os.path.join(submission_dir, "sequence.fsa"), "w+") as f:
			SeqIO.write(records, f, "fasta")
	except PermissionError as e:
		print(f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(1)
	except Exception as e:
		print(f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}", file=sys.stderr)
		print(e, file=sys.stderr)
		sys.exit(2)

x_create_fasta__mutmut_mutants : ClassVar[MutantDict] = { # type: ignore
'x_create_fasta__mutmut_1': x_create_fasta__mutmut_1, 
    'x_create_fasta__mutmut_2': x_create_fasta__mutmut_2, 
    'x_create_fasta__mutmut_3': x_create_fasta__mutmut_3, 
    'x_create_fasta__mutmut_4': x_create_fasta__mutmut_4, 
    'x_create_fasta__mutmut_5': x_create_fasta__mutmut_5, 
    'x_create_fasta__mutmut_6': x_create_fasta__mutmut_6, 
    'x_create_fasta__mutmut_7': x_create_fasta__mutmut_7, 
    'x_create_fasta__mutmut_8': x_create_fasta__mutmut_8, 
    'x_create_fasta__mutmut_9': x_create_fasta__mutmut_9, 
    'x_create_fasta__mutmut_10': x_create_fasta__mutmut_10, 
    'x_create_fasta__mutmut_11': x_create_fasta__mutmut_11, 
    'x_create_fasta__mutmut_12': x_create_fasta__mutmut_12, 
    'x_create_fasta__mutmut_13': x_create_fasta__mutmut_13, 
    'x_create_fasta__mutmut_14': x_create_fasta__mutmut_14, 
    'x_create_fasta__mutmut_15': x_create_fasta__mutmut_15, 
    'x_create_fasta__mutmut_16': x_create_fasta__mutmut_16, 
    'x_create_fasta__mutmut_17': x_create_fasta__mutmut_17, 
    'x_create_fasta__mutmut_18': x_create_fasta__mutmut_18, 
    'x_create_fasta__mutmut_19': x_create_fasta__mutmut_19, 
    'x_create_fasta__mutmut_20': x_create_fasta__mutmut_20, 
    'x_create_fasta__mutmut_21': x_create_fasta__mutmut_21, 
    'x_create_fasta__mutmut_22': x_create_fasta__mutmut_22, 
    'x_create_fasta__mutmut_23': x_create_fasta__mutmut_23, 
    'x_create_fasta__mutmut_24': x_create_fasta__mutmut_24, 
    'x_create_fasta__mutmut_25': x_create_fasta__mutmut_25, 
    'x_create_fasta__mutmut_26': x_create_fasta__mutmut_26, 
    'x_create_fasta__mutmut_27': x_create_fasta__mutmut_27, 
    'x_create_fasta__mutmut_28': x_create_fasta__mutmut_28, 
    'x_create_fasta__mutmut_29': x_create_fasta__mutmut_29, 
    'x_create_fasta__mutmut_30': x_create_fasta__mutmut_30, 
    'x_create_fasta__mutmut_31': x_create_fasta__mutmut_31, 
    'x_create_fasta__mutmut_32': x_create_fasta__mutmut_32, 
    'x_create_fasta__mutmut_33': x_create_fasta__mutmut_33, 
    'x_create_fasta__mutmut_34': x_create_fasta__mutmut_34, 
    'x_create_fasta__mutmut_35': x_create_fasta__mutmut_35, 
    'x_create_fasta__mutmut_36': x_create_fasta__mutmut_36, 
    'x_create_fasta__mutmut_37': x_create_fasta__mutmut_37, 
    'x_create_fasta__mutmut_38': x_create_fasta__mutmut_38, 
    'x_create_fasta__mutmut_39': x_create_fasta__mutmut_39, 
    'x_create_fasta__mutmut_40': x_create_fasta__mutmut_40, 
    'x_create_fasta__mutmut_41': x_create_fasta__mutmut_41, 
    'x_create_fasta__mutmut_42': x_create_fasta__mutmut_42, 
    'x_create_fasta__mutmut_43': x_create_fasta__mutmut_43, 
    'x_create_fasta__mutmut_44': x_create_fasta__mutmut_44, 
    'x_create_fasta__mutmut_45': x_create_fasta__mutmut_45, 
    'x_create_fasta__mutmut_46': x_create_fasta__mutmut_46, 
    'x_create_fasta__mutmut_47': x_create_fasta__mutmut_47, 
    'x_create_fasta__mutmut_48': x_create_fasta__mutmut_48, 
    'x_create_fasta__mutmut_49': x_create_fasta__mutmut_49, 
    'x_create_fasta__mutmut_50': x_create_fasta__mutmut_50, 
    'x_create_fasta__mutmut_51': x_create_fasta__mutmut_51, 
    'x_create_fasta__mutmut_52': x_create_fasta__mutmut_52, 
    'x_create_fasta__mutmut_53': x_create_fasta__mutmut_53, 
    'x_create_fasta__mutmut_54': x_create_fasta__mutmut_54, 
    'x_create_fasta__mutmut_55': x_create_fasta__mutmut_55, 
    'x_create_fasta__mutmut_56': x_create_fasta__mutmut_56, 
    'x_create_fasta__mutmut_57': x_create_fasta__mutmut_57, 
    'x_create_fasta__mutmut_58': x_create_fasta__mutmut_58, 
    'x_create_fasta__mutmut_59': x_create_fasta__mutmut_59, 
    'x_create_fasta__mutmut_60': x_create_fasta__mutmut_60, 
    'x_create_fasta__mutmut_61': x_create_fasta__mutmut_61, 
    'x_create_fasta__mutmut_62': x_create_fasta__mutmut_62, 
    'x_create_fasta__mutmut_63': x_create_fasta__mutmut_63, 
    'x_create_fasta__mutmut_64': x_create_fasta__mutmut_64, 
    'x_create_fasta__mutmut_65': x_create_fasta__mutmut_65, 
    'x_create_fasta__mutmut_66': x_create_fasta__mutmut_66, 
    'x_create_fasta__mutmut_67': x_create_fasta__mutmut_67, 
    'x_create_fasta__mutmut_68': x_create_fasta__mutmut_68, 
    'x_create_fasta__mutmut_69': x_create_fasta__mutmut_69, 
    'x_create_fasta__mutmut_70': x_create_fasta__mutmut_70, 
    'x_create_fasta__mutmut_71': x_create_fasta__mutmut_71, 
    'x_create_fasta__mutmut_72': x_create_fasta__mutmut_72, 
    'x_create_fasta__mutmut_73': x_create_fasta__mutmut_73, 
    'x_create_fasta__mutmut_74': x_create_fasta__mutmut_74, 
    'x_create_fasta__mutmut_75': x_create_fasta__mutmut_75, 
    'x_create_fasta__mutmut_76': x_create_fasta__mutmut_76, 
    'x_create_fasta__mutmut_77': x_create_fasta__mutmut_77, 
    'x_create_fasta__mutmut_78': x_create_fasta__mutmut_78, 
    'x_create_fasta__mutmut_79': x_create_fasta__mutmut_79, 
    'x_create_fasta__mutmut_80': x_create_fasta__mutmut_80, 
    'x_create_fasta__mutmut_81': x_create_fasta__mutmut_81, 
    'x_create_fasta__mutmut_82': x_create_fasta__mutmut_82, 
    'x_create_fasta__mutmut_83': x_create_fasta__mutmut_83, 
    'x_create_fasta__mutmut_84': x_create_fasta__mutmut_84, 
    'x_create_fasta__mutmut_85': x_create_fasta__mutmut_85, 
    'x_create_fasta__mutmut_86': x_create_fasta__mutmut_86, 
    'x_create_fasta__mutmut_87': x_create_fasta__mutmut_87, 
    'x_create_fasta__mutmut_88': x_create_fasta__mutmut_88, 
    'x_create_fasta__mutmut_89': x_create_fasta__mutmut_89, 
    'x_create_fasta__mutmut_90': x_create_fasta__mutmut_90, 
    'x_create_fasta__mutmut_91': x_create_fasta__mutmut_91, 
    'x_create_fasta__mutmut_92': x_create_fasta__mutmut_92, 
    'x_create_fasta__mutmut_93': x_create_fasta__mutmut_93, 
    'x_create_fasta__mutmut_94': x_create_fasta__mutmut_94, 
    'x_create_fasta__mutmut_95': x_create_fasta__mutmut_95, 
    'x_create_fasta__mutmut_96': x_create_fasta__mutmut_96, 
    'x_create_fasta__mutmut_97': x_create_fasta__mutmut_97, 
    'x_create_fasta__mutmut_98': x_create_fasta__mutmut_98, 
    'x_create_fasta__mutmut_99': x_create_fasta__mutmut_99, 
    'x_create_fasta__mutmut_100': x_create_fasta__mutmut_100, 
    'x_create_fasta__mutmut_101': x_create_fasta__mutmut_101, 
    'x_create_fasta__mutmut_102': x_create_fasta__mutmut_102, 
    'x_create_fasta__mutmut_103': x_create_fasta__mutmut_103, 
    'x_create_fasta__mutmut_104': x_create_fasta__mutmut_104, 
    'x_create_fasta__mutmut_105': x_create_fasta__mutmut_105, 
    'x_create_fasta__mutmut_106': x_create_fasta__mutmut_106, 
    'x_create_fasta__mutmut_107': x_create_fasta__mutmut_107, 
    'x_create_fasta__mutmut_108': x_create_fasta__mutmut_108, 
    'x_create_fasta__mutmut_109': x_create_fasta__mutmut_109, 
    'x_create_fasta__mutmut_110': x_create_fasta__mutmut_110, 
    'x_create_fasta__mutmut_111': x_create_fasta__mutmut_111, 
    'x_create_fasta__mutmut_112': x_create_fasta__mutmut_112, 
    'x_create_fasta__mutmut_113': x_create_fasta__mutmut_113, 
    'x_create_fasta__mutmut_114': x_create_fasta__mutmut_114, 
    'x_create_fasta__mutmut_115': x_create_fasta__mutmut_115, 
    'x_create_fasta__mutmut_116': x_create_fasta__mutmut_116
}
x_create_fasta__mutmut_orig.__name__ = 'x_create_fasta'
