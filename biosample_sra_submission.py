#!/usr/bin/env python3

import ftplib
import os
import sys
import pandas as pd
import yaml

pd.set_option("display.max_columns", None)
pd.set_option("max_colwidth", None)
pd.set_option("display.max_rows", None)
config_dict = dict()

#Initialize config file
def initialize_global_variables(config):
    if os.path.isfile(config) == False:
        print("Error: Cannot find submission config at: " + config, file=sys.stderr)
        sys.exit(1)
    else:
        with open(config, "r") as f:
            global config_dict
            config_dict = yaml.safe_load(f)
        if isinstance(config_dict, dict) == False:
            print("Config Error: Config file structure is incorrect.", file=sys.stderr)
            sys.exit(1)

def submit_ftp(unique_name, ncbi_sub_type, config, test, overwrite):
    initialize_global_variables(config)
    submit_file_path: str = os.path.join(config_dict["general"]["submission_directory"], unique_name, "submit.ready")
    if not os.path.isfile(submit_file_path):
        open(submit_file_path, 'w+').close()
    try:
        #Login to ftp
        ftp = ftplib.FTP(config_dict["ncbi"]["hostname"])
        ftp.login(user=config_dict["ncbi"]["username"], passwd = config_dict["ncbi"]["password"])
        if config_dict["ncbi"]["ncbi_ftp_path_to_submission_folders"] != "":
            ftp.cwd(config_dict["ncbi"]["ncbi_ftp_path_to_submission_folders"])
        if test == False:
            ftp.cwd("Production")
        else:
            ftp.cwd("Test")
        dir = unique_name + "_" + ncbi_sub_type
        if dir not in ftp.nlst():
            ftp.mkd(dir)
        ftp.cwd(dir)
        #Check if report.xml exists
        if "report.xml" in ftp.nlst() and overwrite == False:
            print("Submission report exists pulling down.")
            report_file = open(os.path.join(config_dict["general"]["submission_directory"], unique_name, "biosample_sra", unique_name + "_" + ncbi_sub_type + "_report.xml"), 'wb')
            ftp.retrbinary('RETR report.xml', report_file.write, 262144)
            report_file.close()
        else:
            print("Submitting to SRA/BioSample.")
            res = ftp.storlines("STOR " + "submission.xml", open(os.path.join(config_dict["general"]["submission_directory"], unique_name, "biosample_sra", unique_name + "_" + ncbi_sub_type + "_submission.xml"), 'rb'))
            if not res.startswith('226 Transfer complete'):
                print('Submission.xml upload failed.')
            if config_dict["ncbi"]["SRA_file_location"] == "local" and "sra" in ncbi_sub_type:
                with open(os.path.join(config_dict["general"]["submission_directory"], unique_name, "biosample_sra", "sra_file_path.txt"), "r") as file:
                    for line in file:
                        res = ftp.storbinary("STOR " + os.path.basename(line.strip()), open(line.strip(), 'rb'))
                        if not res.startswith('226 Transfer complete'):
                            print('SRA file upload failed. Try again.')
                            sys.exit(1)
            res = ftp.storlines("STOR " + "submit.ready", open(submit_file_path, 'rb'))
            if not res.startswith('226 Transfer complete'):
                print('submit.ready upload failed.')
    except ftplib.all_errors as e:
        print('FTP error:', e)
        sys.exit()
