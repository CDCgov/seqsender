#!/usr/bin/env python3

import ftplib
import argparse
from argparse import RawTextHelpFormatter
import submission_preparation
import gisaid_submission
import genbank_submission
import biosample_sra_submission
import sys
from datetime import datetime
import requests
import json
import os
import pandas as pd
import yaml
from Bio import SeqIO
import xml.etree.ElementTree as ET

config_dict = dict()
_version = "0.1 (Beta)"

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

#Process Biosample Report
def biosample_sra_process_report(unique_name, ncbi_sub_type):
    submission_status = ""
    submission_id = "pending"
    df = pd.read_csv(os.path.join(config_dict["general"]["submission_directory"], unique_name, "accessions.csv"), header = 0, dtype = str, sep = ",")
    if "BioSample_accession" not in df and "biosample" in ncbi_sub_type:
        df["BioSample_accession"] = ""
    if "SRA_accession" not in df and "sra" in ncbi_sub_type:
        df["SRA_accession"] = ""
    with open(os.path.join(config_dict["general"]["submission_directory"], unique_name, "biosample_sra", unique_name + "_" + ncbi_sub_type + "_report.xml"), 'r') as file:
        line = file.readline()
        while line:
            if "<SubmissionStatus status=" in line:
                submission_status = line.split("<SubmissionStatus status=\"")[-1].split("\" submission_id=")[0].split("\">")[0]
                if submission_id == "pending" and "submission_id=" in line:
                    submission_id = line.split("\" submission_id=\"")[-1].split("\" last_update=")[0]
            if "<Object target_db=\"BioSample\"" in line and "accession=" in line:
                df.loc[df.BioSample_sequence == (line.split("spuid=\"")[-1].split("\" spuid_namespace=")[0]), "BioSample_accession"] = line.split("accession=\"")[-1].split("\" spuid=")[0]
            if "<Object target_db=\"SRA\"" in line and "accession=" in line:
                df.loc[df.SRA_sequence == (line.split("spuid=\"")[-1].split("\" spuid_namespace=")[0]), "SRA_accession"] = line.split("accession=\"")[-1].split("\" spuid=")[0]
            line = file.readline()
    df.to_csv(os.path.join(config_dict["general"]["submission_directory"], unique_name, "accessions.csv"), header = True, index = False, sep = ",")
    return submission_id, submission_status

#Update log csv
def update_csv(
    unique_name,
    config,
    upload_log_path,
    type,
    Genbank_submission_id=None,
    Genbank_submission_date=None,
    Genbank_status=None,
    SRA_submission_id=None,
    SRA_status=None,
    SRA_submission_date=None,
    Biosample_submission_id=None,
    Biosample_status=None,
    Biosample_submission_date=None,
    GISAID_submission_date=None,
    GISAID_submitted_total=None,
    GISAID_failed_total=None
):
    curr_time = datetime.now()
    if os.path.isfile(upload_log_path):
        df = _read_upload_log(upload_log_path)
    else:
        df = pd.DataFrame(
            columns=[
                "name",
                "update_date",
                "SRA_submission_id",
                "SRA_submission_date",
                "SRA_status",
                "BioSample_submission_id",
                "BioSample_submission_date",
                "BioSample_status",
                "Genbank_submission_id",
                "Genbank_submission_date",
                "Genbank_status",
                "GISAID_submission_date",
                "GISAID_submitted_total",
                "GISAID_failed_total",
                "directory",
                "config",
                "type"
            ]
        )
    #Check if row exists in log to update instead of write new
    if df['name'].str.contains(unique_name).any():
        df_partial = df.loc[df['name'] == unique_name]
        df.loc[df_partial.index.values, 'update_date'] = curr_time.strftime("%-m/%-d/%Y")
        df.loc[df_partial.index.values, 'directory'] = os.path.join(config_dict["general"]["submission_directory"], unique_name)
        df.loc[df_partial.index.values, 'config'] = config
        df.loc[df_partial.index.values, 'type'] = type
        if Genbank_submission_id is not None:
            df.loc[df_partial.index.values, 'Genbank_submission_id'] = Genbank_submission_id
        if Genbank_submission_date is not None:
            df.loc[df_partial.index.values, 'Genbank_submission_date'] = Genbank_submission_date
        if Genbank_status is not None:
            df.loc[df_partial.index.values, 'Genbank_status'] = Genbank_status
        if SRA_submission_id is not None:
            df.loc[df_partial.index.values, 'SRA_submission_id'] = SRA_submission_id
        if SRA_submission_date is not None:
            df.loc[df_partial.index.values, 'SRA_submission_date'] = SRA_submission_date
        if SRA_status is not None:
            df.loc[df_partial.index.values, 'SRA_status'] = SRA_status
        if Biosample_submission_id is not None:
            df.loc[df_partial.index.values, 'BioSample_submission_id'] = Biosample_submission_id
        if Biosample_submission_date is not None:
            df.loc[df_partial.index.values, 'BioSample_submission_date'] = Biosample_submission_date
        if Biosample_status is not None:
            df.loc[df_partial.index.values, 'BioSample_status'] = Biosample_status
        if GISAID_submission_date is not None:
            df.loc[df_partial.index.values, 'GISAID_submission_date'] = GISAID_submission_date
        if GISAID_submitted_total is not None:
            df.loc[df_partial.index.values, 'GISAID_submitted_total'] = GISAID_submitted_total
        if GISAID_failed_total is not None:
            df.loc[df_partial.index.values, 'GISAID_failed_total'] = GISAID_failed_total
    else:
        new_entry = {
            "name":unique_name,
            "update_date":curr_time.strftime("%-m/%-d/%Y"),
            "Genbank_submission_id":Genbank_submission_id,
            "Genbank_submission_date":Genbank_submission_date,
            "Genbank_status":Genbank_status,
            "directory":os.path.join(config_dict["general"]["submission_directory"], unique_name),
            "config":config,
            "type":type,
            "SRA_submission_id":SRA_submission_id,
            "SRA_submission_date":SRA_submission_date,
            "SRA_status":SRA_status,
            "BioSample_submission_id":Biosample_submission_id,
            "BioSample_submission_date":Biosample_submission_date,
            "BioSample_status":Biosample_status,
            "GISAID_submission_date":GISAID_submission_date,
            "GISAID_submitted_total":GISAID_submitted_total,
            "GISAID_failed_total":GISAID_failed_total
        }
        df = pd.concat([df, pd.Series(new_entry).to_frame().T], ignore_index=True)
    df.to_csv(upload_log_path, header = True, index = False, sep = ",")


def _read_upload_log(upload_log_path: str) -> pd.DataFrame:
    return pd.read_csv(upload_log_path, header=0, dtype=str, sep=",")


#Update log status
#Pulls all entries that do not say processed and updates status
def update_log(upload_log_path: str):
    if os.path.isfile(upload_log_path):
        main_df = _read_upload_log(upload_log_path)
    else:
        error_msg = f"Error: Either a submission has not been made or upload_log.csv has been moved from {upload_log_path}"
        print(error_msg, file=sys.stderr)
        raise RuntimeError(error_msg)
    #Biosample/SRA
    df = main_df.loc[(main_df['BioSample_status'] != None) & (main_df['BioSample_status'] != "processed-ok") & (main_df['SRA_status'] != None) & (main_df['SRA_status'] != "processed-ok") & (main_df['BioSample_status'] != "") & (main_df['SRA_status'] != "")]
    if len(df.index) != 0:
        for index, row in df.iterrows():
            report_generated = False
            try:
                initialize_global_variables(row["config"])
                if config_dict["general"]["submit_BioSample"] == True and config_dict["general"]["submit_SRA"] == True and config_dict["general"]["joint_SRA_BioSample_submission"] == True:
                    print("\nUpdating: " + row["name"] + " BioSample/SRA")
                    #Login to ftp
                    ftp = ftplib.FTP(config_dict["ncbi"]["hostname"])
                    ftp.login(user=config_dict["ncbi"]["username"], passwd = config_dict["ncbi"]["password"])
                    if config_dict["ncbi"]["ncbi_ftp_path_to_submission_folders"] != "":
                        ftp.cwd(config_dict["ncbi"]["ncbi_ftp_path_to_submission_folders"])
                    ftp.cwd(row["type"])
                    submission_folder = row['name'] + "_biosample_sra"
                    #Check if submission folder already exists
                    if submission_folder not in ftp.nlst():
                        print("submission doesn't exist")
                        continue
                    ftp.cwd(submission_folder)
                    #Check if report.xml exists
                    if "report.xml" in ftp.nlst():
                        print("Report exists pulling down")
                        report_file = open(os.path.join(config_dict["general"]["submission_directory"], row["name"], "biosample_sra", row["name"] + "_biosample_sra_report.xml"), 'wb')
                        ftp.retrbinary('RETR report.xml', report_file.write, 262144)
                        report_file.close()
                        report_generated = True
            except ftplib.all_errors as e:
                print('FTP error:', e)
            if report_generated == True:
                submission_id, submission_status = biosample_sra_process_report(row["name"], "biosample_sra")
                update_csv(
                    unique_name=row['name'],
                    config=row["config"],
                    upload_log_path=upload_log_path,
                    type=row["type"],
                    Biosample_submission_id=submission_id,
                    Biosample_status=submission_status,
                    SRA_submission_id=submission_id,
                    SRA_status=submission_status
                )
                if submission_status == "processed-ok" and config_dict["general"]["submit_Genbank"] == True and (row["Genbank_status"] == None or row["Genbank_status"] == "" or pd.isnull(row["Genbank_status"])):
                    print("Submitting to Genbank: " + row["name"])
                    submit_genbank(row["name"], row["config"], upload_log_path, row["type"], False)
    #Check BioSample
    main_df = _read_upload_log(upload_log_path)
    df = main_df.loc[(main_df['BioSample_status'] != None) & (main_df['BioSample_status'] != "processed-ok") & (main_df['BioSample_status'] != "")]
    if len(df.index) != 0:
        for index, row in df.iterrows():
            report_generated = False
            try:
                initialize_global_variables(row["config"])
                if config_dict["general"]["submit_BioSample"] == True:
                    print("\nUpdating: " + row["name"] + " BioSample")
                    #Login to ftp
                    ftp = ftplib.FTP(config_dict["ncbi"]["hostname"])
                    ftp.login(user=config_dict["ncbi"]["username"], passwd = config_dict["ncbi"]["password"])
                    if config_dict["ncbi"]["ncbi_ftp_path_to_submission_folders"] != "":
                        ftp.cwd(config_dict["ncbi"]["ncbi_ftp_path_to_submission_folders"])
                    ftp.cwd(row["type"])
                    submission_folder = row['name'] + "_biosample"
                    #Check if submission folder already exists
                    if submission_folder not in ftp.nlst():
                        print("submission doesn't exist")
                        continue
                    ftp.cwd(submission_folder)
                    #Check if report.xml exists
                    if "report.xml" in ftp.nlst():
                        print("Report exists pulling down")
                        report_file = open(os.path.join(config_dict["general"]["submission_directory"], row["name"], "biosample_sra", row["name"] + "_biosample_report.xml"), 'wb')
                        ftp.retrbinary('RETR report.xml', report_file.write, 262144)
                        report_file.close()
                        report_generated = True
            except ftplib.all_errors as e:
                print('FTP error:', e)
            if report_generated == True:
                submission_id, submission_status = biosample_sra_process_report(row["name"], "biosample")
                update_csv(
                    unique_name=row['name'],
                    config=row["config"],
                    upload_log_path=upload_log_path,
                    type=row["type"],
                    Biosample_submission_id=submission_id,
                    Biosample_status=submission_status
                )
                print("Status: " + submission_status)
                if submission_status == "processed-ok" and config_dict["general"]["submit_Genbank"] == True and (row["Genbank_status"] == None or row["Genbank_status"] == "" or pd.isnull(row["Genbank_status"])):
                    print("Submitting to Genbank: " + row["name"])
                    submit_genbank(row["name"], row["config"], upload_log_path, row["type"], True)
    #Check SRA
    main_df = _read_upload_log(upload_log_path)
    df = main_df.loc[(main_df['SRA_status'] != None) & (main_df['SRA_status'] != "processed-ok") & (main_df['SRA_status'] != "")]
    if len(df.index) != 0:
        for index, row in df.iterrows():
            report_generated = False
            try:
                initialize_global_variables(row["config"])
                if config_dict["general"]["submit_SRA"] == True:
                    print("\nUpdating: " + row["name"] + " SRA")
                    #Login to ftp
                    ftp = ftplib.FTP(config_dict["ncbi"]["hostname"])
                    ftp.login(user=config_dict["ncbi"]["username"], passwd = config_dict["ncbi"]["password"])
                    if config_dict["ncbi"]["ncbi_ftp_path_to_submission_folders"] != "":
                        ftp.cwd(config_dict["ncbi"]["ncbi_ftp_path_to_submission_folders"])
                    ftp.cwd(row["type"])
                    submission_folder = row['name'] + "_sra"
                    #Check if submission folder already exists
                    if submission_folder not in ftp.nlst():
                        print("submission doesn't exist")
                        continue
                    ftp.cwd(submission_folder)
                    #Check if report.xml exists
                    if "report.xml" in ftp.nlst():
                        print("Report exists pulling down")
                        report_file = open(os.path.join(config_dict["general"]["submission_directory"], row["name"], "biosample_sra", row["name"] + "_sra_report.xml"), 'wb')
                        ftp.retrbinary('RETR report.xml', report_file.write, 262144)
                        report_file.close()
                        report_generated = True
            except ftplib.all_errors as e:
                print('FTP error:', e)
            if report_generated == True:
                submission_id, submission_status = biosample_sra_process_report(row["name"], "sra")
                update_csv(
                    unique_name=row['name'],
                    config=row["config"],
                    upload_log_path=upload_log_path,
                    type=row["type"],
                    SRA_submission_id=submission_id,
                    SRA_status=submission_status
                )
                print("Status: " + submission_status)
                if submission_status == "processed-ok" and config_dict["general"]["submit_Genbank"] == True and (row["Genbank_status"] == None or row["Genbank_status"] == "" or pd.isnull(row["Genbank_status"])):
                    print("Submitting to Genbank: " + row["name"])
                    submit_genbank(row["name"], row["config"], upload_log_path, row["type"], False)
    #Check Genbank
    main_df = _read_upload_log(upload_log_path)
    df = main_df.loc[(main_df['Genbank_status'] != None) & (main_df['Genbank_status'] != "processed-ok") & (main_df['Genbank_status'] != "") & (main_df['Genbank_status'].isnull() == False)]
    if len(df.index) != 0:
        for index, row in df.iterrows():
            report_generated = False
            try:
                initialize_global_variables(row["config"])
                if config_dict["general"]["submit_Genbank"] == True:
                    print("\nUpdating: " + row["name"] + " Genbank")
                    #Login to ftp
                    ftp = ftplib.FTP(config_dict["ncbi"]["hostname"])
                    ftp.login(user=config_dict["ncbi"]["username"], passwd = config_dict["ncbi"]["password"])
                    if config_dict["ncbi"]["ncbi_ftp_path_to_submission_folders"] != "":
                        ftp.cwd(config_dict["ncbi"]["ncbi_ftp_path_to_submission_folders"])
                    ftp.cwd(row["type"])
                    submission_folder = row['name'] + "_genbank"
                    #Check if submission folder already exists
                    if submission_folder not in ftp.nlst():
                        print("submission doesn't exist")
                        continue
                    ftp.cwd(submission_folder)
                    #Check if report.xml exists
                    if "report.xml" in ftp.nlst():
                        print("Report exists pulling down")
                        report_file = open(os.path.join(config_dict["general"]["submission_directory"], row["name"], "genbank", row["name"] + "_report.xml"), 'wb')
                        ftp.retrbinary('RETR report.xml', report_file.write, 262144)
                        report_file.close()
                        report_generated = True
            except ftplib.all_errors as e:
                print('FTP error:', e)
            if report_generated == True:
                submission_id, submission_status = genbank_process_report(row["name"])
                update_csv(
                    unique_name=row['name'],
                    config=row["config"],
                    upload_log_path=upload_log_path,
                    type=row["type"],
                    Genbank_submission_id=submission_id,
                    Genbank_status=submission_status
                )
                print("Status: " + submission_status)
                if submission_status == "processed-ok" and config_dict["general"]["submit_GISAID"].lower == "true" and row["type"] != "Test":
                    print("\nSubmitting to GISAID: " + row["name"])
                    submit_gisaid(unique_name=row["name"], config=row["config"], upload_log_path=upload_log_path, test=row["type"])
    #Check GISAID
    main_df = _read_upload_log(upload_log_path)
    df = main_df.loc[(main_df['GISAID_submitted_total'] != None) & (main_df['GISAID_submitted_total'] != "") & (main_df['GISAID_submitted_total'].isnull() == False) & (main_df['type'] != "Test") & (main_df['GISAID_failed_total'] != "0")]
    if len(df.index) != 0:
        for index, row in df.iterrows():
            initialize_global_variables(row["config"])
            if config_dict["general"]["submit_GISAID"] == True:
                print("\nSubmitting to GISAID: " + row["name"])
                submit_gisaid(unique_name=row["name"], config=row["config"], upload_log_path=upload_log_path, test=row["type"])

#Read output log from gisaid submission script
def read_log(unique_name, file):
    if os.path.exists(file):
        with open(file) as f:
            data = json.load(f)
        number_submitted = 0
        number_failed = 0
        already_submitted = []
        for i in data:
            #Sequence successfully uploaded
            if i["code"] == "upload_count":
                number_submitted = int(i["msg"].strip().split("uploaded: ")[1])
            #Sequence failed upload
            elif i["code"] == "failed_count":
                number_failed = int(i["msg"].strip().split("failed: ")[1])
            #Correct number of successfully uploaded for if a sequence fails for already existing
            elif (i["code"] == "validation_error") and ("\"covv_virus_name\": \"already exists\"" in i["msg"]):
                already_submitted.append(i["msg"].split("; validation_error;")[0])
            elif i["code"] == "epi_isl_id":
                continue
        clean_failed_log(unique_name, number_failed, already_submitted)
        number_failed = number_failed - len(already_submitted)
        number_submitted = number_submitted + len(already_submitted)
        return str(number_submitted), str(number_failed)
    else:
        return "error", "error"

#Cleans failed meta log if some of the submissions are just already submitted or
#Removes file if it is empty
def clean_failed_log(unique_name, number_failed, already_submitted):
    if number_failed == 0 or number_failed == len(already_submitted):
        if os.path.exists(os.path.join(config_dict["general"]["submission_directory"], unique_name, "gisaid", unique_name + "_failed_meta.csv")):
            print("No failed sequences.\nCleaning up files.")
            os.remove(os.path.join(config_dict["general"]["submission_directory"], unique_name, "gisaid", unique_name + "_failed_meta.csv"))
    else:
        df = pd.read_csv(os.path.join(config_dict["general"]["submission_directory"], unique_name, "gisaid", unique_name + "_failed_meta.csv"), header = 0, dtype = str)
        clean_df = df[~df.covv_virus_name.isin(already_submitted)]
        clean_df.to_csv(os.path.join(config_dict["general"]["submission_directory"], unique_name, "gisaid", unique_name + "_failed_meta.csv"), header = True, index = False)
        print("Error: Sequences failed please check: " + os.path.join(config_dict["general"]["submission_directory"], unique_name, "gisaid", unique_name + ".log"))

#Pull down report files
def pull_report_files(unique_name, files):
    api_url = config_dict["ncbi"]["api_url"]
    for item in files.keys():
        r = requests.get(api_url.replace("FILE_ID", files[item]), allow_redirects=True)
        open(os.path.join(config_dict["general"]["submission_directory"], unique_name, "genbank", unique_name + "_" + item), 'wb').write(r.content)

def submit_genbank(unique_name, config, upload_log_path, test, overwrite):
    initialize_global_variables(config)
    prepare_genbank(unique_name)
    if test == "Production":
        test_type = False
    else:
        test_type = True
    genbank_submission.submit_ftp(unique_name=unique_name, config=config, test=test_type, overwrite=overwrite)
    curr_time = datetime.now()
    update_csv(
        unique_name=unique_name,
        config=config,
        upload_log_path=upload_log_path,
        type=test,
        Genbank_submission_id="submitted",
        Genbank_submission_date=curr_time.strftime("%-m/%-d/%Y"),
        Genbank_status="submitted"
    )

def submit_gisaid(unique_name, config, upload_log_path, test):
    initialize_global_variables(config)
    if config_dict["gisaid"]["Update_sequences_on_Genbank_auto_removal"] == True and config_dict["ncbi"]["Genbank_auto_remove_sequences_that_fail_qc"] == True:
        prepare_gisaid(unique_name)
    if test == "Production":
        test_type = False
    else:
        test_type = True
    gisaid_submission.run_uploader(unique_name=unique_name, config=config, test=test_type)
    submitted, failed = read_log(unique_name, os.path.join(config_dict["general"]["submission_directory"], unique_name, "gisaid", unique_name + ".log"))
    curr_time = datetime.now()
    update_csv(
        unique_name=unique_name,
        config=config,
        upload_log_path=upload_log_path,
        type=test,
        GISAID_submission_date=curr_time.strftime("%-m/%-d/%Y"),
        GISAID_submitted_total=submitted,
        GISAID_failed_total=failed
    )

#Read xml report and check status of report
def genbank_process_report(unique_name):
    tree = ET.parse(os.path.join(config_dict["general"]["submission_directory"], unique_name, "genbank", unique_name + "_report.xml"))
    root = tree.getroot()
    if root.get('submission_id') == None:
        return "error", "error"
    files = dict()
    submission_id = ""
    status = ""
    for elem in tree.iter():
        if "submission_id" in elem.attrib.keys() and submission_id == "":
            submission_id = elem.attrib["submission_id"]
        if "status" in elem.attrib.keys() and status == "":
            if elem.attrib["status"] == "submitted":
                status = "submitted"
            elif elem.attrib["status"] == "queued":
                status = "queued"
            elif elem.attrib["status"] == "processing":
                status = "processing"
            elif elem.attrib["status"] == "processed-ok":
                status = "processed-ok"
            else:
                print("Possible Error check " + os.path.join(config_dict["general"]["submission_directory"], unique_name, "genbank", unique_name + "_report.xml"))
                status = "error"
        if "file_id" in elem.attrib.keys():
            files[elem.attrib['file_path']] = elem.attrib['file_id']
    if len(files) != 0:
        pull_report_files(unique_name, files)
    return submission_id, status

#Add biosample/SRA data to genbank submissions
def prepare_genbank(unique_name):
    accessions = pd.read_csv(os.path.join(config_dict["general"]["submission_directory"], unique_name, "accessions.csv"), header = 0, dtype=str, sep=',')
    df = pd.read_csv(os.path.join(config_dict["general"]["submission_directory"], unique_name, "genbank", unique_name + "_source.src"), header = 0, dtype=str, sep='\t')
    if config_dict["ncbi"]['BioProject'] != "":
        df['BioProject'] = config_dict["ncbi"]['BioProject']
    df = df.merge(accessions, how='left', left_on='sequence_ID', right_on='Genbank_sequence')
    potential_columns = ["SRA_sequence", "BioSample_sequence", "Genbank_sequence", "Genbank_accession", "GISAID_sequence"]
    drop_columns = []
    for col in potential_columns:
        if col in df:
            drop_columns.append(col)
    df = df.drop(columns=drop_columns)
    if "BioSample_accession" not in df and "SRA_accession" not in df:
        return
    df = df.rename(columns={"BioSample_accession": "BioSample", "SRA_accession": "SRA"})
    if "BioSample" in df:
        df["BioSample"] = df["BioSample"].fillna("")
    if "SRA" in df:
        df["SRA"] = df["SRA"].fillna("")
    col_names = df.columns.values.tolist()
    col_names.remove("sequence_ID")
    col_names.insert(0, "sequence_ID")
    df = df[col_names]
    df.to_csv(os.path.join(config_dict["general"]["submission_directory"], unique_name, "genbank", unique_name + "_source.src"), header = True, index = False, sep = '\t')

#If removing GISAID sequences based on Genbank Auto-remove
def prepare_gisaid(unique_name):
    if config_dict["gisaid"]["Update_sequences_on_Genbank_auto_removal"].lower() != "true":
        return
    accessions = pd.read_csv(os.path.join(config_dict["general"]["submission_directory"], unique_name, "accessions.csv"), header = 0, dtype=str, sep=',')
    if "Genbank_accession" not in accessions:
        return
    df = pd.read_csv(os.path.join(config_dict["general"]["submission_directory"], unique_name, "gisaid", unique_name + "_gisaid.csv"), header = 0, dtype=str, sep=',')
    df = df.merge(accessions, how='left', left_on='GISAID_sequence', right_on='covv_virus_name')
    potential_columns = ["SRA_sequence", "BioSample_sequence", "Genbank_sequence", "BioSample_accession", "SRA_accession", "GISAID_sequence", "Genbank_accession"]
    drop_columns = []
    for col in potential_columns:
        if col in df:
            drop_columns.append(col)
    df = df.dropna(subset=["Genbank_accession"])
    df = df[df.Genbank_accession != ""]
    df = df.drop(columns=[drop_columns])
    shutil.copy2(os.path.join(config_dict["general"]["submission_directory"], unique_name, "gisaid", unique_name + "_gisaid.csv"), os.path.join(config_dict["general"]["submission_directory"], unique_name, "gisaid", "old_" + unique_name + "_gisaid.csv"))
    df.to_csv(os.path.join(config_dict["general"]["submission_directory"],unique_name,"gisaid",unique_name + "_gisaid.csv"), na_rep="Unknown", index = False, header = True, quoting=csv.QUOTE_ALL)
    keep_records = []
    with open(os.path.join(config_dict["general"]["submission_directory"], unique_name, "gisaid", unique_name + "_gisaid.fsa"), "r") as fsa:
        records = SeqIO.parse(fsa, "fasta")
        for record in records:
            if record.id in df["covv_virus_name"]:
                keep_records.append(record)
    shutil.copy2(os.path.join(config_dict["general"]["submission_directory"],unique_name,"gisaid",unique_name + "_gisaid.fsa"), os.path.join(config_dict["general"]["submission_directory"],unique_name,"gisaid", "old_" + unique_name + "_gisaid.fsa"))
    with open(os.path.join(config_dict["general"]["submission_directory"],unique_name,"gisaid",unique_name + "_gisaid.fsa"), "w+") as fasta_file:
        SeqIO.write(keep_records, fasta_file, "fasta")

#For submitting when SRA/Biosample have to be split due to errors
def submit_biosample_sra(unique_name, config, upload_log_path, test, ncbi_sub_type, overwrite):
    initialize_global_variables(config)
    if test == "Production":
        test_type = False
    else:
        test_type = True
    biosample_sra_submission.submit_ftp(unique_name=unique_name, ncbi_sub_type=ncbi_sub_type, config=config, test=test_type, overwrite=overwrite)
    curr_time = datetime.now()
    if ncbi_sub_type == "biosample_sra":
        update_csv(
            unique_name=unique_name,
            config=config,
            upload_log_path=upload_log_path,
            type=test,
            Biosample_submission_id="submitted",
            Biosample_status="submitted",
            Biosample_submission_date=curr_time.strftime("%-m/%-d/%Y"),
            SRA_submission_id="submitted",
            SRA_status="submitted",
            SRA_submission_date=curr_time.strftime("%-m/%-d/%Y")
        )
    elif ncbi_sub_type == "biosample":
        update_csv(
            unique_name=unique_name,
            config=config,
            upload_log_path=upload_log_path,
            type=test,
            Biosample_submission_id="submitted",
            Biosample_status="submitted",
            Biosample_submission_date=curr_time.strftime("%-m/%-d/%Y")
        )
    elif ncbi_sub_type == "sra":
        update_csv(
            unique_name=unique_name,
            config=config,
            upload_log_path=upload_log_path,
            type=test,
            SRA_submission_id="submitted",
            SRA_status="submitted",
            SRA_submission_date=curr_time.strftime("%-m/%-d/%Y")
        )

#Start submission into automated pipeline
def start_submission(unique_name, config, upload_log_path, test, overwrite):
    initialize_global_variables(config)
    if config_dict["general"]["submit_BioSample"] == True and config_dict["general"]["submit_SRA"] == True and config_dict["general"]["joint_SRA_BioSample_submission"] == True:
        submit_biosample_sra(unique_name, config, upload_log_path, test, "biosample_sra", overwrite)
    elif config_dict["general"]["submit_BioSample"] == True and config_dict["general"]["submit_SRA"] == True and config_dict["general"]["joint_SRA_BioSample_submission"] == False:
        submit_biosample_sra(unique_name, config, upload_log_path, test, "biosample", overwrite)
        submit_biosample_sra(unique_name, config, upload_log_path, test, "sra", overwrite)
    elif config_dict["general"]["submit_BioSample"] == True:
        submit_biosample_sra(unique_name, config, upload_log_path, test, "biosample", overwrite)
    elif config_dict["general"]["submit_SRA"] == True:
        submit_biosample_sra(unique_name, config, upload_log_path, test, "sra", overwrite)
    elif config_dict["general"]["submit_Genbank"] == True:
        submit_genbank(unique_name=unique_name, config=config, upload_log_path=upload_log_path, test=test, overwrite=overwrite)
    elif config_dict["general"]["submit_GISAID"] == True:
        submit_gisaid(unique_name=unique_name, config=config, upload_log_path=upload_log_path, test=test)

def test_bioproject(config):
    initialize_global_variables(config)
    submit_file_path: str = os.path.join(config_dict["general"]["submission_directory"], "submit.ready")
    if not os.path.isfile(submit_file_path):
        open(submit_file_path, 'w+').close()
    try:
        #Login to ftp
        ftp = ftplib.FTP(config_dict["ncbi"]["hostname"])
        ftp.login(user=config_dict["ncbi"]["username"], passwd = config_dict["ncbi"]["password"])
        if config_dict["ncbi"]["ncbi_ftp_path_to_submission_folders"] != "":
            ftp.cwd(config_dict["ncbi"]["ncbi_ftp_path_to_submission_folders"])
        ftp.cwd("Test")
        dir = "Test_BioProject"
        if dir not in ftp.nlst():
            ftp.mkd(dir)
        ftp.cwd(dir)
        #Check if report.xml exists
        if "report.xml" in ftp.nlst():
            print("Submission report exists pulling down.")
            report_file = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"test_input", "Test_BioProject_report.xml"), 'wb')
            ftp.retrbinary('RETR report.xml', report_file.write, 262144)
            report_file.close()
            with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"test_input", "Test_BioProject_report.xml"), 'r') as file:
                for line in file:
                    if "Object accession=" in line:
                        print("Test BioProject ID is: " + line.split("Object accession=\"")[1].split("\" last_update=")[0])
        else:
            print("Submitting test BioProject.")
            res = ftp.storlines("STOR " + "submission.xml", open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"test_input", "submission.xml"), 'rb'))
            if not res.startswith('226 Transfer complete'):
                print('Submission.xml upload failed.', file=sys.stderr)
            res = ftp.storlines("STOR " + "submit.ready", open(submit_file_path, 'rb'))
            if not res.startswith('226 Transfer complete'):
                print('Submission.xml upload failed.', file=sys.stderr)
    except ftplib.all_errors as e:
        print('FTP error:', e)
        sys.exit(1)

def version():
    print("Version: " + _version)

def main():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,
        description='Command line upload script.')
    subparsers = parser.add_subparsers(dest='command')

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser_prep = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--unique_name',
        help='unique identifier',
        required=True
    )
    parent_parser.add_argument('--config',
        help='file path for config file for submission',
        required=False,
        default=(os.path.join(os.path.dirname(os.path.abspath(__file__)), "config_files", "default_config.yaml"))
    )
    parent_parser.add_argument('--upload_log',
        help='file path for upload log file',
        required=False,
        default=(os.path.join(os.path.dirname(os.path.abspath(__file__)), "upload_log.csv"))
    )
    parent_parser.add_argument('--test',
        help='Perform test submission.',
        required=False,
        default="Production",
        action="store_const",
        const="Test"
    )
    parent_parser.add_argument('--overwrite',
        help='Overwrite existing submission on NCBI',
        required=False,
        default=False,
        action="store_true"
    )
    parent_parser_prep.add_argument("--metadata",
        help="Metadata file",
        required=True
    )
    parent_parser_prep.add_argument("--fasta",
        help="Fasta file",
        required=True
    )

    runall_parser = subparsers.add_parser('submit',formatter_class=RawTextHelpFormatter,
        help='Creates submission files and begins automated process of submitting to public databases.',
        parents=[parent_parser, parent_parser_prep],
        description='Creates submission files and begins process of submitting to public databases.'
    )

    prep_parser = subparsers.add_parser('prep',formatter_class=RawTextHelpFormatter,
        help='Creates submission files.', parents=[parent_parser, parent_parser_prep],
        description='Creates submission files.'
    )

    log_parser = subparsers.add_parser('update_submissions',formatter_class=RawTextHelpFormatter,
        help='Using submission log, script updates existing process of submissions.',
        parents=[parent_parser],
        description='Using submission log, script updates existing process of submissions.'
    )

    genbank_parser = subparsers.add_parser('genbank',formatter_class=RawTextHelpFormatter,
        help='Performs manual submission to Genbank.', parents=[parent_parser],
        description='Performs submission to Genbank.'
    )

    biosample_parser = subparsers.add_parser('biosample',formatter_class=RawTextHelpFormatter,
        help='Performs manual submission to BioSample.', parents=[parent_parser],
        description='Performs submission to BioSample.'
    )

    biosample_sra_parser = subparsers.add_parser('biosample_sra',formatter_class=RawTextHelpFormatter,
        help='Performs manual joint submission to BioSample and SRA if enabled in config.', parents=[parent_parser],
        description='Performs joint submission to BioSample and SRA if enabled in config.'
    )

    sra_parser = subparsers.add_parser('sra',formatter_class=RawTextHelpFormatter,
        help='Performs manual submission to SRA.', parents=[parent_parser],
        description='Performs submission to SRA.'
    )

    gisaid_parser = subparsers.add_parser('gisaid',formatter_class=RawTextHelpFormatter,
        help='Performs manual submission to GISAID.', parents=[parent_parser],
        description='Performs submission to GISAID.'
    )

    test_bioproject_parser = subparsers.add_parser('test_bioproject',formatter_class=RawTextHelpFormatter,
        help='Create a test BioProject for making test submissions.', parents=[parent_parser],
        description='Create a test BioProject for making test submissions.'
    )

    test_bioproject_parser = subparsers.add_parser('version',formatter_class=RawTextHelpFormatter,
        help='Version info.',
        description='Version info.'
    )

    args = parser.parse_args()

    if args.command == 'submit':
        submission_preparation.process_submission(
            args.unique_name,
            args.fasta,
            args.metadata,
            args.config
        )
        start_submission(
            unique_name=args.unique_name,
            config=args.config,
            upload_log_path=args.upload_log,
            test=args.test,
            overwrite=args.overwrite
        )
    elif args.command == "prep":
        submission_preparation.process_submission(
            args.unique_name,
            args.fasta,
            args.metadata,
            args.config
        )
    elif args.command == "update_submissions":
        update_log(upload_log_path=args.upload_log)
    elif args.command == "genbank":
        submit_genbank(
            unique_name=args.unique_name,
            config=args.config,
            upload_log_path=args.upload_log,
            test=args.test,
            overwrite=args.overwrite
        )
    elif args.command == "gisaid":
        submit_gisaid(
            unique_name=args.unique_name,
            config=args.config,
            upload_log_path=args.upload_log,
            test=args.test
        )
    elif args.command == "biosample_sra":
        submit_biosample_sra(
            unique_name=args.unique_name,
            config=args.config,
            upload_log_path=args.upload_log,
            test=args.test,
            ncbi_sub_type="biosample_sra",
            overwrite=args.overwrite
        )
    elif args.command == "biosample":
        submit_biosample_sra(
            unique_name=args.unique_name,
            config=args.config,
            upload_log_path=args.upload_log,
            test=args.test,
            ncbi_sub_type="biosample",
            overwrite=args.overwrite
        )
    elif args.command == "sra":
        submit_biosample_sra(
            unique_name=args.unique_name,
            config=args.config,
            upload_log_path=args.upload_log,
            test=args.test,
            ncbi_sub_type="sra",
            overwrite=args.overwrite
        )
    elif args.command == "test_bioproject":
        test_bioproject(config=args.config)
    elif args.command == "version":
        version()
    else:
        print ("Invalid option")
        parser.print_help()

if __name__ == "__main__":
    main()
