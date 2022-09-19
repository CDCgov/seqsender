#!/usr/bin/env python3

import ftplib
import argparse
from argparse import RawTextHelpFormatter
import sys
from datetime import datetime
import requests
import json
import os
import pandas as pd
import yaml
from Bio import SeqIO
import xml.etree.ElementTree as ET
import shutil
import re

work_dir = os.getcwd()

def main():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,
        description='Check status of a submission')

    parser.add_argument("--unique_name",
        help="Unique identifier for a submission",
        required=False,
        default="")

    args = parser.parse_args()
    
    unique_name = args.unique_name

    if os.path.isfile(os.path.join(os.path.dirname(os.path.abspath(__file__)), "upload_log.csv")):
        df = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "upload_log.csv"), header = 0, dtype = str)
        df = df.fillna("")
    else:
        print("Error: Either a submission has not been made or upload_log.csv has been moved from current working directory", file=sys.stderr)
        sys.exit(1)

    if df["name"].str.contains(unique_name).any():
        df = df[df["name"] == unique_name]
        for index, row in df.iterrows():
            print("\n" + unique_name + " " + row["type"] + " status:")
            print("BioSample:\n\tSubmission id: " + row["BioSample_submission_id"]  + "\n\tStatus:" + row["BioSample_status"])
            print("SRA:\n\tSubmission id: " + row["SRA_submission_id"]  + "\n\tStatus:" + row["SRA_status"])
            print("Genbank:\n\tSubmission id: " + row["Genbank_submission_id"]  + "\n\tStatus:" + row["Genbank_status"])
            print("GISAID:\n\tSubmitted: " + row["GISAID_submitted_total"]  + "\n\tFailed: :" + row["GISAID_failed_total"])
        print("")
    else:
        print(unique_name + " not in log")

if __name__ == "__main__":
    main()
