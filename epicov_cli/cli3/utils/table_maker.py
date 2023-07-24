"""
    This module does the heavy lifting of building the tables.
    Copyright (C) 2022 Freunde von GISAID e.V.
"""
import pandas as pd

class Table():
    def __init__(self, indata):
        self.indata = indata

    def gisaid_template(self, unknown):
        df = None
        try:
            df = pd.read_excel(self.indata, header=0, sheet_name=1)
        except:
            df = pd.read_csv(self.indata, header=0, sep=",")
        df.set_index("covv_virus_name", inplace=True)
        df = df.drop(df.index[0]) #drops the first row (=duplicated header) only if required ****
        df.replace("unknown", unknown, inplace=True)
        return df