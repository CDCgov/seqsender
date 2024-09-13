from pandera import DataFrameSchema, Column, Check, Index, MultiIndex
import pandas as pd
import re

# Validate date format of "YYYY", "YYYY-MM", or "YYYY-MM-DD", then validate that date is a valid date
def validate_date(date):
	if not date or date.strip() == "":
		return False
	elif re.match(r"(?i)(\W|^)(\d{4}|\d{4}-\d{2}|\d{4}-\d{2}-\d{2})(\W|$)", date):
		for date_format in ["%Y", "%Y-%m", "%Y-%m-%d"]:
			try:
				pd.to_datetime(date, format=date_format, errors="raise")
				return True
			except ValueError:
				pass
	return False

schema = DataFrameSchema(
	columns={
		"organism": Column(
			dtype="object",
			checks=[
				Check.str_matches(r"^(?!\s*$).+"),
			],
			nullable=False,
			unique=False,
			coerce=False,
			required=True,
			description="The most descriptive organism name for the samples.",
			title="organism",
		),
		"authors": Column(
			dtype="object",
			checks=[
				Check.str_matches(r"^(?!\s*$).+"),
			],
			nullable=False,
			unique=False,
			coerce=False,
			required=True,
			description="Citing authors. List of <b>Last, First Middle, suffix</b> separated by a semicolon \";\" E.g.: \"Baker, Howard Henry, Jr.; Powell, Earl Alexander, III.;\"",
			title="authors",
		),
		"collection_date": Column(
			dtype="object",
			checks=[
				Check(validate_date, element_wise=True, error="invalid_date_format")
			],
			nullable=False,
			unique=False,
			coerce=False,
			required=True,
			description="Collection date for sample. Must be in a valid format based on ISO 8601: \"YYYY-MM-DD\", \"YYYY-MM\", or \"YYYY\". Time is not supported.",
			title="collection date",
		),
		"bioproject": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="NCBI BioProject ID.",
			title="bioproject",
		),
	},
     checks=None,
     index=None,
     coerce=False,
     strict="filter",
     name="seqsender_schema",
     ordered=False,
     unique=None,
     report_duplicates="all",
     unique_column_names=True,
     add_missing_columns=False,
     title="seqsender pipeline schema",
     description="Schema validation for the base requirements to operate Seqsender for all databases.",
)
