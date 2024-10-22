from pandera import DataFrameSchema, Column, Check, Index, MultiIndex

schema = DataFrameSchema(
	columns={
		"((?:cmt-StructuredComment)(?:Prefix|Suffix))": Column(
			dtype="object",
			checks=[
				Check.str_length(max_value=50),
			],
			nullable=False,
			unique=False,
			coerce=False,
			required=True,
			regex=True,
			default="Assembly-Data",
			description="Structured comment keyword. For FLU use \"FluData\", HIV use \"HIV-DataBaseData\", and for COV and other organisms use \"Assembly-Data\".",
			title="Structured Comment Prefix and Suffix",
		),
		"Assembly Method": Column(
			dtype="object",
			checks=[
				Check.str_length(max_value=50),
			],
			nullable=False,
			unique=False,
			coerce=False,
			required=True,
			description="Structured comment keyword. For FLU use \"FluData\", HIV use \"HIV-DataBaseData\", and for COV and other organisms use \"Assembly-Data\".",
			title="Structured Comment Prefix and Suffix",
		),
	},
     checks=None,
     index=None,
     coerce=False,
     strict="filter",
     name="genbank_cmt_schema",
     ordered=False,
     unique=None,
     report_duplicates="all",
     unique_column_names=True,
     add_missing_columns=False,
     title="seqsender genbank comment file schema",
     description="Schema validation for GenBank database comment file.",
)
