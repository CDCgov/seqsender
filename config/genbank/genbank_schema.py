from pandera import DataFrameSchema, Column, Check, Index, MultiIndex

schema = DataFrameSchema(
	columns={
		"sequence_name": Column(
			dtype="object",
			checks=[
				Check.str_matches(r"^(?!\s*$).+"),
			],
			nullable=False,
			unique=True,
			coerce=False,
			required=True,
			description="Sequence identifier used in fasta file. This is used to create the fasta file for Genbank and/or GISAID.",
			title="sequence name",
		),
		"gb-sample_name": Column(
			dtype="object",
			checks=[
				Check.str_matches(r"^(?!\s*$).+"),
				Check.str_length(max_value=50),
			],
			nullable=False,
			unique=True,
			coerce=False,
			required=True,
			description="Identifier name used for GenBank. Max length is 50 characters. Fasta modifiers with brackets \"[]\" can be added. They will be added only to the fasta file.",
			title="sample name",
		)
	},
     checks=None,
     index=None,
     coerce=False,
     strict="filter",
     name="genbank_schema",
     ordered=False,
     unique=None,
     report_duplicates="all",
     unique_column_names=True,
     add_missing_columns=False,
     title="seqsender genbank schema",
     description="Schema validation for GenBank database.",
)
