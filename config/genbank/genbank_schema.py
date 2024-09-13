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
			description="Sequence identifier used in fasta file. This is used to create the fasta file for Genbank or GISAID by updating the sequence name in your fasta file to reflect the sample name for the specified database.",
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
			description="Identifier name used for GenBank. Max length is 50 characters.",
			title="sample name",
		),
		"gb-fasta_definition_line_modifiers": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="NCBI fasta definition line modifiers can be added here. As many modifiers as you like can be added, but each must bounded by a set of brackets. Some of the available keywords are listed at \"https://www.ncbi.nlm.nih.gov/genbank/mods_fastadefline/\".",
			title="fasta definition line modifiers",
		),
		"gb-title": Column(
			dtype="object",
			checks=[
				Check(lambda s: s.nunique() == 1),
			],
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Optional internal field for how the GenBank submission should be named when viewed from the NCBI submission portal, . If not provided, when performing submissions <--submission_name> with the suffix \"-GB\" will be used instead.",
			title="genbank submission portal title",
		),
		"sra-comment": Column(
			dtype="object",
			checks=[
				Check(lambda s: s.nunique() == 1),
			],
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Optional internal field explaining the purpose of the submission for when interacting and resolving submission issues with NCBI.",
			title="genbank submission portal description",
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
