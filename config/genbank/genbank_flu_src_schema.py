from pandera import DataFrameSchema, Column, Check, Index, MultiIndex

schema = DataFrameSchema(
	columns={
		"src-Altitude": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Altitude in metres above or below sea level of where the sample was collected.",
			title="Altitude",
		),
		"src-Authority": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="The author or authors of the organism name from which sequence was obtained.",
			title="Authority",
		),
		"src-Bio_material": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="An identifier for the biological material from which the nucleotide sequence was obtained, with optional institution code and collection code for the place where it is currently stored. This should be provided using the following format 'institution-code:collection-code:material_id'. material_id is mandatory, institution-code and collection-code are optional; institution-code is mandatory when collection-code is present. This qualifier should be used to annotate the identifiers of material in biological collections which include zoos and aquaria, stock centers, seed banks, germplasm repositories and DNA banks.",
			title="Bio_material",
		),
		"src-Biotype": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Variety of a species (usually a fungus, bacteria, or virus) characterized by some specific biological property (often geographical, ecological, or physiological). Same as biotype.",
			title="Biotype",
		),
		"src-Biovar": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="See biotype",
			title="Biovar",
		),
		"src-Breed": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="The named breed from which sequence was obtained (usually applied to domesticated mammals).",
			title="Breed",
		),
		"src-Cell_line": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Cell line from which sequence was obtained.",
			title="Cell_line",
		),
		"src-Cell_type": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Type of cell from which sequence was obtained.",
			title="Cell_type",
		),
		"src-Chemovar": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Variety of a species (usually a fungus, bacteria, or virus) characterized by its biochemical properties.",
			title="Chemovar",
		),
		"src-Clone": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Name of clone from which sequence was obtained.",
			title="Clone",
		),
		"src-Collected_by": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Name of person who collected the sample.",
			title="Collected_by",
		),
		"src-Country": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="The country where the sequence's organism was located. May also be an ocean or major sea. Additional region or locality information must be after the country name and separated by a ':'. For example: USA: Riverview Park, Ripkentown, MD",
			title="Country",
		),
		"src-Cultivar": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Cultivated variety of plant from which sequence was obtained.",
			title="Cultivar",
		),
		"src-Culture_collection": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Institution code and identifier for the culture from which the nucleotide sequence was obtained, with optional collection code. This should be provided using the following format 'institution-code:collection-code:culture-id'. culture-id and institution-code are mandatory. This qualifier should be used to annotate live microbial and viral cultures, and cell lines that have been deposited in curated culture collections.",
			title="Culture_collection",
		),
		"src-Dev_stage": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Developmental stage of organism.",
			title="Dev_stage",
		),
		"src-Ecotype": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="The named ecotype (population adapted to a local habitat) from which sequence was obtained (customarily applied to populations of Arabidopsis thaliana).",
			title="Ecotype",
		),
		"src-Forma": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="The forma (lowest taxonomic unit governed by the nomenclatural codes) of organism from which sequence was obtained. This term is usually applied to plants and fungi.",
			title="Forma",
		),
		"src-Forma_specialis": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="The physiologically distinct form from which sequence was obtained (usually restricted to certain parasitic fungi).",
			title="Forma_specialis",
		),
		"src-Fwd_primer_name": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="name of forward PCR primer",
			title="Fwd_primer_name",
		),
		"src-Fwd_primer_seq": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="nucleotide sequence of forward PCR primer",
			title="Fwd_primer_seq",
		),
		"src-Genotype": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Genotype of the organism.",
			title="Genotype",
		),
		"src-Haplogroup": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Name for a group of similar haplotypes that share some sequence variation",
			title="Haplogroup",
		),
		"src-Haplotype": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Haplotype of the organism.",
			title="Haplotype",
		),
		"src-Host": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="When the sequence submission is from an organism that exists in a symbiotic, parasitic, or other special relationship with some second organism, the 'host' modifier can be used to identify the name of the host species.",
			title="Host",
		),
		"src-Identified_by": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="name of the person or persons who identified by taxonomic name the organism from which the sequence was obtained",
			title="Identified_by",
		),
		"src-Isolate": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Identification or description of the specific individual from which this sequence was obtained.",
			title="Isolate",
		),
		"src-Isolation-source": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Describes the local geographical source of the organism from which the sequence was obtained.",
			title="Isolation source",
		),
		"src-Lab_host": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Laboratory host used to propagate the organism from which the sequence was obtained.",
			title="Lab_host",
		),
		"src-Lat_Lon": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Latitude and longitude, in decimal degrees, of where the sample was collected.",
			title="Lat_Lon",
		),
		"src-Note": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Any additional information that you wish to provide about the sequence.",
			title="Note",
		),
		"src-Pathovar": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Variety of a species (usually a fungus, bacteria or virus) characterized by the biological target of the pathogen. Examples include Pseudomonas syringae pathovar tomato and Pseudomonas syringae pathovar tabaci.",
			title="Pathovar",
		),
		"src-Pop_variant": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="name of the population variant from which the sequence was obtained",
			title="Pop_variant",
		),
		"src-Rev_primer_name": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="name of reverse PCR primer",
			title="Rev_primer_name",
		),
		"src-Rev_primer_seq": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="nucleotide sequence of reverse PCR primer",
			title="Rev_primer_seq",
		),
		"src-Segment": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="name of viral or phage segment sequenced",
			title="Segment",
		),
		"src-Serogroup": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Variety of a species (usually a fungus, bacteria, or virus) characterized by its antigenic properties. Same as serogroup and serovar.",
			title="Serogroup",
		),
		"src-Serotype": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="See Serogroup",
			title="Serotype",
		),
		"src-Serovar": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="See Serogroup",
			title="Serovar",
		),
		"src-Sex": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Sex of the organism from which the sequence was obtained.",
			title="Sex",
		),
		"src-Specimen_voucher": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="An identifier of the individual or collection of the source organism and the place where it is currently stored, usually an institution. This should be provided using the following format 'institution-code:collection-code:specimen-id'. specimen-id is mandatory, collection-code is optional; institution-code is mandatory when collection-code is provided. Examples: 99-SRNP UAM:Mamm:52179 personal collection:Joe Smith:99-SRNP AMCC:101706",
			title="Specimen_voucher",
		),
		"src-Strain": Column(
			dtype="object",
			checks=None,
			nullable=False,
			unique=False,
			coerce=False,
			required=True,
			description="Strain of organism from which sequence was obtained.",
			title="Strain",
		),
		"src-Sub_species": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Subspecies of organism from which sequence was obtained.",
			title="Sub_species",
		),
		"src-Subclone": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Name of subclone from which sequence was obtained.",
			title="Subclone",
		),
		"src-Subtype": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Subtype of organism from which sequence was obtained.",
			title="Subtype",
		),
		"src-Substrain": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Sub-strain of organism from which sequence was obtained.",
			title="Substrain",
		),
		"src-Tissue_lib": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Tissue library from which the sequence was obtained.",
			title="Tissue_lib",
		),
		"src-Tissue_type": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Type of tissue from which sequence was obtained.",
			title="Tissue_type",
		),
		"src-Type": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Type of organism from which sequence was obtained.",
			title="Type",
		),
		"src-Variety": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Variety of organism from which sequence was obtained.",
			title="Variety",
		),
	},
     checks=None,
     index=None,
     coerce=False,
     strict="filter",
     name="genbank_src_schema",
     ordered=False,
     unique=None,
     report_duplicates="all",
     unique_column_names=True,
     add_missing_columns=False,
     title="seqsender genbank source file schema",
     description="Schema validation for GenBank database source file.",
)