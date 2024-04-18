from pandera import DataFrameSchema, Column, Check, Index, MultiIndex

schema = DataFrameSchema(
	columns={
		"bs-isolate": Column(
			dtype="object",
			checks=None,
			nullable=False,
			unique=False,
			coerce=False,
			required=True,
			description="identification or description of the specific individual from which this sample was obtained",
			title="isolate",
		),
		"bs-env_broad_scale": Column(
			dtype="object",
			checks=None,
			nullable=False,
			unique=False,
			coerce=False,
			required=True,
			description="Add terms that identify the major environment type(s) where your sample was collected. Recommend subclasses of biome [ENVO:00000428]. Multiple terms can be separated by one or more pipes e.g.: Â mangrove biome [ENVO:01000181]|estuarine biome [ENVO:01000020]",
			title="broad-scale environmental context",
		),
		"bs-env_local_scale": Column(
			dtype="object",
			checks=None,
			nullable=False,
			unique=False,
			coerce=False,
			required=True,
			description="Add terms that identify environmental entities having causal influences upon the entity at time of sampling, multiple terms can be separated by pipes, e.g.: Â shoreline [ENVO:00000486]|intertidal zone [ENVO:00000316]",
			title="local-scale environmental context",
		),
		"bs-env_medium": Column(
			dtype="object",
			checks=None,
			nullable=False,
			unique=False,
			coerce=False,
			required=True,
			description="Add terms that identify the material displaced by the entity at time of sampling. Recommend subclasses of environmental material [ENVO:00010483]. Multiple terms can be separated by pipes e.g.: estuarine water [ENVO:01000301]|estuarine mud [ENVO:00002160]",
			title="environmental medium",
		),
		"bs-geo_loc_name": Column(
			dtype="object",
			checks=None,
			nullable=False,
			unique=False,
			coerce=False,
			required=True,
			description="Geographical origin of the sample; use the appropriate name from this list http://www.insdc.org/documents/country-qualifier-vocabulary. Use a colon to separate the country or ocean from more detailed information about the location, eg \"Canada: Vancouver\" or \"Germany: halfway down Zugspitze, Alps\"",
			title="geographic location",
		),
		"bs-host_dependence": Column(
			dtype="object",
			checks=None,
			nullable=False,
			unique=False,
			coerce=False,
			required=True,
			description="Type of host dependence for the symbiotic host organism to its host., e.g., facultative, obligate",
			title="host dependence",
		),
		"bs-host_life_stage": Column(
			dtype="object",
			checks=None,
			nullable=False,
			unique=False,
			coerce=False,
			required=True,
			description="description of host life stage",
			title="host life stage",
		),
		"bs-isolation_source": Column(
			dtype="object",
			checks=None,
			nullable=False,
			unique=False,
			coerce=False,
			required=True,
			description="Describes the physical, environmental and/or local geographical source of the biological sample from which the sample was derived.",
			title="isolation source",
		),
		"bs-lat_lon": Column(
			dtype="object",
			checks=None,
			nullable=False,
			unique=False,
			coerce=False,
			required=True,
			description="The geographical coordinates of the location where the sample was collected. Specify as degrees latitude and longitude in format \"d[d.dddd] N|S d[dd.dddd] W|E\", eg, 38.98 N 77.11 W",
			title="latitude and longitude",
		),
		"bs-sym_life_cycle_type": Column(
			dtype="object",
			checks=None,
			nullable=False,
			unique=False,
			coerce=False,
			required=True,
			description="Type of life cycle of the symbiotic host species (the thing being sampled). Simple life cycles occur within a single host, complex ones within multiple different hosts over the course of their normal life cycle, e.g., complex life cycle, simple life cycle",
			title="symbiotic host organism life cycle type",
		),
		"bs-altitude": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="The altitude of the sample is the vertical distance between Earth's surface above Sea Level and the sampled position in the air.",
			title="altitude",
		),
		"bs-association_duration": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Time spent in host of the symbiotic organism at the time of sampling; relevant scale depends on symbiotic organism and study",
			title="duration of association with the host",
		),
		"bs-chem_administration": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="list of chemical compounds administered to the host or site where sampling occurred, and when (e.g. antibiotics, N fertilizer, air filter); can include multiple compounds. For Chemical Entities of Biological Interest ontology (CHEBI) (v1.72), please see http://bioportal.bioontology.org/visualize/44603",
			title="chemical administration",
		),
		"bs-collection_method": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Process used to collect the sample, e.g., bronchoalveolar lavage (BAL)",
			title="collection method",
		),
		"bs-depth": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Depth is defined as the vertical distance below surface, e.g. for sediment or soil samples depth is measured from sediment or soil surface, respectivly. Depth can be reported as an interval for subsurface samples.",
			title="depth",
		),
		"bs-elev": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="The elevation of the sampling site as measured by the vertical distance from mean sea level.",
			title="elevation",
		),
		"bs-experimental_factor": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Variable aspect of experimental design",
			title="experimental factor",
		),
		"bs-gravidity": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="whether or not subject is gravid, and if yes date due or date post-conception, specifying which is used",
			title="gravidity",
		),
		"bs-host_age": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Age of host at the time of sampling",
			title="host age",
		),
		"bs-host_body_habitat": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="original body habitat where the sample was obtained from",
			title="host body habitat",
		),
		"bs-host_body_product": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="substance produced by the host, e.g. stool, mucus, where the sample was obtained from",
			title="host body product",
		),
		"bs-host_cellular_loc": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="The localization of the symbiotic host organism within the host from which it was sampled, e.g., intracellular, extracellular",
			title="host cellular location",
		),
		"bs-host_color": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="the color of host",
			title="host color",
		),
		"bs-host_common_name": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="The natural language (non-taxonomic) name of the host organism, e.g., mouse",
			title="host common name",
		),
		"bs-host_dry_mass": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="measurement of dry mass",
			title="host dry mass",
		),
		"bs-host_family_relationship": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			title="host family relationship",
		),
		"bs-host_genotype": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			title="host genotype",
		),
		"bs-host_growth_cond": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="literature reference giving growth conditions of the host",
			title="host growth conditions",
		),
		"bs-host_height": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="the height of subject",
			title="host height",
		),
		"bs-host_infra_specific_name": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="taxonomic information subspecies level",
			title="host infra specific name",
		),
		"bs-host_infra_specific_rank": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="taxonomic rank information below subspecies level, such as variety, form, rank etc.",
			title="host infra specific rank",
		),
		"bs-host_length": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="the length of subject",
			title="host length",
		),
		"bs-host_number": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Number of symbiotic host individuals pooled at the time of collection",
			title="host number individual",
		),
		"bs-host_of_host_coinf": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="The taxonomic name of any coinfecting organism observed in a symbiotic relationship with the host of the sampled host organism, e.g. where a sample collected from a host trematode species (A) which was collected from a host_of_host fish (B) that was also infected with a nematode (C), the value here would be (C) the nematode {species name} or {common name}. Multiple co-infecting species may be added in a comma-separated list. For listing symbiotic organisms associated with the host (A) use the term Observed host symbiont",
			title="observed coinfecting organisms in host of host",
		),
		"bs-host_of_host_disease": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="List of diseases with which the host of the symbiotic host organism has been diagnosed; can include multiple diagnoses. The value of the field depends on host; for humans the terms should be chosen from the DO (Human Disease Ontology) at https://www.disease-ontology.org, non-human host diseases are free text",
			title="host of the symbiotic host disease status",
		),
		"bs-host_of_host_env_loc": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="For a symbiotic host organism the local anatomical environment within its host may have causal influences. Report the anatomical entity(s) which are in the direct environment of the symbiotic host organism being sampled and which you believe have significant causal influences on your sample or specimen. For example, if the symbiotic host organism being sampled is an intestinal worm, its local environmental context will be the term for intestine from UBERON (http://uberon.github.io/)",
			title="host of the symbiotic host local environmental context",
		),
		"bs-host_of_host_env_med": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Report the environmental material(s) immediately surrounding the symbiotic host organism at the time of sampling. This usually will be a tissue or substance type from the host, but may be another material if the symbiont is external to the host. We recommend using classes from the UBERON ontology, but subclasses of 'environmental material' (http://purl.obolibrary.org/obo/ENVO_00010483) may also be used. EnvO documentation about how to use the field: https://github.com/EnvironmentOntology/envo/wiki/Using-ENVO-with-MIxS. 
Terms from other OBO ontologies are permissible as long as they reference mass/volume nouns (e.g., air, water, blood) and not discrete, countable entities (e.g., intestines, heart)",
			title="host of the symbiotic host environemental medium",
		),
		"bs-host_of_host_fam_rel": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Familial relationship of the host of the symbiotic host organisms to other hosts of symbiotic host organism in the same study; can include multiple relationships",
			title="host of the symbiotic host family relationship",
		),
		"bs-host_of_host_geno": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Observed genotype of the host of the symbiotic host organism",
			title="host of the symbiotic host genotype",
		),
		"bs-host_of_host_gravid": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Whether or not the host of the symbiotic host organism is gravid, and if yes date due or date post-conception, specifying which is used",
			title="host of the symbiotic host gravidity",
		),
		"bs-host_of_host_infname": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Taxonomic name information of the host of the symbiotic host organism below subspecies level",
			title="host of the symbiotic host infra-specific name",
		),
		"bs-host_of_host_infrank": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Taxonomic rank information about the host of the symbiotic host organism below subspecies level, such as variety, form, rank",
			title="host of the symbiotic host infra-specific rank",
		),
		"bs-host_of_host_name": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Common name of the host of the symbiotic host organism",
			title="host of the symbiotic host common name",
		),
		"bs-host_of_host_pheno": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Phenotype of the host of the symbiotic host organism. For phenotypic quality ontology (PATO) terms, see http://purl.bioontology.org/ontology/pato",
			title="host of the symbiotic host phenotype",
		),
		"bs-host_of_host_sub_id": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="A unique identifier by which each host of the symbiotic host organism subject can be referred to, de-identified, e.g. #H14",
			title="host of the symbiotic host subject id",
		),
		"bs-host_of_host_taxid": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="NCBI taxon id of the host of the symbiotic host organism",
			title="host of the symbiotic host taxon id",
		),
		"bs-host_of_host_totmass": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Total mass of the host of the symbiotic host organism at collection, the unit depends on the host",
			title="host of the symbiotic host total mass",
		),
		"bs-host_phenotype": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			title="host phenotype",
		),
		"bs-host_sex": Column(
			dtype="object",
			checks=Check.str_matches(r"(?i)(\W|^)(male|female|pooled male and female|neuter|hermaphrodite|intersex|not determined|missing|not applicable|not collected|not provided|restricted access)(\W|$)"),
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Gender or physical sex of the host",
			title="host sex",
		),
		"bs-host_shape": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="morphological shape of host",
			title="host shape",
		),
		"bs-host_specificity": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Level of specificity of symbiont-host interaction, e.g., family-specific, generalist, genus-specific, species-specific",
			title="host specificity",
		),
		"bs-host_subject_id": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="a unique identifier by which each subject can be referred to, de-identified, e.g. #131",
			title="host subject id",
		),
		"bs-host_substrate": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="the growth substrate of the host",
			title="host substrate",
		),
		"bs-host_symbiont": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="The taxonomic name of the organism(s) found living in mutualistic, commensalistic, or parasitic symbiosis with the specific host",
			title="observed host symbionts",
		),
		"bs-host_taxid": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="NCBI taxonomy ID of the host, e.g. 9606",
			title="host taxonomy ID",
		),
		"bs-host_tissue_sampled": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="name of body site where the sample was obtained from, such as a specific organ or tissue, e.g., tongue, lung. For foundational model of anatomy ontology (fma) (v 4.11.0) or Uber-anatomy ontology (UBERON) (v releases/2014-06-15) terms, please see http://purl.bioontology.org/ontology/FMA or http://purl.bioontology.org/ontology/UBERON",
			title="host tissue sampled",
		),
		"bs-host_tot_mass": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="total mass of the host at collection, the unit depends on host",
			title="host total mass",
		),
		"bs-misc_param": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="any other measurement performed or parameter collected, that is not listed here",
			title="miscellaneous parameter",
		),
		"bs-mode_transmission": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="The process through which the symbiotic host organism entered the host from which it was sampled, e.g., horizontal:castrator, horizontal:directly transmitted, horizontal:micropredator, horizontal:parasitoid, horizontal:trophically transmitted, horizontal:vector transmitted, vertical",
			title="mode of transmission",
		),
		"bs-neg_cont_type": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="The substance or equipment used as a negative control in an investigation, e.g., distilled water, phosphate buffer, empty collection device, empty collection tube, DNA-free PCR mix, sterile swab, sterile syringe",
			title="negative control type",
		),
		"bs-omics_observ_id": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="A unique identifier of the omics-enabled observatory (or comparable time series) your data derives from. This identifier should be provided by the OMICON ontology; if you require a new identifier for your time series, contact the ontology's developers. Information is available here:Â https://github.com/GLOMICON/omicon. This field is only applicable to records which derive from an omics time-series or observatory.",
			title="Omics Observatory ID",
		),
		"bs-organism_count": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="total count of any organism per gram or volume of sample,should include name of organism followed by count; can include multiple organism counts",
			title="organism count",
		),
		"bs-oxy_stat_samp": Column(
			dtype="object",
			checks=Check.str_matches(r"(?i)(\W|^)(aerobic|anaerobic)(\W|$)"),
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="oxygenation status of sample",
			title="oxygenation status of sample",
		),
		"bs-perturbation": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="type of perturbation, e.g. chemical administration, physical disturbance, etc., coupled with time that perturbation occurred; can include multiple perturbation types",
			title="perturbation",
		),
		"bs-pos_cont_type": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="The substance, mixture, product, or apparatus used to verify that a process which is part of an investigation delivers a true positive",
			title="positive control type",
		),
		"bs-ref_biomaterial": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Primary publication or genome report",
			title="reference for biomaterial",
		),
		"bs-rel_to_oxygen": Column(
			dtype="object",
			checks=Check.str_matches(r"(?i)(\W|^)(aerobe|anaerobe|facultative|microaerophilic|microanaerobe|obligate aerobe|obligate anaerobe|missing|not applicable|not collected|not provided|restricted access)(\W|$)"),
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Is this organism an aerobe, anaerobe? Please note that aerobic and anaerobic are valid descriptors for microbial environments, eg, aerobe, anaerobe, facultative, microaerophilic, microanaerobe, obligate aerobe, obligate anaerobe, missing, not applicable, not collected, not provided, restricted access",
			title="relationship to oxygen",
		),
		"bs-route_transmission": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Description of path taken by the symbiotic host organism being sampled in order to establish a symbiotic relationship with the host (with which it was observed at the time of sampling) via a mode of transmission (specified in mode_transmission), e.g., environmental:faecal-oral, transplacental, vector-borne:vector penetration",
			title="route of transmission",
		),
		"bs-samp_collect_device": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Method or device employed for collecting sample",
			title="sample collection device or method",
		),
		"bs-samp_mat_process": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Processing applied to the sample during or after isolation",
			title="sample material processing",
		),
		"bs-samp_salinity": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			title="sample salinity",
		),
		"bs-samp_size": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Amount or size of sample (volume, mass or area) that was collected",
			title="sample size",
		),
		"bs-samp_store_dur": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			title="sample storage duration",
		),
		"bs-samp_store_loc": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			title="sample storage location",
		),
		"bs-samp_store_sol": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Solution within which sample was stored, if any",
			title="sample storage solution",
		),
		"bs-samp_store_temp": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			title="sample storage temperature",
		),
		"bs-samp_vol_we_dna_ext": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="volume (mL) or weight (g) of sample processed for DNA extraction",
			title="sample volume or weight for DNA extraction",
		),
		"bs-size_frac": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Filtering pore size used in sample preparation, e.g., 0-0.22 micrometer",
			title="size fraction selected",
		),
		"bs-source_material_id": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="unique identifier assigned to a material sample used for extracting nucleic acids, and subsequent sequencing. The identifier can refer either to the original material collected or to any derived sub-samples.",
			title="source material identifiers",
		),
		"bs-symbiont_host_role": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Role of the host in the life cycle of the symbiotic organism, e.g., accidental, dead-end, definitive, intermediate, paratenic, reservoir, single host",
			title="host of the symbiont role",
		),
		"bs-temp": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="temperature of the sample at time of sampling",
			title="temperature",
		),
		"bs-type_of_symbiosis": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Type of biological interaction established between the symbiotic host organism being sampled and its respective host, e.g., commensalistic, mutualistic, parasitic",
			title="type of symbiosis",
		),
		"bs-wga_amp_appr": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="method used to amplify genomic DNA in preparation for sequencing. Examples: MDR, WGA-X, MDA",
			title="WGA amplification approach",
		),
	},
	checks=None,
	index=None,
	dtype=None,
	coerce=False,
	strict="filter",
	name="biosample_package_MISAG.symbiont-associated.6.0_schema",
	ordered=False,
	unique=None,
	report_duplicates="all",
	unique_column_names=True,
	add_missing_columns=True,
	title="BioSample package MISAG.symbiont-associated.6.0 schema",
	description="Schema validation for BioSample database using MISAG.symbiont-associated.6.0 package.",
)