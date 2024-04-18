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
		"bs-air_temp_regm": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="information about treatment involving an exposure to varying temperatures; should include the temperature, treatment duration, interval and total experimental duration; can include different temperature regimens",
			title="air temperature regimen",
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
		"bs-ances_data": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Information about either pedigree or other ancestral information description, e.g., parental variety in case of mutant or selection, A/3*B (meaning [(A x B) x B] x B)",
			title="ancestral data",
		),
		"bs-antibiotic_regm": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="information about treatment involving antibiotic administration; should include the name of antibiotic, amount administered, treatment duration, interval and total experimental duration; can include multiple antibiotic regimens",
			title="antibiotic regimen",
		),
		"bs-biol_stat": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="The level of genome modification, e.g., wild, natural, semi-natural, inbred line, breeder's line, hybrid, clonal selection, mutant",
			title="biological status",
		),
		"bs-biotic_regm": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Information about treatment(s) involving use of biotic factors, such as bacteria, viruses or fungi",
			title="biotic regimen",
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
		"bs-chem_mutagen": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="treatment involving use of mutagens; should include the name of mutagen, amount administered, treatment duration, interval and total experimental duration; can include multiple mutagen regimens",
			title="chemical mutagen",
		),
		"bs-climate_environment": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="treatment involving an exposure to a particular climate; can include multiple climates",
			title="climate environment",
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
		"bs-cult_root_med": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Name or reference for the hydroponic or in vitro culture rooting medium; can be the name of a commonly used medium or reference to a specific medium, e.g., Murashige and Skoog medium. If the medium has not been formally published, use the rooting medium descriptors",
			title="culture rooting medium",
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
		"bs-fertilizer_regm": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="information about treatment involving the use of fertilizers; should include the name fertilizer, amount administered, treatment duration, interval and total experimental duration; can include multiple fertilizer regimens",
			title="fertilizer regimen",
		),
		"bs-fungicide_regm": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="information about treatment involving use of fungicides; should include the name of fungicide, amount administered, treatment duration, interval and total experimental duration; can include multiple fungicide regimens",
			title="fungicide regimen",
		),
		"bs-gaseous_environment": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="use of conditions with differing gaseous environments; should include the name of gaseous compound, amount administered, treatment duration, interval and total experimental duration; can include multiple gaseous environment regimens",
			title="gaseous environment",
		),
		"bs-genetic_mod": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Genetic modifications of the genome of an organism, which may occur naturally by spontaneous mutation, or be introduced by some experimental means, e.g. specification of a transgene or the gene knocked-out or details of transient transfection",
			title="genetic modification",
		),
		"bs-gravity": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="information about treatment involving use of gravity factor to study various types of responses in presence, absence or modified levels of gravity; can include multiple treatments",
			title="gravity",
		),
		"bs-growth_facil": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Type of facility where the sampled plant was grown; controlled vocabulary: growth chamber, open top chamber, glasshouse, experimental garden, field. Alternatively use Crop Ontology (CO) terms, see http://www.cropontology.org/ontology/CO_715/Crop%20Research",
			title="growth facility",
		),
		"bs-growth_habit": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Characteristic shape, appearance or growth form of a plant species, e.g., erect, semi-erect, spreading, prostrate",
			title="growth habit",
		),
		"bs-growth_hormone_regm": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="information about treatment involving use of growth hormones; should include the name of growth hormone, amount administered, treatment duration, interval and total experimental duration; can include multiple growth hormone regimens",
			title="growth hormone regimen",
		),
		"bs-herbicide_regm": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="information about treatment involving use of herbicides; information about treatment involving use of growth hormones; should include the name of herbicide, amount administered, treatment duration, interval and total experimental duration; can include multiple regimens",
			title="herbicide regimen",
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
		"bs-host_disease": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Name of relevant disease, e.g. Salmonella gastroenteritis. Controlled vocabulary, http://bioportal.bioontology.org/ontologies/1009 or http://www.ncbi.nlm.nih.gov/mesh",
			title="host disease",
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
		"bs-host_genotype": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			title="host genotype",
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
		"bs-host_life_stage": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="description of host life stage",
			title="host life stage",
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
		"bs-host_subspecf_genlin": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Information about the genetic distinctness of the host organism below the subspecies level e.g., serovar, serotype, biotype, ecotype, variety, cultivar, or any relevant genetic typing schemes like Group I plasmid. Subspecies should not be recorded in this term, but in the NCBI taxonomy. Supply both the lineage name and the lineage rank separated by a colon, e.g., biovar:abc123",
			title="host subspecific genetic lineage",
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
		"bs-host_wet_mass": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="measurement of wet mass",
			title="host wet mass",
		),
		"bs-humidity_regm": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="information about treatment involving an exposure to varying degree of humidity; information about treatment involving use of growth hormones; should include amount of humidity administered, treatment duration, interval and total experimental duration; can include multiple regimens",
			title="humidity regimen",
		),
		"bs-light_regm": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Information about treatment(s) involving exposure to light, including both light intensity and quality",
			title="light regimen",
		),
		"bs-mechanical_damage": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="information about any mechanical damage exerted on the plant; can include multiple damages and sites",
			title="mechanical damage",
		),
		"bs-mineral_nutr_regm": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="information about treatment involving the use of mineral supplements; should include the name of mineral nutrient, amount administered, treatment duration, interval and total experimental duration; can include multiple mineral nutrient regimens",
			title="mineral nutrient regimen",
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
		"bs-non_mineral_nutr_regm": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="information about treatment involving the exposure of plant to non-mineral nutrient such as oxygen, hydrogen or carbon; should include the name of non-mineral nutrient, amount administered, treatment duration, interval and total experimental duration; can include multiple non-mineral nutrient regimens",
			title="non mineral nutrient regimen",
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
		"bs-pesticide_regm": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="information about treatment involving use of insecticides; should include the name of pesticide, amount administered, treatment duration, interval and total experimental duration; can include multiple pesticide regimens",
			title="pesticide regimen",
		),
		"bs-ph_regm": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="information about treatment involving exposure of plants to varying levels of pH of the growth media; can include multiple regimen",
			title="pH regimen",
		),
		"bs-plant_growth_med": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Type of the media used for growing sampled plants, e.g., soil [ENVO:00001998]",
			title="plant growth medium",
		),
		"bs-plant_product": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="substance produced by the plant, where the sample was obtained from",
			title="plant product",
		),
		"bs-plant_sex": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Sex of the reproductive parts on the whole plant, e.g., androdioecious, androecious, androgynous, androgynomonoecious, andromonoecious, bisexual, dichogamous, diclinous, dioecious, gynodioecious, gynoecious, gynomonoecious, hermaphroditic, imperfect, monoclinous, monoecious, perfect, polygamodioecious, polygamomonoecious, polygamous, protandrous, protogynous, subandroecious, subdioecious, subgynoecious, synoecious, trimonoecious, trioecious, unisexual",
			title="plant sex",
		),
		"bs-plant_struc": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Name of plant structure the sample was obtained from; for Plant Ontology (PO) (v releases/2017-12-14) terms, see http://purl.bioontology.org/ontology/PO, e.g., petiole epidermis (PO_0000051). If an individual flower is sampled, the sex of it can be recorded here",
			title="plant structure",
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
		"bs-radiation_regm": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="information about treatment involving exposure of plant or a plant part to a particular radiation regimen; should include the radiation type, amount or intensity administered, treatment duration, interval and total experimental duration; can include multiple radiation regimens",
			title="radiation regimen",
		),
		"bs-rainfall_regm": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="information about treatment involving an exposure to a given amount of rainfall; can include multiple regimens",
			title="rainfall regimen",
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
		"bs-root_cond": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Relevant rooting conditions such as field plot size, sowing density, container dimensions, number of plants per container",
			title="rooting conditions",
		),
		"bs-root_med_carbon": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Source of organic carbon in the culture rooting medium, e.g., sucrose",
			title="rooting medium carbon",
		),
		"bs-root_med_macronutr": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Measurement of the culture rooting medium macronutrients (N,P, K, Ca, Mg, S), e.g., KH2PO4 (170mg/L)",
			title="rooting medium macronutrients",
		),
		"bs-root_med_micronutr": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Measurement of the culture rooting medium micronutrients (Fe, Mn, Zn, B, Cu, Mo), e.g., H3BO3 (6.2mg/L)",
			title="rooting medium micronutrients",
		),
		"bs-root_med_ph": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="pH measurement of the culture rooting medium, e.g., 5.5",
			title="rooting medium pH",
		),
		"bs-root_med_regl": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Growth regulators in the culture rooting medium such as cytokinins, auxins, gybberellins, abscisic acid, e.g., 0.5 mg/L NAA",
			title="rooting medium regulators",
		),
		"bs-root_med_solid": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Specification of the solidifying agent in the culture rooting medium, e.g., agar",
			title="rooting medium solidifier",
		),
		"bs-root_med_suppl": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Organic supplements of the culture rooting medium, such as vitamins, amino acids, organic acids, antibiotics activated charcoal, e.g., nicotinic acid (0.5 mg/L)",
			title="rooting medium organic supplements",
		),
		"bs-salt_regm": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="information about treatment involving use of salts as supplement to liquid and soil growth media; should include the name of salt, amount administered, treatment duration, interval and total experimental duration; can include multiple salt regimens",
			title="salt regimen",
		),
		"bs-samp_capt_status": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Reason for the sample, e.g., active surveillance in response to an outbreak, active surveillance not initiated by an outbreak, farm sample, market sample",
			title="sample capture status",
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
		"bs-samp_dis_stage": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="Stage of the disease at the time of sample collection, e.g., dissemination, growth and reproduction, infection, inoculation, penetration",
			title="sample disease stage",
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
		"bs-season_environment": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="treatment involving an exposure to a particular season (e.g. winter, summer, rabi, rainy etc.)",
			title="seasonal environment",
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
		"bs-standing_water_regm": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="treatment involving an exposure to standing water during a plant's life span, types can be flood water or standing water; can include multiple regimens",
			title="standing water regimen",
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
		"bs-tiss_cult_growth_med": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="description of plant tissue culture growth media used",
			title="tissue culture growth media",
		),
		"bs-water_temp_regm": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="information about treatment involving an exposure to water with varying degree of temperature; can include multiple regimens",
			title="water temperature regimen",
		),
		"bs-watering_regm": Column(
			dtype="object",
			checks=None,
			nullable=True,
			unique=False,
			coerce=False,
			required=False,
			description="information about treatment involving an exposure to watering frequencies; can include multiple regimens",
			title="watering regimen",
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
	name="biosample_package_MISAG.plant-associated.6.0_schema",
	ordered=False,
	unique=None,
	report_duplicates="all",
	unique_column_names=True,
	add_missing_columns=True,
	title="BioSample package MISAG.plant-associated.6.0 schema",
	description="Schema validation for BioSample database using MISAG.plant-associated.6.0 package.",
)