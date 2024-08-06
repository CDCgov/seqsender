from shiny import App, render, ui

test_data = [
    ui.h4("Generate test data for usage with SeqSender"),
    ui.p("Test data for each database can be generated for multiple available organisms. ",
         "To generate test data, use the ", ui.code(ui.strong("test_data")), " command. "
    ),
    ui.p("Based on the selected ", ui.code(ui.strong("databases")), " and ", ui.code(ui.strong("organism")),
         ", the applicable config file, raw reads, assembled fasta, and metadata will be generated at the specified ", ui.code(ui.strong("submission_dir")), ". ",

    ),
]

seqsender_requirements = [
    ui.h4("Config File"),
    ui.p("SeqSender uses a config file to simplify, keep consistent, and automate the submission file creation/upload. ",
         "Fields required in your config file can change based on your database and organism selection. ",
         "Use the tab: ", ui.strong("Submission Wizard"), " to generate a config file based on your database selection. ",
    ),
    ui.h4("Metadata Specific Fields"),
    ui.p("SeqSender has specific required fields depending on which database you submit to that must be present in your metadata. ",
        "All fields below are required unless specified otherwise and will be used repeatedly for when submitting to multiple databases. "
    ),
    ui.tags.ul(
        ui.tags.li(ui.strong("sequence_name")),
        ui.tags.ul(
            ui.p("Required for when submitting to ", ui.strong("GenBank"), " or ", ui.strong("GISAID"), ". ",
                "This field should contain the sequence header name of the fasta file you are using with SeqSender. ",
            ),
            ui.p("SeqSender will use this to join your metadata with your fasta file. ",
                "This is to create the corresponding fasta file's for submission to GISAID/GenBank with their updated submission name based on the name used for ", ui.code(ui.strong("gb-sample_name")), " or ", ui.code(ui.strong("gs-sample_name")), ". ",
            ),
        ),
        ui.tags.li(ui.strong("organism")),
        ui.tags.ul(
            ui.p("This should be the full name of the organism for the sample. ",
                 "Be sure that you have correctly spelled the full name of your organism to prevent submission issues. ",
            ),
        ),
        ui.tags.li(ui.strong("authors")),
        ui.tags.ul(
            ui.p("This should be a one line, list of names associated to the sample you are submitting.",
                 "The names should be structured in the format of ", ui.code(ui.strong("<first name> <middle name> <last name>")), " with each name separated by a semi-colon. Middle names can be excluded. ",
                 "Your list of names should look something like this: ",
            ),
            ui.code(ui.strong("John Doe; Jane Doe; John Quincy Adams; Franklin Delano Roosevelt")),
            ui.p(ui.strong("Note: "), "SeqSender uses a natural language processing library to sort names into first, middle, and last. ",
                 "While SeqSender does handle international and hyphenated names fairly well, it can still have issues correctly identifying every name's structure. ",
                 "If you find this is an issue for your name, try either removing/initializing your middle name or by removing the hypen from your name if present. ",
            ),
        ),
        ui.tags.li(ui.strong("collection_date")),
        ui.tags.ul(
            ui.p("This should be an ISO-8601 formatted collection date for your sample. ",
                 "The format's accepted are: ", ui.code(ui.strong("YYYY-MM-DD")), ", ", ui.code(ui.strong("YYYY-MM")), ", or ", ui.code(ui.strong("YYYY")), ". ",
                 "Time should not be included as it could create issues during submission. "
            ),
        ),
        ui.tags.li(ui.strong("bioproject")),
        ui.tags.ul(
            ui.p("This field is required when submitting to SRA or BioSample but is recommended to be used for all NCBI database submissions. "
                "Visit the ", ui.a(ui.strong("BioProject Web Portal"), href="https://submit.ncbi.nlm.nih.gov/subs/bioproject/"), " to create a BioProject for your samples. ",
            ),
            ui.p("To add your BioProject to your metadata, all you need is the accession ID ", ui.code(ui.strong("(PRJNAxxxxxx, i.e. PRJNA123456)")), " created once you complete the BioProject."),
        ),
    ),
]

sra_database_submission = [
    ui.h4("Raw Reads Locations"),
    ui.p("The raw reads for your samples can be stored in a variety of different locations and still be used in submission with SeqSender. ",
         "Raw reads can either be stored locally or in the cloud when uploading to SRA, as well as, multiple files can be associated to one sample. ",
         "The metadata column ", ui.code(ui.strong("sra-file_#")), ", can be repeated where # is the numeric value of that file. ",
         "(i.e. If you were uploading two files, you would have the first file path in ", ui.code("sra-file_1"), " and the second file path in ", ui.code("sra-file_2"), ".)"
    ),
    ui.p("If your raw reads are stored in the cloud you will want to include the full URL path to the file for the metadata column. ",
         "To signify to SeqSender that this is a cloud URL, in the metadata column ", ui.code(ui.strong("sra-file_location")), ", the keyword ", ui.code(ui.strong("cloud")), " should be included. ",
    ),
    ui.p("If your raw reads are locally stored locally they can be located in two locations to be accessed by SeqSender. ",
         "If just the file name is included in your metadata column, then your raw reads should be stored in a folder named ", ui.code(ui.strong("raw_reads")), " at the location of your ", ui.code(ui.strong("--submission_dir")), " when using SeqSender. ",
         "Otherwise, the full file path can be included in your submission. (", ui.strong("Note: "), "Be sure if using docker/singularity that SeqSender has access to the location of your raw reads to prevent submission errors. ",
         "For local files, the metadata column ", ui.code(ui.strong("sra-file_location")), " should use the keyword ", ui.code(ui.strong("local")), " to signify that the file is stored locally. ",
    ),
    ui.h4("Raw Reads Processing Error"),
    ui.p("If your raw reads fail to process on SRA and it could be because they are incorrectly formatted for SRA. ",
         "If the SRA team notifies you of this and recommends adding a ", ui.strong("sra loader"), " to correct the problem, that can be easily added to SeqSender. ",
         "Using the metadata column", ui.code(ui.strong("sra-loader")), " add the name of the specified loader by the SRA team to correct this issue. ",
    ),
    ui.h4("Metadata Fields"),
    ui.p("For all required fields for SeqSender to submit data to SRA, refer to the tab: ", ui.strong("Submission Wizard"), ". ",
        "For all other optional fields, refer to NCBI's ", ui.a(ui.strong("SRA portal"), href="https://www.ncbi.nlm.nih.gov/sra/docs/submitportal/"), ", any metadata field can be easily added by using the specified column name with the prefix ", ui.code("sra-"), ". ",
    ),
]

biosample_database_submission = [
    ui.h4("Metadata Fields"),
    ui.p("Required metadata columns change based on the specified ", ui.a(ui.strong("BioSample Package"), href=""), " used in your submission. ",
         "Refer to the tab: ", ui.strong("Submission Wizard"), " to get all the required metadata columns for your chosen package. ",
         "If you want to add a supported column that is not present in the submission wizard, you can easily add it to your metadata using the specified column name with the prefix ", ui.code("bs-"), ". ",
    ),
]

first_submission_body = [
    ui.h2("My First Submission"),
    ui.hr(),
    ui.navset_pill_list(
        ui.nav_panel("Generate Test Data", test_data),
        ui.nav_panel("SeqSender Requirements", seqsender_requirements),
        ui.nav_menu("Database Specific Requirements",
            ui.nav_panel("SRA", sra_database_submission),
            ui.nav_panel("BioSample", biosample_database_submission),
        ),
        id="tab",
    ),
]
