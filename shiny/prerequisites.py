from shiny import App, render, ui

ncbi_prereq = [
    ui.h4("NCBI Submissions:"),
    ui.p("SeqSender utilizes a UI-Less Data Submission Protocol to bulk upload submission files to NCBI databases. SeqSender also uses ", ui.a(ui.strong("table2asn"), href="https://www.ncbi.nlm.nih.gov/genbank/table2asn/"), ", to upload GenBank submissions for organisms other than Influenza and COVID-19, by emailing it to NCBI via a Simple Mail Transfer Protocol (SMTP). ",
        "To begin submitting to NCBI with SeqSender, you must: ",
    ),
    ui.tags.ul(
        ui.tags.li("Have a ", ui.a(ui.strong("NCBI account"), href="https://account.ncbi.nlm.nih.gov/"), ". Sign up or login in to continue. "),
        ui.p("To generate the files for submitting data to NCBI via their ", ui.a(ui.strong("submission portal"), href = "https://submit.ncbi.nlm.nih.gov/"), " or if you're only submitting data via table2asn, no extra steps must be taken. Use the tab: ", ui.strong("Submission Wizard"), " to get the required metadata columns for the database you're uploading to. ",
            "If you're wanting SeqSender to handle uploading your samples for you, then follow the next steps: ",
        ),
        ui.tags.li("Contact NCBI at: ", ui.code(ui.strong("gb-admin@ncbi.nlm.nih.gov")), " to create your institution/group/lab's UI-less submission account. This will create an account for your institution/group/lab's and allow you to specify NCBI users of your institution/group/lab's access to the uploads via NCBI's ", ui.a(ui.strong("web portal"), href="https://submit.ncbi.nlm.nih.gov/subs/"), ". Make sure to answer the following questions below when contacting NCBI: "),
        ui.tags.ul(
            ui.tags.li(ui.strong("MyNCBI account email of the primary submitter")),
            ui.tags.li(ui.strong("center/account abbreviation")),
            ui.tags.li(ui.strong("full center/account name")),
            ui.tags.li(ui.strong("names and email addresses of all additional users")),
            ui.tags.li(ui.strong("postal address of institute (including postal code and country)")),
        ),
        ui.tags.li("Once created, the UI-less submission account credentials and spuid_namespace can be added to your config file to begin submitting your samples to NCBI. Be sure to use the SeqSender submission flag ", ui.code(ui.strong("--test")), " to perform an initial test submission to NCBI, to ensure everything is setup correctly. "),
    )
]
gisaid_prereq = [
    ui.h4("GISAID Submissions"),
    ui.p("SeqSender makes use of GISAID's Command Line Interface tools to bulk upload sample organism's metadata and assembly sequence to the GISAID databases. Currently the following databases are supported for GISAID: EpiFLU, EpiCoV, EpiArbo, EpiRSV. ",
        "To begin submitting data to GISAID using SeqSender, you must: ",
    ),
    ui.tags.ul(
        ui.tags.li(
            ui.p("Have a GISAID account and access to the Epi database you're wanting to upload to. ",
                ui.a(ui.strong("Sign up"), href="https://gisaid.org/register/"), " for a GISAID account if you do not have one and be sure you can access the database on your ", ui.a(ui.strong("submission portal"), href="https://www.epicov.org/epi3/frontend#"), " to ensure you have met the data access agreements for GISAID. ",
            ),
        ),
        ui.p("To generate the files for submitting data to GISAID via their ", ui.a(ui.strong("submission portal"), href = "https://www.epicov.org/epi3/frontend#"), ", no extra steps must be taken. Use the tab: ", ui.strong("Submission Wizard"), " to get the required metadata columns for the database you're uploading to. ",
            "If you're wanting SeqSender to handle uploading your sequences for you, then follow the next steps: ",
        ),
        ui.tags.li(
            ui.p("Request a ", ui.strong("Client-ID"), " from GISAID for the database you're wanting to upload to by emailing them at: ", ui.code(ui.strong("clisupport@gisaid.org")), ". "),
        ),
        ui.tags.li(
            ui.p("Download the CLI tool for your chosen database from GISAID's ", ui.a(ui.strong("submission portal"), href="https://www.epicov.org/epi3/frontend#"), " and place the ", ui.strong("CLI File"), " in one of two locations to be accessed by SeqSender: "),
            ui.tags.ul(
                ui.tags.li("Place it in SeqSender's local directory in an a folder called: ", ui.code(ui.strong("seqsender/gisaid_cli/")), ". i.e.(seqsender/gisaid_cli/epiFLU) or (seqsender/gisaid_cli/epiPOX/epiPOX)"),
                ui.tags.li("Place it in the directory ", ui.code(ui.strong("--submission_dir")), " when running SeqSender in a folder called: ", ui.code(ui.strong("<--submission_dir>/gisaid_cli/")), ". i.e.(my_folder/gisaid_cli/epiCOV) or (my_folder/gisaid_cli/epiFLU/EpiFLU)")
            ),
        ),
        ui.tags.li(
            ui.p("Make a test submission first with the ", ui.strong("test Client-ID"), " provided by GISAID when you requested your client ID. Use the tab: ", ui.strong("Submission Wizard"), " to input your client-ID and login credentials for GISAID to upload your sequences via their CLI tool. "),
        ),
    ),
]
prerequisites_body = [
    ui.h2("Prerequisites"),
    ui.navset_tab(
        ui.nav_panel("NCBI", ncbi_prereq),
        ui.nav_panel("GISAID", gisaid_prereq),
        id="prerequisites_tab",
    ),
]
