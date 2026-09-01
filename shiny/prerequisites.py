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
prerequisites_body = [
    ui.h2("Prerequisites"),
    ui.navset_tab(
        ui.nav_panel("NCBI", ncbi_prereq),
        id="prerequisites_tab",
    ),
]
