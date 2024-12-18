from shiny import App, Inputs, Outputs, Session, render, req, ui, reactive
from htmltools import TagList, div


terminal_css = ""

def software_requirements(version):
    requirements = """<br>
<p><strong>SOFTWARE REQUIREMENTS:</strong></p>
<ul>
<li>Linux (64-bit) or Mac OS X (64-bit)</li>
<li>Git version 2.25.1 or later</li>
<li>Standard utilities: curl, tar, unzip</li>
"""
    if version == "Local":
        return ui.HTML(requirements + """</ul>""")
    elif version == "Docker":
        return ui.HTML(
            requirements
            + """
<li>Docker version 20.10.14 or later</li>
<li>[Optional] Docker Compose version 2.21 or later</li>
</ul>
"""
        )
    elif version == "Singularity":
        return ui.HTML(
            requirements
            + """
<li>Singularity version 3.8.7 or later</li>
</ul>
"""
        )
    else:
        return ""

def seqsender_help_output_msg(version):
    message = """
usage: seqsender.py [-h] {prep,submit,check_submission_status,template,update_biosample,version} ...</p>
<p>Automate the process of batch uploading consensus sequences and metadata to databases of your choices</p>
<p>positional arguments:<br>
&nbsp;|&nbsp;{prep,submit,check_submission_status,template,update_biosample,version}</p>
<p>optional arguments:<br>
&nbsp;|&nbsp;-h, --help&nbsp;|&nbsp;&nbsp;show this help message and exit</p>
"""
    if version == "Local":
        command = (
            """<p>(seqsender) [USER@server seqsender]$ python seqsender.py --help<br>"""
        )
    elif version == "Docker":
        command = """<p>docker exec -it seqsender bash seqsender-kickoff --help<br>"""
    elif version == "Singularity":
        command = """<p>singularity exec ~/singularity/seqsender.sif seqsender-kickoff --help<br>"""
    else:
        command = """<p>"""
    return ui.card(ui.HTML(command + message), style=terminal_css)


databases_parameter = [
    ui.strong("Databases"),
    ui.tags.ul(
        ui.tags.li(
            ui.strong(
                ui.code("--biosample"),
                " | ",
                ui.code("-b"),
            ),
            ui.tags.ul(
                "  Generate files required for submission to ",
                ui.strong("NCBI Database: BioSample"),
                ".",
                inline=True,
            ),
        ),
        ui.tags.li(
            ui.strong(
                ui.code("--sra"),
                " | ",
                ui.code("-s"),
            ),
            ui.tags.ul(
                "  Generate files required for submission to ",
                ui.strong("NCBI Database: SRA"),
                ".",
                inline=True,
            ),
        ),
        ui.tags.li(
            ui.strong(
                ui.code("--genbank"),
                " | ",
                ui.code("-n"),
            ),
            ui.tags.ul(
                "  Generate files required for submission to ",
                ui.strong("NCBI Database: GenBank"),
                ".",
                inline=True,
            ),
        ),
        ui.tags.li(
            ui.strong(
                ui.code("--gisaid"),
                " | ",
                ui.code("-g"),
            ),
            ui.tags.ul(
                "  Generate files required for submission to ",
                ui.strong("GISAID"),
                ". Not available as an option when using the organism parameter ",
                ui.code("OTHER"),
                ".",
                inline=True,
            ),
        ),
    ),
]

organisms_parameter = [
    ui.strong("Organism"),
    ui.div(
        ui.strong(ui.code("--organism "), "{'FLU', 'COV', 'POX', 'ARBO', 'RSV', 'OTHER'}"),
        inline=True,
    ),
    ui.p(
        "Organism selection adjusts seqsender submission process for organism specific criteria."
    ),
    ui.tags.ul(
        ui.tags.li(
            ui.strong(
                ui.code("FLU"),
            ),
            ui.tags.ul(
                " For ",
                ui.strong("Influenza A virus"),
                ", enables FTP GenBank submission and GISAID CLI submission.",
            ),
        ),
        ui.tags.li(
            ui.strong(
                ui.code("COV"),
            ),
            ui.tags.ul(
                " For ",
                ui.strong("SARS-CoV-2"),
                ", enables FTP GenBank submission and GISAID CLI submission.",
            ),
        ),
        ui.tags.li(
            ui.strong(
                ui.code("POX"),
            ),
            ui.tags.ul(
                " For ",
                ui.strong("Mpox"),
                " (monkeypox), enables GISAID CLI submission.",
            ),
        ),
        ui.tags.li(
            ui.strong(
                ui.code("ARBO"),
            ),
            ui.tags.ul(
                " For ", ui.strong("Arbovirus"), ", enables GISAID CLI submission."
            ),
        ),
        ui.tags.li(
            ui.strong(
                ui.code("RSV"),
            ),
            ui.tags.ul(
                " For ", ui.strong("Respiratory syncytial virus"), ", enables GISAID CLI submission."
            ),
        ),
        ui.tags.li(
            ui.strong(
                ui.code("OTHER"),
            ),
            ui.tags.ul(
                " For any organism that does not have special submission options."
            ),
        ),
    ),
]

submission_directory_parameter = ui.div(
    ui.strong(ui.code("--submission_dir")),
    ui.tags.ul("Directory where all files required for submission are stored."),
)

submission_name_parameter = ui.div(
    ui.strong(ui.code("--submission_name")),
    ui.tags.ul("Unique name to use for submission. Used as a unique identifier for submission to databases and to name the files created for submission."),
)

config_parameter = ui.div(
    ui.strong(ui.code("--config_file")),
    ui.tags.ul("Full path to config file if not stored in ", ui.code("--submission_dir"), " location."),
)

metadata_parameter = ui.div(
    ui.strong(ui.code("--metadata_file")),
    ui.tags.ul("Full path to metadata file if not stored in ", ui.code("--submission_dir"), " location."),
)

fasta_parameter = ui.div(
    ui.strong(ui.code("--fasta_file")),
    ui.tags.ul("Full path to fasta file if not stored in ", ui.code("--submission_dir"), " location."),
)

table2asn_parameter = ui.div(
    ui.strong(ui.code("--table2asn")),
    ui.tags.ul("Flag to enable a \"table2asn\" submission for GenBank. Automatically enabled for organisms other than \"FLU\" and \"COV\"."),
)

gff_parameter = ui.div(
    ui.strong(ui.code("--gff_file")),
    ui.tags.ul("Full path to gff file if not stored in ", ui.code("--submission_dir"), " location."),
)

test_parameter = ui.div(
    ui.strong(ui.code("--test")),
    ui.tags.ul("Flag to enable a \"test\" submission to each database selected."),
)

validation_parameter = ui.div(
    ui.strong(ui.code("--skip_validation")),
    ui.tags.ul("Flag to skip pandera validation for ", ui.code("--metadata_file"), ". Warning, skipping validation can cause unexpected errors when missing required data."),
)

def seqsender_submit_help_output_msg(version):
    message = """
usage: seqsender.py submit [-h] [--biosample] [--sra] [--genbank] [--gisaid] --organism {FLU,COV,POX,ARBO,RSV,OTHER}
--submission_name SUBMISSION_NAME --submission_dir SUBMISSION_DIR --config_file CONFIG_FILE --metadata_file METADATA_FILE
--fasta_file FASTA_FILE [--table2asn] [--gff_file GFF_FILE] [--test]</p>
<p>Create submission files and then batch uploading them to databases of choices.</p>
<p>optional arguments:<br>
&nbsp;-h, --help &nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;
show this help message and exit<br>
&nbsp;--biosample, -b &nbsp;|&nbsp;&nbsp;|&nbsp;
Submit to BioSample database. (default: )<br>
&nbsp;--sra, -s &nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;
Submit to SRA database. (default: )<br>
&nbsp;--genbank, -n &nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;
Submit to Genbank database. (default: )<br>
&nbsp;--gisaid, -g &nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;
Submit to GISAID database. (default: )<br>
&nbsp;--organism {FLU,COV,POX,ARBO,RSV,OTHER} <br>&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;
Type of organism data (default: FLU)<br>
&nbsp;--submission_name SUBMISSION_NAME <br>&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;
Name of the submission (default: None)<br>
&nbsp;--submission_dir SUBMISSION_DIR  <br>&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;
Directory to where all required files (such as metadata, fasta, etc.) are stored (default: None)<br>
&nbsp;--config_file CONFIG_FILE  <br>&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;
Config file stored in submission directory (default: None)<br>
&nbsp;--metadata_file METADATA_FILE  <br>&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;
Metadata file stored in submission directory (default: None)<br>
&nbsp;--fasta_file FASTA_FILE  <br>&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;
Fasta file stored in submission directory (default: None)<br>
&nbsp;--table2asn &nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;
Whether to prepare a Table2asn submission. (default: False)<br>
&nbsp;--gff_file GFF_FILE &nbsp;
An annotation file to add to a Table2asn submission (default: None)<br>
&nbsp;--test &nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;
Whether to perform a test submission. (default: False)<br>
"""
    if version == "Local":
        command = """<p>(seqsender) [USER@server seqsender]$ python seqsender.py submit --help<br>"""
    elif version == "Docker":
        command = (
            """<p>docker exec -it seqsender bash seqsender-kickoff submit --help<br>"""
        )
    elif version == "Singularity":
        command = """<p>singularity exec ~/singularity/seqsender.sif seqsender-kickoff submit --help<br>"""
    else:
        command = """"""
    return ui.card(ui.HTML(command + message), style=terminal_css)

# Returns parameters in format for code example and descriptive list of parameter
def get_parameters(command, output):
    command_string = ""
    parameters = []
    # Databases parameter
    if command in ["prep", "submit", "test_data"]:
        command_string += " [&lt;databases&gt;]"
        parameters.append(databases_parameter)
    # Organisms parameter
    if command in ["prep", "submit", "test_data"]:
        command_string += " [--organism &lt;organism&gt;]"
        parameters.append(organisms_parameter)
    # Submission_dir parameter
    if command in ["prep", "submit", "test_data", "submission_status"]:
        command_string += " [--submission_dir &lt;directory_path&gt;]"
        parameters.append(submission_directory_parameter)
    # Add submission_name parameter
    if command in ["prep", "submit", "submission_status"]:
        command_string += " [--submission_name &lt;name&gt;]"
        parameters.append(submission_name_parameter)
    # Add config_file parameter
    if command in ["prep", "submit"]:
        command_string += " [--config_file &lt;config_file_path&gt;]"
        parameters.append(config_parameter)
    # Add metadata_file parameter
    if command in ["prep", "submit"]:
        command_string += " [--metadata_file &lt;metadata_file_path&gt;]"
        parameters.append(metadata_parameter)
    # Add fasta_file parameter
    if command in ["prep", "submit"]:
        command_string += " [--fasta_file &lt;fasta_file_path&gt;]"
        parameters.append(fasta_parameter)
    # Add table2asn parameter
    if command in ["prep", "submit"]:
        command_string += " [--table2asn]"
        parameters.append(table2asn_parameter)
    # Add gff_file parameter
    if command in ["prep", "submit"]:
        command_string += " [--gff_file &lt;gff_file_path&gt;]"
        parameters.append(gff_parameter)
    # Add test parameter
    if command in ["submit"]:
        command_string += " [--test]"
        parameters.append(test_parameter)
    # Add skip pandera validation parameter
    if command in ["prep", "submit"]:
        command_string += " [--skip_validation]"
        parameters.append(validation_parameter)
    if output == "string":
        return command_string
    elif output == "description":
        return parameters

def seqsender_example_command(platform, command, flags):
    # Platform specific layout
    message = "<pre><code>"
    if platform == "local":
        message += "python seqsender.py "
    elif platform == "docker":
        message += "docker exec -it seqsender bash seqsender-kickoff "
    elif platform == "singularity":
        message += "singularity exec ~/singularity/seqsender.sif seqsender-kickoff "
    else:
        message += ""
    # Seqsender command
    message += command
    ### Add parameters
    message += get_parameters(command, "string") + "</pre></code>"
    return ui.HTML(message)

def command_submenu_accordion(command):
    submenu = []
    open_tabs = []
    parameters = get_parameters(command, "description")
    if len(parameters) > 0:
        submenu.append(ui.accordion_panel("Parameters", parameters))
        open_tabs.append("Parameters")
    notes = [] # notes function
    if len(notes) > 0:
        submenu.append(ui.accordion_panel("Notes", notes))
        open_tabs.append("Notes")
    examples = [] # examples function
    if len(examples) > 0:
        submenu.append(ui.accordion_panel("Examples", examples))
    if len(submenu) > 0:
        return ui.accordion(*submenu, open = open_tabs)
    else:
        return ""

def command_accordion_panel(command, description):
    return ui.accordion_panel(
        command.capitalize(),
        ui.p(
            "The ",
            ui.code(command),
            description,
        ),
        seqsender_example_command("local", command, ""),
        command_submenu_accordion(command),)

def file_output_column_info(column_name, description, controlled_fields):
    list_item = ui.tags.li(
        ui.strong(
            ui.code(column_name),
        ),
        ui.tags.ul(
            ui.p(description),
            file_output_controlled_field_info(controlled_fields),
        )
    )
    return list_item

def file_output_controlled_field_info(controlled_fields):
    keywords = []
    if controlled_fields:
        for key, description in controlled_fields:
            keywords.append(ui.tags.li(ui.p(ui.strong(key, ": "), description), style="margin-top:10px;"))
    return keywords

def create_help_tooltip(id, description, position):
    return ui.tooltip(ui.p("❓"),description,id=id,style="display:inline-block;float:" + position + ";",)
