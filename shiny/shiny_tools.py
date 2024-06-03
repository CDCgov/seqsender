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


def seqsender_submit_help_output_msg(version):
    message = """
usage: seqsender.py submit [-h] [--biosample] [--sra] [--genbank] [--gisaid] --organism {FLU,COV,POX,ARBO,OTHER}
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
&nbsp;--organism {FLU,COV,POX,ARBO,OTHER} <br>&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;
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
    # Add databases parameter
    if command in ["prep", "submit", "test_data"]:
        message += " [&lt;databases&gt;]"
    # Add organisms parameter
    if command in ["prep", "submit", "test_data"]:
        message += " [--organism &lt;organism&gt;]"
    # Add submission_dir parameter
    if command in ["prep", "submit", "test_data"]:
        message += " [--submission_dir &lt;directory_path&gt;]"
    # Add submission_name parameter
    if command in ["prep", "submit"]:
        message += " [--submission_name &lt;name&gt;]"
    # Add config_file parameter
    if command in ["prep", "submit"]:
        message += " [--config_file &lt;config_file_path&gt;]"
    # Add metadata_file parameter
    if command in ["prep", "submit"]:
        message += " [--metadata_file &lt;metadata_file_path&gt;]"
    # Add fasta_file parameter
    if command in ["prep", "submit"]:
        message += " [--fasta_file &lt;fasta_file_path&gt;]"
    # Add gff_file parameter
    if command in ["prep", "submit"]:
        message += " [--gff_file &lt;gff_file_path&gt;]"
    # Add test parameter
    if command in ["prep", "submit"]:
        message += " [--test]"
    message += "</pre></code>"
    return ui.HTML(message)
