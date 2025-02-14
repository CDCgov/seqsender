from shiny import App, Inputs, Outputs, Session, render, req, ui, reactive
from htmltools import TagList, div

###################### CSS #######################
import pandas as pd
from shinyswatch import theme
import shiny_tools
from index import index_body
from setup import setup_body
from about import about_body
from prerequisites import prerequisites_body
from my_first_submission import first_submission_body
import pathlib
import yaml
from yaml import SafeDumper

terminal_css = ""
yaml_css = "background-color: #F0F0F0;white-space: nowrap; font-size: 20px ;margin-top:-15px;font-family: Consolas,Monaco,Lucida Console,Liberation Mono,DejaVu Sans Mono,Bitstream Vera Sans Mono,Courier New, monospace;-webkit-user-select: none; -ms-user-select: none; user-select: none;"
############ HEADER #######################
header = (
    ui.card_header(
        ui.HTML(
            """<p><strong>Beta Version</strong>: 1.3.2. This pipeline is currently in Beta testing, and issues could appear during submission. Please use it at your own risk. Feedback and suggestions are welcome!</p>"""
        )
    ),
)
############## FOOTER #######################
footer = (
    ui.card_footer(
        ui.HTML(
            """<p><strong>General Disclaimer</strong>: This repository was created for use by CDC programs to collaborate on public health related projects in support of the <a href="https://www.cdc.gov/about/organization/mission.htm">CDC mission</a>. GitHub is not hosted by the CDC, but is a third party website used by CDC and its partners to share information and collaborate on software. CDC use of GitHub does not imply an endorsement of any one particular service, product, or enterprise.</p>"""
        ),
    ),
)
############### FUNCTIONS #########################
ui.input_checkbox_group(
    "select_databases",
    "Select Database(s):",
    {
        "GenBank": ui.span("GenBank"),
        "SRA": ui.span("SRA"),
        "BioSample": ui.span("BioSample"),
        "blue": ui.span("GISAID"),
    },
),

####################### INSTALLATION PAGE ###########################
#### LOCAL INSTALLATION
local_installation_content = [
    shiny_tools.software_requirements("Local"),
    ui.HTML(
        """
<h2>Micromamba Installation</h2>
<p>Here we recommend using <strong>micromamba</strong> to set up a virtual environment to run <code>seqsender</code>. <strong>Micromamba</strong> is a tiny, statically linked C++ reimplementation of mamba which is an alternative to conda. The tool works as a standalone package manager that supports a subset of all mamba or conda commands, but it also has its own separate command line interfaces. For more information, visit <a href="https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html">micromamba documentation</a>.</p>
<p>To manually install, download and unzip the executable from the official <strong>conda-forge</strong> package to your <code>$HOME</code> directory using <code>tar</code>.</p>
<pre><code>cd $HOME</code></pre>
<strong>LINUX</strong>
<pre><code># Linux Intel (x86_64):
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
# Linux ARM64:
curl -Ls https://micro.mamba.pm/api/micromamba/linux-aarch64/latest | tar -xvj bin/micromamba
# Linux Power:
curl -Ls https://micro.mamba.pm/api/micromamba/linux-ppc64le/latest | tar -xvj bin/micromamba</code></pre>
<strong>macOS</strong>
<pre><code># macOS Intel (x86_64):
curl -Ls https://micro.mamba.pm/api/micromamba/osx-64/latest | tar -xvj bin/micromamba
# macOS Silicon/M1 (ARM64):
curl -Ls https://micro.mamba.pm/api/micromamba/osx-arm64/latest | tar -xvj bin/micromamba</code></pre>
<p>After the extraction is completed, you can find the executable at <code>$HOME/bin/micromamba</code></p>
<ul>
<li>To quickly use <code>micromamba</code>, you can simply run</li>
</ul>
<pre><code>export MAMBA_ROOT_PREFIX=&quot;$HOME/micromamba&quot;
eval &quot;$($HOME/bin/micromamba shell hook -s posix)&quot;</code></pre>
<ul>
<li>To persist using <code>micromamba</code>, you can append the following script to your <code>.bashrc</code> (or <code>.zshrc</code>)</li>
</ul>
<pre><code># &gt;&gt;&gt; mamba initialize &gt;&gt;&gt;
export MAMBA_EXE=&quot;$HOME/bin/micromamba&quot;;
export MAMBA_ROOT_PREFIX=&quot;$HOME/micromamba&quot;;
__mamba_setup=&quot;$(&quot;$MAMBA_EXE&quot; shell hook --shell bash --root-prefix &quot;$MAMBA_ROOT_PREFIX&quot; 2&gt; /dev/null)&quot;
if [ $? -eq 0 ]; then
    eval &quot;$__mamba_setup&quot;
else
    alias micromamba=&quot;$MAMBA_EXE&quot;  # Fallback on help from mamba activate
fi
unset __mamba_setup
# &lt;&lt;&lt; mamba initialize &lt;&lt;&lt;</code></pre>
<ul>
<li>To check the current version of <code>micromamba</code></li>
</ul>
<pre><code>micromamba --version
1.5.6</code></pre>
<h2>Set up a <code>micromamba</code> environment</h2>
<ol style="list-style-type: decimal">
<li>Clone this repository to your <code>$HOME</code> directory</li>
</ol>
<pre><code>cd $HOME
git clone https://github.com/CDCgov/seqsender.git</code></pre>
<ol start="2" style="list-style-type: decimal">
<li><code>CD</code> to <strong>seqsender</strong> folder where the <code>env.yaml</code> file is stored. Let’s create a virtual environment named <strong>mamba</strong> that contains all dependencies needed to run <code>seqsender</code> from the source file.</li>
</ol>
<pre><code>cd seqsender
micromamba create --name seqsender --file env.yaml</code></pre>
<ol start="3" style="list-style-type: decimal">
<li>Activate the named environment – <strong>seqsender</strong></li>
</ol>
<pre><code>micromamba activate seqsender</code></pre>
<h2>Run <code>seqsender</code> within the mamba environment</h2>
<pre><code>python seqsender.py --help</code></pre>"""
    ),
    shiny_tools.seqsender_help_output_msg("Local"),
    ui.HTML(
        """
<p>To see the arguments required for each command, for example, the <code>submit</code> command, run</p>
<pre><code>python seqsender.py submit --help</code></pre>
"""
    ),
    shiny_tools.seqsender_submit_help_output_msg("Local"),
]
#### DOCKER INSTALLATION
docker_installation_content = [
    shiny_tools.software_requirements("Docker"),
    ui.HTML(
        """
<h2>Clone <code>seqsender</code> repo to your $HOME directory and navigate to the repo</h2>
<pre><code>cd $HOME
git clone https://github.com/CDCgov/seqsender.git
cd seqsender</code></pre>"""
    ),
    ui.navset_card_pill(
        ui.nav_panel(
            "Docker-Build",
            [
                ui.HTML(
                    """
<h2>In the directory where the <code>Dockerfile</code> is stored, build docker image</h2>
<pre><code>docker build -t seqsender:latest .</code></pre>
<p><strong>-t</strong>: add a tag to an image, e.g. <em>seqsender:1.0.0</em> or <em>seqsender:latest</em></p>
<br>
<h2>After the build is complete, you can check if the image is created successfully.</h2>
<pre><code>docker images</code></pre>"""
                ),
                ui.card(
                    ui.HTML(
                        """
docker images<br>
REPOSITORY&nbsp;TAG&nbsp;|&nbsp;&nbsp;IMAGE ID&nbsp;|&nbsp;&nbsp;CREATED&nbsp;|&nbsp;SIZE<br>
seqsender&nbsp;latest&nbsp;d9e2578d2211&nbsp;2 weeks ago&nbsp;581GB
"""
                    ),
                    style=terminal_css,
                ),
                ui.HTML(
                    """
<h2>Run <code>seqsender</code> container</h2>
<pre><code>docker run \
-v $HOME:/data \
-t -d seqsender:latest \
--name seqsender</code></pre>
<p><strong><code>-t</code></strong>: allocate a pseudo-tty <br> <strong><code>-d</code></strong>: run the container in detached mode <br> <strong><code>-v</code></strong>: mount data files from host directory to container directory <strong>[host_div]:[container_dir]</strong>. By exposing the host directory to docker container, docker will be able to access data files within that mounted directory and use it to fire up the <code>seqsender</code>workflows. <strong>NOTE:</strong> Here we are mounting the local <code>$HOME</code> directory to <code>/data</code> directory inside the container.<br> <strong><code>--name</code></strong>: give an identity to the container <br></p>
<p>For more information about the Docker syntax, see <a href="https://docs.docker.com/engine/reference/run/" target="_blank">Docker run reference</a></p>"""
                ),
            ],
        ),
        ui.nav_panel("Docker-Compose",
            ui.h2("Update docker-compose.yaml for your local storage"),
            ui.p(ui.strong("Note: "), ui.code("source"), " is the storage location of your local machine. This location will be mapped to ", ui.code("/data"), " directory inside of the container. Here we are mounting the local variable ", ui.code("$HOME"), " directory to ", ui.code("/data"), " inside of the container. Change the value ", ui.code("$HOME"), " to reflect your local storage system."),
            ui.card(ui.HTML(
                """
<code>version: "3.9"<br>
x-data-volumes:<br>
&nbsp;&nbsp;&data-volume<br>
&nbsp;&nbsp;type: bind<br>
&nbsp;&nbsp;source: $HOME<br>
&nbsp;&nbsp;target: /data<br>
<br>
services:<br>
&nbsp;&nbsp;seqsender:<br>
&nbsp;&nbsp;&nbsp;&nbsp;container_name: seqsender<br>
&nbsp;&nbsp;&nbsp;&nbsp;image: seqsender:latest<br>
&nbsp;&nbsp;&nbsp;&nbsp;build:<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;context: .<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;dockerfile: Dockerfile<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;args:<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- micromamba_version=1.5.3<br>
&nbsp;&nbsp;&nbsp;&nbsp;restart: always<br>
&nbsp;&nbsp;&nbsp;&nbsp;volumes:<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- *data-volume<br>
&nbsp;&nbsp;&nbsp;&nbsp;command: tail -f /dev/null<br>
</code>"""
                ),
            style=terminal_css,
            ),
            ui.h2("Start ", ui.code("seqsender"), " container"),
            ui.card(
                ui.p(ui.code("docker-compose up -d")),
            style=terminal_css,
            ),
            ui.p("Running docker-compose with the flag ", ui.code("-d"), " runs the container in detached mode, allowing the container to continue to run in the background. For more information refer to ", ui.a("docker-compose documentation", href="https://docs.docker.com/reference/cli/docker/compose/up/"), "."),
        ),
        id="docker_options",
    ),
    ui.HTML(
        """
<h2>Check if the container is created successfully</h2>
<pre><code>docker container ps</code></pre>
"""
    ),
    ui.card(
        ui.HTML(
            """
docker container ps<br>
CONTAINER ID&nbsp;|&nbsp;&nbsp;|&nbsp;IMAGE&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;COMMAND&nbsp;|&nbsp;&nbsp;|&nbsp;CREATED&nbsp;|&nbsp;&nbsp;STATUS&nbsp;|&nbsp;PORTS&nbsp;|&nbsp;NAMES<br>
b37b6b19c4e8&nbsp;seqsender:latest&nbsp;&quot;/bin/bash&quot;&nbsp;5 hours ago&nbsp;Up 5 hours&nbsp;|&nbsp;&nbsp;|&nbsp;seqsender
"""
        ),
        style=terminal_css,
    ),
    ui.HTML(
        """
<h2>See a list of commands in <code>seqsender</code> container</h2>
<pre><code>docker exec -it seqsender bash seqsender-kickoff --help</code></pre>
<p><strong><code>-t</code></strong>: allocate a pseudo-tty <br> <strong><code>-i</code></strong>: keep STDIN open even if not attached <br> <strong><code>-h</code></strong>, <strong><code>--help</code></strong>: show help messages and exit</p>
"""
    ),
    shiny_tools.seqsender_help_output_msg("Docker"),
]
#### SINGULARITY INSTALLATION
singularity_installation_content = [
    shiny_tools.software_requirements("Singularity"),
    ui.HTML(
        """
<h2>Convert seqsender Docker image into a Singularity image</h2>
<p>There is a <code>seqsender</code> Docker image already built and stored on our DockerHub registry: <strong>cdcgov/seqsender-dev:latest</strong>. You can directly pull the Docker Image down from the registry, convert it into a Singularity image, and store it in a destination of your choice.</p>
<pre><code>singularity build ~/singularity/seqsender.sif docker://cdcgov/seqsender-dev:latest</code></pre>
<h2>After the Singularity image is built successfully, we can go ahead and use it to run <code>seqsender</code></h2>
<p>Here is the command that shows the help messages of <code>seqsender</code></p>
<pre><code>mkdir ~/singularity
singularity exec ~/singularity/seqsender.sif seqsender-kickoff --help</code></pre>
"""
    ),
    shiny_tools.seqsender_help_output_msg("Singularity"),
    ui.HTML(
        """
<p>To see the arguments required for each command, for example, the <code>submit</code> command, run</p>
<pre><code>singularity exec ~/singularity/seqsender.sif seqsender-kickoff submit --help</code></pre>
"""
    ),
    shiny_tools.seqsender_submit_help_output_msg("Singularity"),
]
installation_body = [
    ui.h2("Installation"),
    ui.navset_tab(
        ui.nav_panel("Local", local_installation_content),
        ui.nav_panel("Docker", docker_installation_content),
        ui.nav_panel("Singularity",singularity_installation_content),
        id="installation_tab",
    ),
]


####################### TEST SUBMISSIONS PAGE ##############################
user_data_content = (
    [
        ui.p("Content Panel A"),
    ],
)

example_data_content = (
    [
        ui.br(),
        ui.p("Access the example data by:"),
        ui.p(
            "First, use the ",
            ui.code("prep"),
            " command to generate all the necessary files for submission. This step is not required for submitting to databases, as this command is ran when running the ",
            ui.code("submit"),
            " command. However, it is helpful to test/debug your submission without submitting to a database and allows you to visualize all of the files generated for submission.",
        ),
        shiny_tools.seqsender_example_command("local", "prep", ""),
        ui.p(
            "After you've validated that your data is being generated correctly using command ",
            ui.code("prep"),
            ", you can now make a test submission to the database. Database submissions are made using the ",
            ui.code("submit"),
            " command, be sure to always include the flag ",
            ui.code("--test"),
            " when submitting test submissions.",
        ),
        shiny_tools.seqsender_example_command("local", "submit", "--test"),
        ui.p(
            "Now that your submission has been made you can check its submission progress using the ",
            ui.code("check_submission_status"),
            " command.",
        ),
    ],
)

testing_body = [
    ui.h2("Test Submissions"),
    ui.p(
        "Test submissions can be made with either the provided test files or your own data."
    ),
    ui.navset_tab(
        ui.nav_panel("I have my own data", user_data_content),
        ui.nav_panel("I don't have my own data", example_data_content),
        id="testing_tab",
        selected="I don't have my own data",
    ),
]

####################### OUTPUT FILE PAGE ###################

output_body = [
    ui.h2("SeqSender Output Files"),
    ui.p("Output files described below are split into two categories: ", ui.strong("SeqSender output"), " or ", ui.strong("database specific"), " which lists all the files that can appear in your database submission directory (", ui.code("<--submission_dir>/<--submission_name>/submission_files/<Database Name>"), "). Files generated will also change for each directory based on the SeqSender submission process for the database, metadata provided, and database submission response files. "),
    ui.h4("SeqSender script output:"),
    ui.hr(),
    ui.h6("submission_log.csv"),
    ui.tags.ul(
        ui.p("The ", ui.strong("submission_log.csv"), " is used by Seqender to track the submission status of each of your submissions made. SeqSender updates each ", ui.strong("Submission_Name"), " that does not have a complete ", ui.strong("Submission_Status"), " by loading the relevant ", ui.strong("Submission_Directory"), " and ", ui.strong("Config_File"), "."),
        ui.tags.ul(
            shiny_tools.file_output_column_info(column_name="Submission_Name",
                description=("Unique ", ui.code("--submission_name"), " used when making your submission via the ", ui.code("submit"), " command. Used by SeqSender to locate the currect submission directory and to name your batch submission during the upload process."),
                controlled_fields=None
            ),
            shiny_tools.file_output_column_info(column_name="Organism",
                description=("Submission organism option ", ui.code("--organism"), " when making your submission which can enable certain additional submission options."),
                controlled_fields=[((ui.code("FLU"), "|", ui.code("COV")), ("For ", ui.strong("Influenza Virus A"), " or ", ui.strong("Severe Acute Respiratory Syndrome Coronavirus 2"), ", it enables GISAID and GenBank (via FTP) as submission options.")),
                    ((ui.code("POX"), "|", ui.code("ARBO"), "|", ui.code("RSV")), ("For ", ui.strong("Mpox"), ", ", ui.strong("Arbovirus"), ", or ", ui.strong("Respiratory syncytial virus"), ", it enables GISAID as a submission option.")),
                    ((ui.code("OTHER")), ("For any organism without additional submission options. It provides access to the default available databases: BioSample, SRA, and GenBank (table2asn via email)."))]
            ),
            shiny_tools.file_output_column_info(column_name="Database",
                description="Database being submitted to.",
                controlled_fields=[((ui.code("BIOSAMPLE"), "|", ui.code("SRA"), "|", ui.code("GENBANK-TBL2ASN"), "|", ui.code("GENBANK-FTP"), "|", ui.code("GISAID")), ("Specifies the database being submitted to. For GenBank, it also includes the submission method for table2asn as \"-TBL2ASN\" or FTP as \"-FTP\"."))]
            ),
            shiny_tools.file_output_column_info(column_name="Submission_Type",
                description="Whether the submission you're making is live or a test.",
                controlled_fields=[((ui.code("TEST"), "|", ui.code("PRODUCTION")), ("Test or live (production) submission."))]
            ),
            shiny_tools.file_output_column_info(column_name="Submission_Date",
                description=("Date submission was started with the SeqSender ", ui.code("submit"), " command."),
                controlled_fields=[((ui.code("YYYY-MM-DD")), ("ISO-8601 standard date format. (i.e. 2024-01-01, 2024-12-31)"))]
            ),
            shiny_tools.file_output_column_info(column_name="Submission_Status",
                description=("Current status of the submission to the database specified. \"GISAID\" only uses ", ui.strong("WAITING"), ", ", ui.strong("PROCESSING"), ", ", ui.strong("PROCESSED"), ", and ", ui.strong("ERROR"), ". \"GENBANK-TBL2ASN\" uses ", ui.strong("EMAILED"), " instead of ", ui.strong("PROCESSED"), " to designate the submission is complete."),
                controlled_fields=[((ui.code("SUBMITTED")), ("Submission has been uploaded to NCBI database.")),
                    ((ui.code("CREATED")), ("NCBI is currently loading the submission files.")),
                    ((ui.code("QUEUED")), ("NCBI has queued your submission for processing.")),
                    ((ui.code("PROCESSING")), ("Submission is currently processing for NCBI database. For GISAID this means that when SeqSender attempted to upload your samples, it was unable to completely upload them all. When the ", ui.code("submission_status"), " command is ran again it will attempt to submit the rest of the samples.")),
                    ((ui.code("FAILED")), ("Submission failed processing for NCBI database.")),
                    ((ui.code("PROCESSED")), ("Submission has been successfully uploaded to NCBI or GISAID database.")),
                    ((ui.code("ERROR")), ("SeqSender failed to process report for NCBI database. For GISAID, SeqSender is failing to upload or samples are unable to be submitted to GISAID.")),
                    ((ui.code("WAITING")), ("Submission is waiting on other submissions to complete processing, to correctly link information.")),
                    ((ui.code("DELETED")), ("NCBI has deleted your submission. This could be because your submission remained errored for too long without resolution or you requested to have the submission removed.")),
                    ((ui.code("RETRIED")), ("NCBI has attempted retrying processing of your submission.")),
                    ((ui.code("EMAILED")), ("SeqSender has successfully emailed your table2asn submission to upload it to GenBank."))]
            ),
            shiny_tools.file_output_column_info(column_name="Submission_ID",
                description=("NCBI submission ID tied to the ", ui.code("Submission_Name"), "."),
                controlled_fields=[((ui.code("PENDING"), "|", ui.code("SUB#"), "|", ui.code("SUBMITTED")), ("Until submission ID is generated it will use ", ui.strong("PENDING"), ". GISAID uses ", ui.strong("SUBMITTED"), " to designate a submission is complete in lieu of a submission ID."))]
            ),
            shiny_tools.file_output_column_info(column_name="Submission_Directory",
                description=("Full file path for the specified ", ui.strong("Submission_Name"), ", ", ui.strong("Database"), ", submission directory."),
                controlled_fields=None
            ),
            shiny_tools.file_output_column_info(column_name="Config_File",
                description=("Full file path for the ", ui.code("--config_file"), " input provided during the ", ui.code("submit"), " command."),
                controlled_fields=None
            ),
            shiny_tools.file_output_column_info(column_name="Update_Date",
                description=("Date submission was last updated with SeqSender via the ", ui.code("submission_status"), " command."),
                controlled_fields=[((ui.code("YYYY-MM-DD")), ("ISO-8601 standard date format. (i.e. 2024-01-01, 2024-12-31)"))]
            ),
        ),
    ),
    ui.h6("submission_status_report.csv"),
    ui.tags.ul(
        ui.p("The ", ui.strong("submission_status_report.csv"), " is used by SeqSender to record/link the samples submitted to each database. When SeqSender is linking accessions between NCBI databases, it will use this file to retrieve the recorded accessions and update the relevant submisison files."),
        ui.tags.ul(
            shiny_tools.file_output_column_info(column_name="bs-sample_name",
                description=("Sample name used for submission to BioSample."),
                controlled_fields=None,
            ),
            shiny_tools.file_output_column_info(column_name="biosample_status",
                description=("Status of the sample submitted to BioSample."),
                controlled_fields=None,
            ),
            shiny_tools.file_output_column_info(column_name="biosample_accession",
                description=("Accession assigned to sample by BioSample."),
                controlled_fields=None,
            ),
            shiny_tools.file_output_column_info(column_name="biosample_message",
                description=("Message related to sample from BioSample."),
                controlled_fields=None,
            ),
            shiny_tools.file_output_column_info(column_name="sra-sample_name",
                description=("Sample name used for submission to SRA."),
                controlled_fields=None,
            ),
            shiny_tools.file_output_column_info(column_name="sra_status",
                description=("Status of the sample submitted to SRA."),
                controlled_fields=None,
            ),
            shiny_tools.file_output_column_info(column_name="sra_accession",
                description=("Accession assigned to sample by SRA."),
                controlled_fields=None,
            ),
            shiny_tools.file_output_column_info(column_name="sra_message",
                description=("Message related to sample from SRA."),
                controlled_fields=None,
            ),
            shiny_tools.file_output_column_info(column_name="gb-sample_name",
                description=("Sample name used for submission to GenBank."),
                controlled_fields=None,
            ),
            shiny_tools.file_output_column_info(column_name="genbank_status",
                description=("Status of the sample submitted to GenBank."),
                controlled_fields=None,
            ),
            shiny_tools.file_output_column_info(column_name="genbank_accession",
                description=("Accession assigned to sample by GenBank."),
                controlled_fields=None,
            ),
            shiny_tools.file_output_column_info(column_name="genbank_message",
                description=("Message related to sample from GenBank."),
                controlled_fields=None,
            ),
            shiny_tools.file_output_column_info(column_name="gs-sample_name",
                description=("Sample name used for submission to GISAID."),
                controlled_fields=None,
            ),
            shiny_tools.file_output_column_info(column_name="gs-segment_name",
                description=("Sample name of influenza genome segment used for submission to GISAID."),
                controlled_fields=None,
            ),
            shiny_tools.file_output_column_info(column_name="gisaid_accession_epi_isl_id",
                description=("Accession assigned to sample by GISAID."),
                controlled_fields=None,
            ),
            shiny_tools.file_output_column_info(column_name="gisaid_accession_epi_isl_id",
                description=("Accession assigned by GISAID to sample influenza genome segment."),
                controlled_fields=None,
            ),
            shiny_tools.file_output_column_info(column_name="gisaid_message",
                description=("Message related to sample from GISAID."),
                controlled_fields=None,
            ),
        ),
    ),
    ui.hr(),
    ui.h4("BioSample Directory:"),
    ui.hr(),
    ui.h6("metadata.tsv"),
    ui.tags.ul(
        ui.p("BioSample submission metadata file. Can be used to submit your BioSample data to NCBI via their submission website instead of FTP."),
    ),
    ui.h6("report.xml"),
    ui.tags.ul(
        ui.p("Report generated by BioSample regarding the status of the current submission and sample information."),
    ),
    ui.h6("submission.xml"),
    ui.tags.ul(
        ui.p("BioSample submission xml file. Used by SeqSender to submit your BioSample data to NCBI via their FTP submission option."),
    ),
    ui.hr(),
    ui.h4("SRA Directory:"),
    ui.hr(),
    ui.h6("metadata.tsv"),
    ui.tags.ul(
        ui.p("SRA submission metadata file. Can be used to submit your SRA data to NCBI via their submission website instead of FTP."),
    ),
    ui.h6("raw_reads_location.txt"),
    ui.tags.ul(
        ui.p("SeqSender file, used to record the location of all raw reads files to be uploaded to SRA via FTP."),
    ),
    ui.h6("report.xml"),
    ui.tags.ul(
        ui.p("Report generated by SRA regarding the status of the current submission and sample information."),
    ),
    ui.h6("submission.xml"),
    ui.tags.ul(
        ui.p("SRA submission xml file. Used by SeqSender to submit your SRA data to NCBI via FTP."),
    ),
    ui.hr(),
    ui.h4("GenBank Directory:"),
    ui.hr(),
    ui.h6("AccessionReport.tsv"),
    ui.tags.ul(
        ui.p("GenBank FTP output file. Contains the status of each sample and its corrsponding accession.."),
    ),
    ui.h6("authorset.sbt"),
    ui.tags.ul(
        ui.p("GenBank submission file for author/organization information."),
    ),
    ui.h6("comment.cmt"),
    ui.tags.ul(
        ui.p("Optional tab-delimited, GenBank submission file. Contains additional sequencing and assembly information."),
    ),
    ui.h6("email.txt"),
    ui.tags.ul(
        ui.p("GenBank FTP output file. This is a copy of the email sent by GenBank to every member of your group account after your GenBank sequences have been ", ui.strong("PROCESSED"), "."),
    ),
    ui.h6("flatfile.txt"),
    ui.tags.ul(
        ui.p("GenBank FTP output file. This is a copy of the annotated sample GenBank record, in GenBank flatfile format. This is how your sample looks like live on GenBank. "),
    ),
    ui.h6("report.xml"),
    ui.tags.ul(
        ui.p("Report generated by GenBank regarding the status of the current submission."),
    ),
    ui.h6("seq-edit-report.html"),
    ui.tags.ul(
        ui.p("Optional GenBank FTP output file. If your samples were automatically modified by GenBank during the submission process to adhere to GenBank standards then the changes to each sample will be listed here. This does not mean your sample failed to upload, this usually is to just a minor typo or incorrect syntax. If your config file field ", ui.code("GenBank_Auto_Remove_Failed_Samples"), " is set to ", ui.code("True"), " then the samples automatically removed by GenBank for failing to meet their submission criteria will also be listed here with the reason."),
    ),
    ui.h6("sequence.fsa"),
    ui.tags.ul(
        ui.p("GenBank fasta file. Can be used to submit your data to GenBank via their website or it is used by SeqSender to upload via FTP or to create the tabl2asn sqn file."),
    ),
    ui.h6("source.src"),
    ui.tags.ul(
        ui.p("GenBank metadata file. Can be used to submit your data to GenBank via their website or it is used by SeqSender to upload via FTP or to create the tabl2asn sqn file."),
    ),
    ui.h6("submission.xml"),
    ui.tags.ul(
        ui.p("GenBank FTP submission xml file. Used to identify submission information and GenBank submission zip file."),
    ),
    ui.h6("<submission_name>.zip"),
    ui.tags.ul(
        ui.p("GenBank FTP submission zip file. When submitting to GenBank via FTP, this contains all of your submission files: ", ui.strong("sequence.fsa"), ", ", ui.strong("source.src"), ", ", ui.strong("authorset.sbt"), ", ", ui.strong("comment.cmt"), "."),
    ),
    ui.h6("<submission_name>.gff"),
    ui.tags.ul(
        ui.p("Copy of gff file provided via command ", ui.code("--gff_file"), " when submitting via table2asn. It is used by SeqSender to create the table2asn sqn file."),
    ),
    ui.hr(),
    ui.h4("GISAID Directory:"),
    ui.hr(),
    ui.h6("gisaid_upload_log_#.txt"),
    ui.tags.ul(
        ui.p("GISAID CLI output log where \"#\" is the numbered attempt. When SeqSender is attempting to uplo5ad to GISAID, if it fails during the submission it will attempt the submission again ", ui.strong("3"), " times. When SeqSender command ", ui.code("submission_status"), ", is ran, if the GISAID submission is still incomplete, it will attempt the process again, overwriting the previous log files until it reaches ", ui.strong("3"), " again or completes the submission."),
    ),
    ui.h6("metadata.csv"),
    ui.tags.ul(
        ui.p("GISAID metadata file. Can be used to submit your data to GISAID via their website or it is used by SeqSender to upload it to GISAID."),
    ),
    ui.h6("orig_metadata.csv"),
    ui.tags.ul(
        ui.p("Unmodified copy of the GISAID metadata file. As SeqSender uploads to GISAID, if it fails during the submission process, the samples successfully loaded to GISAID need to be removed from the ", ui.strong("metadata.csv"), " in order to reattempt uploading the remaining samples.")
    ),
    ui.h6("orig_sequence.fsa"),
    ui.tags.ul(
        ui.p("Unmodified copy of the GISAID fasta file. As SeqSender uploads to GISAID, if it fails during the submission process, the samples successfully loaded to GISAID need to be removed from the ", ui.strong("sequence.fsa"), " in order to reattempt uploading the remaining samples.")
    ),
    ui.h6("sequence.fsa"),
    ui.tags.ul(
        ui.p("GISAID fasta file. Can be used to submit your data to GISAID via their website or it is used by SeqSender to upload it to GISAID.")
    ),
]

####################### COMMANDS PAGE ###################

commands_body = [
    ui.h2("SeqSender Commands"),
    ui.p(
        "This is a list of all available commands for seqsender, usage, and available options."
    ),
    ui.accordion(
        # Prep command
        shiny_tools.command_accordion_panel("prep",
            description=" command is used to generate all of the required files for submission to the databases specified. This command is not required for performing submissions, however, it is useful for visualizing your data before submission/troubleshooting submission issues."
        ),
        # Submit command
        shiny_tools.command_accordion_panel("submit",
            description=" command is used to generate all of the required files for submission to the databases specified. Then it performs batch upload to each database based on provided config file."
        ),
        # Submission status command
        shiny_tools.command_accordion_panel("submission_status",
            description=" command is used to update the progress of files submitted via the submit command. If database selection choices require submission in a sequential order then this command will also submit files when ready."
        ),
        # Test data command
        shiny_tools.command_accordion_panel("test_data", description=" command is used generate test data for seqsender, to be used for testing the prep and submit commands."),
        # Update biosample command
        shiny_tools.command_accordion_panel("update_biosample", description=" command is used to update biosample schema options based on available BioSample Packages."),
        # version command
        shiny_tools.command_accordion_panel("version", description=" command prints the current seqsender version."),
    ),
]

###################### FAQ PAGE ##########################

faq_body = [
    ui.h2("Frequently Asked Questions:"),
]

######################## APP UI ##########################

app_ui = ui.page_fluid(
    theme.lumen(),
    ui.head_content(ui.include_css(pathlib.Path(__file__).parent / "seqsender.css")),
    ui.page_navbar(
        ui.nav_panel("SeqSender", index_body),
        ui.nav_panel("About", about_body),
        ui.nav_panel("Installation", installation_body),
        ui.nav_panel("Prerequisites", prerequisites_body),
        ui.nav_panel("Submission Wizard", setup_body),
        ui.nav_panel("My First Submission", first_submission_body),
        ui.nav_panel("Output Files", output_body),
        ui.nav_panel("Commands", commands_body),
        # ui.nav_panel("FAQ", faq_body),
        selected="SeqSender",
        header=header,
        footer=footer,
    ),
)

####################### SERVER #########################
import pathlib

dir = pathlib.Path(__file__).parent

@reactive.file_reader(dir / "templates/config.seqsender.seqsender.schema_template.csv")
def read_file():
    df = pd.read_csv(dir / "templates/config.seqsender.seqsender.schema_template.csv", index_col = "column_name")
    df = df.fillna("")
    df = df.transpose()
    return df

@reactive.file_reader(dir / "templates/config.sra.sra.schema_template.csv")
def read_sra_file():
    df = pd.read_csv(dir / "templates/config.sra.sra.schema_template.csv", index_col = "column_name")
    df = df.fillna("")
    df = df.transpose()
    return df

# Function to style metadata shiny index
def metadata_index_css(index):
    return "background-color: #f0f0f0;"

# Function to style metadata shiny column headers
def metadata_database_css(column):
    if "gb-" in column:
        # Set color of GenBank table
        return "background-color: #cce6ff;"
    elif column.startswith("gs-"):
        # Set color of GISAID table
        return "background-color: #deede8;"
    elif column.startswith("bs-"):
        # Set color of BioSample table
        return "background-color: #ffffb3;"
    elif column.startswith("sra-"):
        # Set color of SRA table
        return "background-color: #ffcce6;"
    else:
        return "background-color: #dcdcdc;"


def server(input, output, session):
    @reactive.Calc
    @reactive.file_reader(dir / "templates/")
    def read_biosample_file():
        df = pd.read_csv(dir / ("templates/config.biosample." + input.BioSample_packages() + "_template.csv"), index_col = "column_name")
        df = df.fillna("")
        df = df.transpose()
        return df

    @reactive.file_reader(dir / "templates/")
    def read_genbank_file():
        df = pd.read_csv(dir / "templates/config.genbank.genbank.schema_template.csv", index_col = "column_name")
        df = df.fillna("")
        df = df.transpose()
        return df

    @reactive.file_reader(dir / "templates/")
    def read_gisaid_file():
        df = pd.read_csv(dir / ("templates/config.gisaid.gisaid." + input.GISAID_databases() + ".schema_template.csv"), index_col = "column_name")
        df = df.fillna("")
        df = df.transpose()
        return df

    @render.text
    def BioSample_Package_Name():
        return input.BioSample_packages()

    @reactive.effect
    @reactive.event(input.SRA_checkbox)
    def sra_submission_requires_biosample():
        if input.SRA_checkbox() == True and input.BioSample_checkbox() == False:
            with reactive.isolate():
                ui.update_checkbox("BioSample_checkbox", value = True)


    @reactive.effect
    @reactive.event(input.ncbi_submission_position)
    def gisaid_submission_position():
        if input.ncbi_submission_position() == "1":
            value = "2"
            with reactive.isolate():
                ui.update_radio_buttons(
                    "gisaid_submission_position",
                    selected=value,
                )
        elif input.ncbi_submission_position() == "2":
            value = "1"
            with reactive.isolate():
                ui.update_radio_buttons(
                    "gisaid_submission_position",
                    selected=value,
                )
        else:
            value = ""
            with reactive.isolate():
                ui.update_radio_buttons(
                    "gisaid_submission_position",
                    selected=value,
                )

    @reactive.effect
    @reactive.event(input.gisaid_submission_position)
    def ncbi_submission_position():
        if input.gisaid_submission_position() == "1":
            value = "2"
            with reactive.isolate():
                ui.update_radio_buttons(
                    "ncbi_submission_position",
                    selected=value,
                )
        elif input.gisaid_submission_position() == "2":
            value = "1"
            with reactive.isolate():
                ui.update_radio_buttons(
                    "ncbi_submission_position",
                    selected=value,
                )
        else:
            value = ""
            with reactive.isolate():
                ui.update_radio_buttons(
                    "ncbi_submission_position",
                    selected=value,
                )


    @reactive.Calc
    def initialize_base_dataframe():
        return read_file()

    @reactive.Calc
    def initialize_biosample_dataframes():
        return read_biosample_file()

    @reactive.Calc
    def initialize_sra_dataframes():
        return read_sra_file()

    @reactive.Calc
    def initialize_genbank_dataframes():
        return read_genbank_file()

    @reactive.Calc
    def initialize_gisaid_dataframes():
        return read_gisaid_file()

    @output
    @render.table
    @reactive.Calc
    def load_database_metadata_dataframe():
        database_df = initialize_base_dataframe()
        if input.GenBank_checkbox():
            genbank_df = initialize_genbank_dataframes()
            database_df = pd.concat([database_df, genbank_df], axis=1)
        if input.GISAID_checkbox():
            gisaid_df = initialize_gisaid_dataframes()
            database_df = pd.concat([database_df, gisaid_df], axis=1)
        if input.BioSample_checkbox():
            biosample_df = initialize_biosample_dataframes()
            database_df = pd.concat([database_df, biosample_df], axis=1)
            database_df.loc["required_column", "bioproject"] = "Required"
        if input.SRA_checkbox():
            sra_df = initialize_sra_dataframes()
            database_df = pd.concat([database_df, sra_df], axis=1)
            database_df.loc["required_column", "bioproject"] = "Required"
        database_df = database_df.loc[:, ~database_df.columns.duplicated()]
        database_df = (
            database_df.style.set_table_styles(
                [
                    dict(
                        selector="th", props=[("border", "3px black solid !important"), ("white-space", "nowrap"), ("vertical-align", "top")]
                    ),
                    dict(selector="ti", props=[("background-color", "#DCDCDC"), ("vertical-align", "top")]),
                    dict(
                        selector="td", props=[("border", "1px black solid !important"), ("vertical-align", "top")]
                    ),
                ]
            )
            .applymap_index(metadata_database_css, axis="columns")
            .applymap_index(metadata_index_css, axis="index")
        )
        return database_df

    @reactive.Calc
    @render.download(filename="seqsender_config.yaml")
    def download_config():
        config_file = initialize_config()
        yield yaml.safe_dump(config_file, sort_keys=False).replace(r"''", '')

    @reactive.Calc
    def initialize_config():
        config_file = {"Submission": {
            **({"NCBI":{
                "Username": input.ncbi_config_username() or "",
                "Password": input.ncbi_config_password() or "",
                "Spuid_Namespace": input.ncbi_config_spuid_namespace() or "",
                **({"BioSample_Package": input.BioSample_packages() or ""} if input.BioSample_checkbox() else {}),
                **({"GenBank_Auto_Remove_Failed_Samples": input.ncbi_config_auto_remove_genbank() or ""} if input.GenBank_checkbox() else {}),
                "Publication_Title": input.ncbi_config_publication_title() or "",
                "Publication_Status": input.ncbi_config_publication_status() or "",
                **({"Submission_Position": input.ncbi_submission_position() or ""} if input.GenBank_checkbox() and input.GISAID_checkbox() else {}),
                "Specified_Release_Date": input.ncbi_config_release_date() or "",
                "Link_Sample_Between_NCBI_Databases": input.ncbi_config_link_samples() or "",
                "Description": {
                    "Organization": {
                        "Role": input.ncbi_config_role() or "",
                        "Type": input.ncbi_config_type() or "",
                        "Name": input.ncbi_config_org_name() or "",
                        "Address": {
                            "Affil": input.ncbi_config_affil() or "",
                            "Div": input.ncbi_config_div() or "",
                            "Street": input.ncbi_config_street() or "",
                            "City": input.ncbi_config_city() or "",
                            "Sub": input.ncbi_config_state() or "",
                            "Postal_Code": input.ncbi_config_postal() or "",
                            "Country": input.ncbi_config_country() or "",
                            "Email": input.ncbi_config_email() or "",
                            "Phone": input.ncbi_config_phone() or "",
                            "Submitter": {
                                "Email": input.ncbi_sub_email_one() or "",
                                "Alt_Email": input.ncbi_sub_email_two() or "",
                                "Name": {
                                    "First": input.ncbi_config_first_name() or "",
                                    "Last": input.ncbi_config_last_name() or "",
                                    }
                                }
                            }
                        }
                    }
                }} if input.BioSample_checkbox() or input.SRA_checkbox() or input.GenBank_checkbox() else {}),
            **({"GISAID": {
                "Client-Id": input.gisaid_config_client() or "",
                "Username": input.gisaid_config_username() or "",
                "Password": input.gisaid_config_password() or "",
                **({"Submission_Position": input.gisaid_submission_position() or ""} if input.GenBank_checkbox() and input.GISAID_checkbox() else {}),
                }} if input.GISAID_checkbox() else {})
            }
        }
        config_file = {key: (None if value == "" else value) for key, value in config_file.items()}
        return config_file

    @reactive.Calc
    @render.download(filename="metadata_template.csv")
    def download_metadata():
        database_df = initialize_base_dataframe()
        if input.GenBank_checkbox():
            genbank_df = initialize_genbank_dataframes()
            database_df = pd.concat([database_df, genbank_df], axis=1)
        if input.GISAID_checkbox():
            gisaid_df = initialize_gisaid_dataframes()
            database_df = pd.concat([database_df, gisaid_df], axis=1)
        if input.BioSample_checkbox():
            biosample_df = initialize_biosample_dataframes()
            database_df = pd.concat([database_df, biosample_df], axis=1)
            database_df.loc["required_column", "bioproject"] = "Required"
        if input.SRA_checkbox():
            sra_df = initialize_sra_dataframes()
            database_df = pd.concat([database_df, sra_df], axis=1)
            database_df.loc["required_column", "bioproject"] = "Required"
        database_df = database_df.loc[:,~database_df.columns.duplicated()].copy()
        yield database_df.to_csv()


app = App(app_ui, server)
