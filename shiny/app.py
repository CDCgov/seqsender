from shiny import App, Inputs, Outputs, Session, render, req, ui, reactive
from htmltools import TagList, div

###################### CSS #######################
import pandas as pd
from shinyswatch import theme
import shiny_tools
from index import index_body
from setup import setup_body
import pathlib

terminal_css = ""
yaml_css = "background-color: #F0F0F0;white-space: nowrap; font-size: 20px ;margin-top:-15px;font-family: Consolas,Monaco,Lucida Console,Liberation Mono,DejaVu Sans Mono,Bitstream Vera Sans Mono,Courier New, monospace;-webkit-user-select: none; -ms-user-select: none; user-select: none;"
############ HEADER #######################
header = (
    ui.card_header(
        ui.HTML(
            """<p><strong>Beta Version</strong>: 1.2.0. This pipeline is currently in Beta testing, and issues could appear during submission. Please use it at your own risk. Feedback and suggestions are welcome!</p>"""
        ),
        ui.card(
            ui.strong("Documentation Page is Under Development. Not all buttons are functioning correctly/some documentation is not yet available. This should be resolved shortly."),
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


####################### PREREQUISITES PAGE ###########################
ncbi_prereq = [
    ui.HTML(
        """
<strong>NCBI Submissions</strong>
<p><code>seqsender</code> utilizes an UI-Less Data Submission Protocol to bulk upload submission files (e.g., <em>submission.xml</em>, <em>submission.zip</em>, etc.) to NCBI archives. The submission files are uploaded to the NCBI server via FTP on the command line. Before attempting to submit a submission using <code>seqsender</code>, submitter will need to</p>
<ol style="list-style-type: decimal">
<li><p>Have a NCBI account. To sign up, visit <a href="https://account.ncbi.nlm.nih.gov/">NCBI website</a>.</p></li>
<li><p>Required for CDC users and highly recommended for others is creating a center account for your institution/lab <a href="https://submit.ncbi.nlm.nih.gov/sarscov2/sra/#step6">NCBI Center Account Instructions</a>. Center accounts allow you to perform submissions UI-less submissions as your institution/lab.</p></li>
<li><p>Required for CDC users and also recommended is creating a submission group in <a href="https://submit.ncbi.nlm.nih.gov">NCBI Submission Portal</a>. A group should include all individuals who need access to UI-less submissions through the web interface with your center account. Each member of the group must also have an individual NCBI account. <a href="https://account.ncbi.nlm.nih.gov/">NCBI website</a>.</p></li>
<li><p>Refer to this page for information regarding requirements for GenBank submissions via FTP only. This page applies only for COVID and Influenza <a href="https://submit.ncbi.nlm.nih.gov/sarscov2/genbank/#step5">NCBI GenBank FTP Submissions</a> For further questions contact <a href="mailto:gb-admin@ncbi.nlm.nih.gov">gb-admin@ncbi.nlm.nih.gov</a> to discuss requirements for submissions.</p></li>
<li><p>Coordinate a NCBI namespace name (<strong>spuid_namespace</strong>)</p></li>
</ol>
"""
    ),
]
gisaid_prereq = [
    ui.HTML(
        """
<strong>GISAID Submissions</strong>
<p><code>seqsender</code> makes use of GISAID’s Command Line Interface tools to bulk uploading meta- and sequence-data to GISAID databases. Presently, the pipeline only allows upload to EpiFlu (<strong>Influenza A Virus</strong>) and EpiCoV (<strong>SARS-COV-2</strong>) databases. Before uploading, submitter needs to</p>
<ol style="list-style-type: decimal">
<li><p>Have a GISAID account. To sign up, visit <a href="https://gisaid.org/">GISAID Platform</a>.</p></li>
<li><p>Request a client-ID for EpiFlu or EpiCoV database in order to use its CLI tool. The CLI utilizes the client-ID along with the username and password to authenticate the database prior to make a submission. To obtain a client-ID, please email <a href="mailto:clisupport@gisaid.org">clisupport@gisaid.org</a> to request. <em><strong>Important note</strong>: If submitter would like to upload a “test” submission first to familiarize themselves with the submission process prior to make a real submission, one should additionally request a test client-id to perform such submissions.</em></p></li>
<li><p>Download the <a href="https://cdcgov.github.io/seqsender/articles/images/fluCLI_download.png" target="_blank">EpiFlu</a> or <a href="https://cdcgov.github.io/seqsender/articles/images/covCLI_download.png" target="_blank">EpiCoV</a> CLI from the <strong>GISAID platform</strong> and stored them in the destination of choice prior to perform a batch upload.</p></li>
</ol>
"""
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
####################### INSTALLATION PAGE ###########################
#### LOCAL INSTALLATION
gisaid_installation_content = [
    ui.HTML(
        """<p>Here is a quick look of where to store the downloaded <strong>GISAID CLI</strong> package.</p>
"""
    ),
]
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
        ui.nav_panel("Docker-Compose", "Panel A content"),
        id="docker_options",
    ),
    ui.HTML(
        """
<h2>check if the container is created successfully</h2>
<pre><code>docker container ps</code></pre>
"""
    ),
    ui.card(
        ui.HTML(
            """
docker container ps<br>
CONTAINER ID&nbsp;|&nbsp;&nbsp;|&nbsp;IMAGE&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;COMMAND&nbsp;|&nbsp;&nbsp;|&nbsp;CREATED&nbsp;|&nbsp;&nbsp;STATUS&nbsp;|&nbsp;PORTS&nbsp;|&nbsp;NAMES
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
        ui.nav_panel("Local", local_installation_content + gisaid_installation_content),
        ui.nav_panel(
            "Docker", docker_installation_content + gisaid_installation_content
        ),
        ui.nav_panel(
            "Singularity",
            singularity_installation_content + gisaid_installation_content,
        ),
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
    ),
]

####################### COMMANDS PAGE ###################

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
        ui.strong(ui.code("--organism "), "{'FLU', 'COV', 'POX', 'ARBO', 'OTHER'}"),
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

gff_parameter = ui.div(
    ui.strong(ui.code("--gff_file")),
    ui.tags.ul("Full path to gff file if not stored in ", ui.code("--submission_dir"), " location."),
)

test_parameter = ui.div(
    ui.strong(ui.code("--test")),
    ui.tags.ul("Flag to enable a \"test\" submission to each database selected."),
)

commands_body = [
    ui.h2("SeqSender Commands"),
    ui.p(
        "This is a list of all available commands for seqsender, usage, and available options."
    ),
    ui.accordion(
        ui.accordion_panel(
            "Prep",
            ui.p(
                "The ",
                ui.code("prep"),
                " command is used to generate all of the required files for submission to the databases specified. This command is not required for performing submissions, however, it is useful for visualizing your data before submission/troubleshooting submission issues.",
            ),
            shiny_tools.seqsender_example_command("local", "prep", ""),
            ui.accordion(
                ui.accordion_panel(
                    "Parameters",
                    databases_parameter,
                    organisms_parameter,
                    submission_directory_parameter,
                    submission_name_parameter,
                    config_parameter,
                    metadata_parameter,
                    fasta_parameter,
                    gff_parameter,
                    test_parameter,
                ),
                ui.accordion_panel(
                    "Notes"
                ),
                ui.accordion_panel(
                    "Examples"
                ),
                open=["Parameters", "Notes"],
            ),
        ),
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
        ui.nav_panel("Prerequisites", prerequisites_body),
        ui.nav_panel("Installation", installation_body),
        ui.nav_panel("Setup", setup_body),
        ui.nav_panel("First Submission", testing_body),
        ui.nav_panel("Commands", commands_body),
        ui.nav_panel("FAQ", faq_body),
        selected="SeqSender",
        header=header,
        footer=footer,
    ),
)

####################### SERVER #########################
import pathlib

dir = pathlib.Path(__file__).parent


@reactive.file_reader(dir / "templates/config.seqsender.schema_template.csv")
def read_file():
    df = pd.read_csv(dir / "templates/config.seqsender.schema_template.csv", index_col = "column_name")
    df = df.fillna("")
    df = df.transpose()
    return df

@reactive.file_reader(dir / "templates/config.biosample.Pathogen.cl.1.0_template.csv")
def read_biosample_file():
    df = pd.read_csv(dir / "templates/config.biosample.Pathogen.cl.1.0_template.csv", index_col = "column_name")
    df = df.fillna("")
    df = df.transpose()
    return df

@reactive.file_reader(dir / "templates/config.sra.schema_template.csv")
def read_sra_file():
    df = pd.read_csv(dir / "templates/config.sra.schema_template.csv", index_col = "column_name")
    df = df.fillna("")
    df = df.transpose()
    return df

@reactive.file_reader(dir / "templates/config.genbank.genbank.schema_template.csv")
def read_genbank_file():
    df = pd.read_csv(dir / "templates/config.genbank.genbank.schema_template.csv", index_col = "column_name")
    df = df.fillna("")
    df = df.transpose()
    return df

@reactive.file_reader(dir / "templates/config.gisaid.gisaid.FLU.schema_template.csv")
def read_gisaid_file():
    df = pd.read_csv(dir / "templates/config.gisaid.gisaid.FLU.schema_template.csv", index_col = "column_name")
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
    @render.text
    def BioSample_Package_Name():
        return input.BioSample_packages()

    @reactive.effect
    @reactive.event(input.ncbi_submission_position)
    def gisaid_submission_position():
        if "First" in input.ncbi_submission_position():
            value = "Second"
        else:
            value = "First"
        with reactive.isolate():
            ui.update_radio_buttons(
                "gisaid_submission_position",
                selected=value,
            )

    @reactive.effect
    @reactive.event(input.gisaid_submission_position)
    def ncbi_submission_position():
        if "First" in input.gisaid_submission_position():
            value = "Second"
        else:
            value = "First"
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
        if input.BioSample_checkbox():
            biosample_df = initialize_biosample_dataframes()
            database_df = pd.concat([database_df, biosample_df], axis=1)
        if input.SRA_checkbox():
            sra_df = initialize_sra_dataframes()
            database_df = pd.concat([database_df, sra_df], axis=1)
        if input.GenBank_checkbox():
            genbank_df = initialize_genbank_dataframes()
            database_df = pd.concat([database_df, genbank_df], axis=1)
        if input.GISAID_checkbox():
            gisaid_df = initialize_gisaid_dataframes()
            database_df = pd.concat([database_df, gisaid_df], axis=1)
        database_df = (
            database_df.style.set_table_styles(
                [
                    dict(
                        selector="th", props=[("border", "3px black solid !important")]
                    ),
                    dict(selector="ti", props=[("background-color", "#DCDCDC")]),
                    dict(
                        selector="td", props=[("border", "1px black solid !important")]
                    ),
                ]
            )
            .applymap_index(metadata_database_css, axis="columns")
            .applymap_index(metadata_index_css, axis="index")
        )
        return database_df

    @reactive.Calc
    @render.download(filename="metadata_template.csv")
    def download_metadata():
        database_df = initialize_base_dataframe()
        if input.BioSample_checkbox():
            biosample_df = initialize_biosample_dataframes()
            database_df = pd.concat([database_df, biosample_df], axis=1)
        if input.SRA_checkbox():
            sra_df = initialize_sra_dataframes()
            database_df = pd.concat([database_df, sra_df], axis=1)
        if input.GenBank_checkbox():
            genbank_df = initialize_genbank_dataframes()
            database_df = pd.concat([database_df, genbank_df], axis=1)
        if input.GISAID_checkbox():
            gisaid_df = initialize_gisaid_dataframes()
            database_df = pd.concat([database_df, gisaid_df], axis=1)
        yield database_df.to_csv()


app = App(app_ui, server)
