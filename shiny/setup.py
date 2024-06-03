from shiny import App, Inputs, Outputs, Session, render, req, ui, reactive
from htmltools import TagList, div
####################### SETUP PAGE ###############################
yaml_css = "background-color: #F0F0F0;white-space: nowrap; font-size: 20px ;margin-top:-15px;font-family: Consolas,Monaco,Lucida Console,Liberation Mono,DejaVu Sans Mono,Bitstream Vera Sans Mono,Courier New, monospace;-webkit-user-select: none; -ms-user-select: none; user-select: none;"

setup_body = [
    ui.h2("Setup"),
    ui.layout_columns(
        ui.card(
            ui.card_header(
                ui.h4("Select Database:", style="display:inline-block;"),
                ui.tooltip(
                    ui.p("❓"),
                    "Additional info",
                    id="select_database_tooltip",
                    style="display:inline-block;float:right;",
                ),
            ),
            ui.input_checkbox(
                "BioSample_checkbox", "BioSample", width=None, value=True
            ),
            ui.input_checkbox("SRA_checkbox", "SRA", width=None),
            ui.input_checkbox("GenBank_checkbox", "GenBank", width=None),
            ui.input_checkbox("GISAID_checkbox", "GISAID", width=None, value=True),
        ),
        ui.card(
            ui.card_header(
                ui.h3("Database Options:", style="display:inline-block;"),
                ui.tooltip(
                    ui.p("❓"),
                    "Additional info",
                    id="database_options_tooltip",
                    style="display:inline-block;float:right;",
                ),
            ),
            # If BioSample checkbox checked then load BioSample Package Options
            ui.panel_conditional(
                "input.BioSample_checkbox",
                ui.input_select(
                    "BioSample_packages",
                    label="Select BioSample Package:",
                    choices=["SARS", "PATHOGEN"],
                    selected="SARS",
                ),
            ),
            # If GenBank checkbox checked then load Genbank Schema Options
            ui.panel_conditional(
                "input.GenBank_checkbox",
                ui.input_select(
                    "GenBank_schemas",
                    label="Select GenBank Schema:",
                    choices=["COV", "FLU", "OTHER"],
                ),
            ),
            # If GISAID checkbox checked then load GISAID Database Options
            ui.panel_conditional(
                "input.GISAID_checkbox",
                ui.input_select(
                    "GISAID_databases",
                    label="Select GISAID Database:",
                    choices=["COV", "FLU", "OTHER"],
                ),
            ),
        ),
    ),
    # Based on choices display config and metadata template
    ui.panel_conditional(
        "input.BioSample_checkbox || input.SRA_checkbox || input.GenBank_checkbox || input.GISAID_checkbox",
        ui.card(
            ui.card_header(
                ui.h4("Create Config File:", style="display:inline-block;"),
                ui.tooltip(
                    ui.p("❓"),
                    "Additional info",
                    id="config_tooltip",
                    style="display:inline-block;float:right;",
                ),
            ),
            ui.card(
                ui.page_fluid(
                    div(ui.HTML("Submission:"), style="margin-top:-20px"),
                    ui.panel_conditional(
                        "input.BioSample_checkbox || input.SRA_checkbox || input.GenBank_checkbox",
                        div(ui.HTML("&nbsp;|&nbsp;NCBI:")),
                        div(
                            ui.HTML("&nbsp;|&nbsp;&nbsp;|&nbsp;Username:"),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_text(
                                "ncbi_config_username",
                                label=None,
                                placeholder="NCBI FTP Username",
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        div(
                            ui.HTML("&nbsp;|&nbsp;&nbsp;|&nbsp;Password:"),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_password(
                                "ncbi_config_password",
                                label=None,
                                placeholder="NCBI FTP Password",
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        div(
                            ui.HTML("&nbsp;|&nbsp;&nbsp;|&nbsp;BioSample_Package:"),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.output_text(
                                "BioSample_Package_Name",
                                inline=True,
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        div(
                            ui.HTML("&nbsp;|&nbsp;&nbsp;|&nbsp;Submission_Position:"),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_radio_buttons(
                                "ncbi_submission_position",
                                label=None,
                                choices=["First", "Second"],
                                inline=True,
                            ),
                            style="display:inline-block;height:5px;font-size:medium;",
                        ),
                        div(ui.HTML("&nbsp;|&nbsp;&nbsp;|&nbsp;Description:")),
                        div(
                            ui.HTML("&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;Title:"),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_text(
                                "ncbi_config_title",
                                label=None,
                                placeholder="Submission Title",
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        div(
                            ui.HTML("&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;Comment:"),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_text(
                                "ncbi_config_comment",
                                label=None,
                                placeholder="Submission Comment",
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        div(
                            ui.HTML(
                                "&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;Organization:"
                            )
                        ),
                        div(
                            ui.HTML(
                                "&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;@role:"
                            ),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_text(
                                "ncbi_config_role",
                                label=None,
                                placeholder="Organization Role",
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        div(
                            ui.HTML(
                                "&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;@type:"
                            ),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_text(
                                "ncbi_config_type",
                                label=None,
                                placeholder="Organization Type",
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        div(
                            ui.HTML(
                                "&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;Name:"
                            ),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_text(
                                "ncbi_config_org_name",
                                label=None,
                                placeholder="Organization Name",
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        div(
                            ui.HTML(
                                "&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;Address:"
                            )
                        ),
                        div(
                            ui.HTML(
                                "&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;Affil:"
                            ),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_text(
                                "ncbi_config_affil",
                                label=None,
                                placeholder="Affiliation",
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        div(
                            ui.HTML(
                                "&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;Div:"
                            ),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_text(
                                "ncbi_config_div", label=None, placeholder="Division"
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        div(
                            ui.HTML(
                                "&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;Street:"
                            ),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_text(
                                "ncbi_config_street",
                                label=None,
                                placeholder="Organization Street Name",
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        div(
                            ui.HTML(
                                "&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;City:"
                            ),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_text(
                                "ncbi_config_city",
                                label=None,
                                placeholder="Organization City",
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        div(
                            ui.HTML(
                                "&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;Sub:"
                            ),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_text(
                                "ncbi_config_state",
                                label=None,
                                placeholder="Organization State/Providence",
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        div(
                            ui.HTML(
                                "&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;Postal_code:"
                            ),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_text(
                                "ncbi_config_postal",
                                label=None,
                                placeholder="Organization Postal Code",
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        div(
                            ui.HTML(
                                "&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;Country:"
                            ),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_text(
                                "ncbi_config_country",
                                label=None,
                                placeholder="Organization Country",
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        div(
                            ui.HTML(
                                "&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;Email:"
                            ),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_text(
                                "ncbi_config_email",
                                label=None,
                                placeholder="Organization Contact Email",
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        div(
                            ui.HTML(
                                "&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;Phone:"
                            ),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_text(
                                "ncbi_config_phone",
                                label=None,
                                placeholder="Organization Contact Phone",
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        div(
                            ui.HTML(
                                "&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;Submitter:"
                            )
                        ),
                        div(
                            ui.HTML(
                                "&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;@email:"
                            ),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_text(
                                "ncbi_config_sub_email1",
                                label=None,
                                placeholder="Submitter Contact Email",
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        div(
                            ui.HTML(
                                "&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;@alt_email:"
                            ),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_text(
                                "ncbi_config_sub_email2",
                                label=None,
                                placeholder="Alternate Submitter Contact Email",
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        div(
                            ui.HTML(
                                "&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;Name:"
                            )
                        ),
                        div(
                            ui.HTML(
                                "&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;First:"
                            ),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_text(
                                "ncbi_config_first_name",
                                label=None,
                                placeholder="Submitter First Name",
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        div(
                            ui.HTML(
                                "&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;&nbsp;|&nbsp;Last:"
                            ),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_text(
                                "ncbi_config_last_name",
                                label=None,
                                placeholder="Submitter Last Name",
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                    ),
                    ui.panel_conditional(
                        "input.GISAID_checkbox",
                        div(ui.HTML("&nbsp;|&nbsp;GISAID:<br>")),
                        div(
                            ui.HTML("&nbsp;|&nbsp;&nbsp;|&nbsp;Client-Id:"),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_text(
                                "gisaid_config_client",
                                label=None,
                                placeholder="GISAID Client ID",
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        div(
                            ui.HTML("&nbsp;|&nbsp;&nbsp;|&nbsp;Username:"),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_text(
                                "gisaid_config_username",
                                label=None,
                                placeholder="GISAID Username",
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        div(
                            ui.HTML("&nbsp;|&nbsp;&nbsp;|&nbsp;Password:"),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_text(
                                "gisaid_config_password",
                                label=None,
                                placeholder="GISAID Password",
                            ),
                            style="display:inline-block;height:5px;",
                        ),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        div(
                            ui.HTML("&nbsp;|&nbsp;&nbsp;|&nbsp;Submission_Position:"),
                            style="display:inline-block;",
                        ),
                        div(
                            ui.input_radio_buttons(
                                "gisaid_submission_position",
                                label=None,
                                choices=["First", "Second"],
                                inline=True,
                                selected="Second",
                            ),
                            style="display:inline-block;height:5px;font-size:medium;",
                        ),
                    ),
                ),
                style=yaml_css,
            ),
            ui.download_button("download_config", "Download Config File"),
        ),
        ui.card(
            ui.card_header(
                ui.h4("Metadata Template:", style="display:inline-block;"),
                ui.tooltip(
                    ui.p("❓"),
                    "Additional info",
                    id="metadata_tooltip",
                    style="display:inline-block;float:right;",
                ),
            ),
            ui.HTML(
                "<p><strong>Note:</strong> More columns may be available not displayed in this list. To add a column not here during submission, add it matching the style of the other columns related to the database: <code>database_prefix-column_name</code>."
            ),
            ui.output_table("load_database_metadata_dataframe"),
            ui.HTML(
                """
<h8>Metadata Sections Legend:</h8>
<p><span style="color:#696969"> SeqSender: No prefix</span> | <span style="color:#b3b300">BioSample: bs-</span> | <span style="color:#b3005c">SRA: sra-</span> | <span style="color:#005cb3">GenBank: gb-</span> | <span style="color:#3f7362">GISAID: gs-</span></p>
"""
            ),
            ui.download_button("download_metadata", "Download Metadata Template"),
        ),
    ),
]
