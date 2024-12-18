from shiny import App, Inputs, Outputs, Session, render, req, ui, reactive
from htmltools import TagList, div
import shiny_tools
####################### SETUP PAGE ###############################
yaml_css = "background-color: #F0F0F0;white-space: nowrap; font-size: 20px ;margin-top:-15px;font-family: Consolas,Monaco,Lucida Console,Liberation Mono,DejaVu Sans Mono,Bitstream Vera Sans Mono,Courier New, monospace;-webkit-user-select: none; -ms-user-select: none; user-select: none;"

def config_indent(tabs, string, custom_style = "display:inline-block;"):
    return div(ui.HTML(("&nbsp;|&nbsp;" * tabs) + string), style=custom_style)

def config_text_input(id, placeholder, help_msg, custom_style = "display:inline-block;height:5px;"):
    return div(ui.tooltip(ui.input_text(id, label=None, placeholder=placeholder), help_msg, id=(id + "_tooltip")), style=custom_style)

setup_body = [
    ui.h2("Setup"),
    ui.layout_columns(
        ui.card(
            ui.card_header(
                ui.h4("Select Database:", style="display:inline-block;"),
                shiny_tools.create_help_tooltip(
                    id = "db_select_tooltip",
                    description = "Select all databases you plan to submit your data to. The databases you plan to submit to will change the requirements for the config and metadata file.",
                    position="right",
                ),
            ),
            ui.input_checkbox("BioSample_checkbox", "BioSample", width=None),
            ui.tooltip(ui.input_checkbox("SRA_checkbox", "SRA", width=None),
                "When submitting to SRA, a BioSample submission is also required.",
                style="display:inline-block;float:right;"),
            ui.input_checkbox("GenBank_checkbox", "GenBank", width=None),
            ui.input_checkbox("GISAID_checkbox", "GISAID", width=None),
        ),
        ui.card(
            ui.card_header(
                ui.h3("Database Options:", style="display:inline-block;"),
                shiny_tools.create_help_tooltip(id = "db_option_tooltip",
                    description = "Select database dependent options. This changes metadata requirements based on database submission criteria.",
                    position = "right",
                ),
            ),
            # If BioSample checkbox checked then load BioSample Package Options
            ui.panel_conditional(
                "input.BioSample_checkbox",
                ui.input_select(
                    "BioSample_packages",
                    label="Select BioSample Package:",
                    choices={"STANDARD": {"SARS-CoV-2.cl.1.0":"SARS-CoV-2.cl.1.0","SARS-CoV-2.wwsurv.1.0":"SARS-CoV-2.wwsurv.1.0","Virus.1.0":"Virus.1.0","Beta-lactamase.1.0":"Beta-lactamase.1.0","Human.1.0":"Human.1.0","Invertebrate.1.0":"Invertebrate.1.0","Metagenome.environmental.1.0":"Metagenome.environmental.1.0","Microbe.1.0":"Microbe.1.0","Model.organism.animal.1.0":"Model.organism.animal.1.0","OneHealthEnteric.1.0":"OneHealthEnteric.1.0","Plant.1.0":"Plant.1.0"},
                        "PATHOGEN": {"Pathogen.cl.1.0":"Pathogen.cl.1.0","Pathogen.env.1.0":"Pathogen.env.1.0"},
                        "MIGS.ba": {"MIGS.ba.6.0":"MIGS.ba.6.0","MIGS.ba.agriculture.6.0":"MIGS.ba.agriculture.6.0","MIGS.ba.air.6.0":"MIGS.ba.air.6.0","MIGS.ba.built.6.0":"MIGS.ba.built.6.0","MIGS.ba.food-animal.6.0":"MIGS.ba.food-animal.6.0","MIGS.ba.food-farm.env.6.0":"MIGS.ba.food-farm.env.6.0","MIGS.ba.food-human.foods.6.0":"MIGS.ba.food-human.foods.6.0","MIGS.ba.food-prod.facility.6.0":"MIGS.ba.food-prod.facility.6.0","MIGS.ba.host-associated.6.0":"MIGS.ba.host-associated.6.0","MIGS.ba.human-associated.6.0":"MIGS.ba.human-associated.6.0","MIGS.ba.human-gut.6.0":"MIGS.ba.human-gut.6.0","MIGS.ba.human-oral.6.0":"MIGS.ba.human-oral.6.0","MIGS.ba.human-skin.6.0":"MIGS.ba.human-skin.6.0","MIGS.ba.human-vaginal.6.0":"MIGS.ba.human-vaginal.6.0","MIGS.ba.hydrocarbon-cores.6.0":"MIGS.ba.hydrocarbon-cores.6.0","MIGS.ba.hydrocarbon-fluids.swabs.6.0":"MIGS.ba.hydrocarbon-fluids.swabs.6.0","MIGS.ba.microbial.6.0":"MIGS.ba.microbial.6.0","MIGS.ba.miscellaneous.6.0":"MIGS.ba.miscellaneous.6.0","MIGS.ba.plant-associated.6.0":"MIGS.ba.plant-associated.6.0","MIGS.ba.sediment.6.0":"MIGS.ba.sediment.6.0","MIGS.ba.soil.6.0":"MIGS.ba.soil.6.0","MIGS.ba.symbiont-associated.6.0":"MIGS.ba.symbiont-associated.6.0","MIGS.ba.wastewater.6.0":"MIGS.ba.wastewater.6.0","MIGS.ba.water.6.0":"MIGS.ba.water.6.0"},
                        "MIGS.eu": {"MIGS.eu.6.0":"MIGS.eu.6.0","MIGS.eu.agriculture.6.0":"MIGS.eu.agriculture.6.0","MIGS.eu.air.6.0":"MIGS.eu.air.6.0","MIGS.eu.built.6.0":"MIGS.eu.built.6.0","MIGS.eu.food-animal.6.0":"MIGS.eu.food-animal.6.0","MIGS.eu.food-farm.env.6.0":"MIGS.eu.food-farm.env.6.0","MIGS.eu.food-human.foods.6.0":"MIGS.eu.food-human.foods.6.0","MIGS.eu.food-prod.facility.6.0":"MIGS.eu.food-prod.facility.6.0","MIGS.eu.host-associated.6.0":"MIGS.eu.host-associated.6.0","MIGS.eu.human-associated.6.0":"MIGS.eu.human-associated.6.0","MIGS.eu.human-gut.6.0":"MIGS.eu.human-gut.6.0","MIGS.eu.human-oral.6.0":"MIGS.eu.human-oral.6.0","MIGS.eu.human-skin.6.0":"MIGS.eu.human-skin.6.0","MIGS.eu.human-vaginal.6.0":"MIGS.eu.human-vaginal.6.0","MIGS.eu.hydrocarbon-cores.6.0":"MIGS.eu.hydrocarbon-cores.6.0","MIGS.eu.hydrocarbon-fluids.swabs.6.0":"MIGS.eu.hydrocarbon-fluids.swabs.6.0","MIGS.eu.microbial.6.0":"MIGS.eu.microbial.6.0","MIGS.eu.miscellaneous.6.0":"MIGS.eu.miscellaneous.6.0","MIGS.eu.plant-associated.6.0":"MIGS.eu.plant-associated.6.0","MIGS.eu.sediment.6.0":"MIGS.eu.sediment.6.0","MIGS.eu.soil.6.0":"MIGS.eu.soil.6.0","MIGS.eu.symbiont-associated.6.0":"MIGS.eu.symbiont-associated.6.0","MIGS.eu.wastewater.6.0":"MIGS.eu.wastewater.6.0","MIGS.eu.water.6.0":"MIGS.eu.water.6.0"},
                        "MIGS.vi": {"MIGS.vi.6.0":"MIGS.vi.6.0","MIGS.vi.agriculture.6.0":"MIGS.vi.agriculture.6.0","MIGS.vi.air.6.0":"MIGS.vi.air.6.0","MIGS.vi.built.6.0":"MIGS.vi.built.6.0","MIGS.vi.food-animal.6.0":"MIGS.vi.food-animal.6.0","MIGS.vi.food-farm.env.6.0":"MIGS.vi.food-farm.env.6.0","MIGS.vi.food-human.foods.6.0":"MIGS.vi.food-human.foods.6.0","MIGS.vi.food-prod.facility.6.0":"MIGS.vi.food-prod.facility.6.0","MIGS.vi.host-associated.6.0":"MIGS.vi.host-associated.6.0","MIGS.vi.human-associated.6.0":"MIGS.vi.human-associated.6.0","MIGS.vi.human-gut.6.0":"MIGS.vi.human-gut.6.0","MIGS.vi.human-oral.6.0":"MIGS.vi.human-oral.6.0","MIGS.vi.human-skin.6.0":"MIGS.vi.human-skin.6.0","MIGS.vi.human-vaginal.6.0":"MIGS.vi.human-vaginal.6.0","MIGS.vi.hydrocarbon-cores.6.0":"MIGS.vi.hydrocarbon-cores.6.0","MIGS.vi.hydrocarbon-fluids.swabs.6.0":"MIGS.vi.hydrocarbon-fluids.swabs.6.0","MIGS.vi.microbial.6.0":"MIGS.vi.microbial.6.0","MIGS.vi.miscellaneous.6.0":"MIGS.vi.miscellaneous.6.0","MIGS.vi.plant-associated.6.0":"MIGS.vi.plant-associated.6.0","MIGS.vi.sediment.6.0":"MIGS.vi.sediment.6.0","MIGS.vi.soil.6.0":"MIGS.vi.soil.6.0","MIGS.vi.symbiont-associated.6.0":"MIGS.vi.symbiont-associated.6.0","MIGS.vi.wastewater.6.0":"MIGS.vi.wastewater.6.0","MIGS.vi.water.6.0":"MIGS.vi.water.6.0"},
                        "MIMAG": {"MIMAG.6.0":"MIMAG.6.0","MIMAG.agriculture.6.0":"MIMAG.agriculture.6.0","MIMAG.air.6.0":"MIMAG.air.6.0","MIMAG.built.6.0":"MIMAG.built.6.0","MIMAG.food-animal.6.0":"MIMAG.food-animal.6.0","MIMAG.food-farm.env.6.0":"MIMAG.food-farm.env.6.0","MIMAG.food-human.foods.6.0":"MIMAG.food-human.foods.6.0","MIMAG.food-prod.facility.6.0":"MIMAG.food-prod.facility.6.0","MIMAG.host-associated.6.0":"MIMAG.host-associated.6.0","MIMAG.human-associated.6.0":"MIMAG.human-associated.6.0","MIMAG.human-gut.6.0":"MIMAG.human-gut.6.0","MIMAG.human-oral.6.0":"MIMAG.human-oral.6.0","MIMAG.human-skin.6.0":"MIMAG.human-skin.6.0","MIMAG.human-vaginal.6.0":"MIMAG.human-vaginal.6.0","MIMAG.hydrocarbon-cores.6.0":"MIMAG.hydrocarbon-cores.6.0","MIMAG.hydrocarbon-fluids.swabs.6.0":"MIMAG.hydrocarbon-fluids.swabs.6.0","MIMAG.microbial.6.0":"MIMAG.microbial.6.0","MIMAG.miscellaneous.6.0":"MIMAG.miscellaneous.6.0","MIMAG.plant-associated.6.0":"MIMAG.plant-associated.6.0","MIMAG.sediment.6.0":"MIMAG.sediment.6.0","MIMAG.soil.6.0":"MIMAG.soil.6.0","MIMAG.symbiont-associated.6.0":"MIMAG.symbiont-associated.6.0","MIMAG.wastewater.6.0":"MIMAG.wastewater.6.0","MIMAG.water.6.0":"MIMAG.water.6.0"},
                        "MIMARKS.specimen": {"MIMARKS.specimen.6.0":"MIMARKS.specimen.6.0","MIMARKS.specimen.agriculture.6.0":"MIMARKS.specimen.agriculture.6.0","MIMARKS.specimen.air.6.0":"MIMARKS.specimen.air.6.0","MIMARKS.specimen.built.6.0":"MIMARKS.specimen.built.6.0","MIMARKS.specimen.food-animal.6.0":"MIMARKS.specimen.food-animal.6.0","MIMARKS.specimen.food-farm.env.6.0":"MIMARKS.specimen.food-farm.env.6.0","MIMARKS.specimen.food-human.foods.6.0":"MIMARKS.specimen.food-human.foods.6.0","MIMARKS.specimen.food-prod.facility.6.0":"MIMARKS.specimen.food-prod.facility.6.0","MIMARKS.specimen.host-associated.6.0":"MIMARKS.specimen.host-associated.6.0","MIMARKS.specimen.human-associated.6.0":"MIMARKS.specimen.human-associated.6.0","MIMARKS.specimen.human-gut.6.0":"MIMARKS.specimen.human-gut.6.0","MIMARKS.specimen.human-oral.6.0":"MIMARKS.specimen.human-oral.6.0","MIMARKS.specimen.human-skin.6.0":"MIMARKS.specimen.human-skin.6.0","MIMARKS.specimen.human-vaginal.6.0":"MIMARKS.specimen.human-vaginal.6.0","MIMARKS.specimen.hydrocarbon-cores.6.0":"MIMARKS.specimen.hydrocarbon-cores.6.0","MIMARKS.specimen.hydrocarbon-fluids.swabs.6.0":"MIMARKS.specimen.hydrocarbon-fluids.swabs.6.0","MIMARKS.specimen.microbial.6.0":"MIMARKS.specimen.microbial.6.0","MIMARKS.specimen.miscellaneous.6.0":"MIMARKS.specimen.miscellaneous.6.0","MIMARKS.specimen.plant-associated.6.0":"MIMARKS.specimen.plant-associated.6.0","MIMARKS.specimen.sediment.6.0":"MIMARKS.specimen.sediment.6.0","MIMARKS.specimen.soil.6.0":"MIMARKS.specimen.soil.6.0","MIMARKS.specimen.symbiont-associated.6.0":"MIMARKS.specimen.symbiont-associated.6.0","MIMARKS.specimen.wastewater.6.0":"MIMARKS.specimen.wastewater.6.0","MIMARKS.specimen.water.6.0":"MIMARKS.specimen.water.6.0"},
                        "MIMARKS.survey": {"MIMARKS.survey.agriculture.6.0":"MIMARKS.survey.agriculture.6.0","MIMARKS.survey.air.6.0":"MIMARKS.survey.air.6.0","MIMARKS.survey.built.6.0":"MIMARKS.survey.built.6.0","MIMARKS.survey.food-animal.6.0":"MIMARKS.survey.food-animal.6.0","MIMARKS.survey.food-farm.env.6.0":"MIMARKS.survey.food-farm.env.6.0","MIMARKS.survey.food-human.foods.6.0":"MIMARKS.survey.food-human.foods.6.0","MIMARKS.survey.food-prod.facility.6.0":"MIMARKS.survey.food-prod.facility.6.0","MIMARKS.survey.host-associated.6.0":"MIMARKS.survey.host-associated.6.0","MIMARKS.survey.human-associated.6.0":"MIMARKS.survey.human-associated.6.0","MIMARKS.survey.human-gut.6.0":"MIMARKS.survey.human-gut.6.0","MIMARKS.survey.human-oral.6.0":"MIMARKS.survey.human-oral.6.0","MIMARKS.survey.human-skin.6.0":"MIMARKS.survey.human-skin.6.0","MIMARKS.survey.human-vaginal.6.0":"MIMARKS.survey.human-vaginal.6.0","MIMARKS.survey.hydrocarbon-cores.6.0":"MIMARKS.survey.hydrocarbon-cores.6.0","MIMARKS.survey.hydrocarbon-fluids.swabs.6.0":"MIMARKS.survey.hydrocarbon-fluids.swabs.6.0","MIMARKS.survey.microbial.6.0":"MIMARKS.survey.microbial.6.0","MIMARKS.survey.miscellaneous.6.0":"MIMARKS.survey.miscellaneous.6.0","MIMARKS.survey.plant-associated.6.0":"MIMARKS.survey.plant-associated.6.0","MIMARKS.survey.sediment.6.0":"MIMARKS.survey.sediment.6.0","MIMARKS.survey.soil.6.0":"MIMARKS.survey.soil.6.0","MIMARKS.survey.symbiont-associated.6.0":"MIMARKS.survey.symbiont-associated.6.0","MIMARKS.survey.wastewater.6.0":"MIMARKS.survey.wastewater.6.0","MIMARKS.survey.water.6.0":"MIMARKS.survey.water.6.0"},
                        "MIMS.me": {"MIMS.me.agriculture.6.0":"MIMS.me.agriculture.6.0","MIMS.me.air.6.0":"MIMS.me.air.6.0","MIMS.me.built.6.0":"MIMS.me.built.6.0","MIMS.me.food-animal.6.0":"MIMS.me.food-animal.6.0","MIMS.me.food-farm.env.6.0":"MIMS.me.food-farm.env.6.0","MIMS.me.food-human.foods.6.0":"MIMS.me.food-human.foods.6.0","MIMS.me.food-prod.facility.6.0":"MIMS.me.food-prod.facility.6.0","MIMS.me.host-associated.6.0":"MIMS.me.host-associated.6.0","MIMS.me.human-associated.6.0":"MIMS.me.human-associated.6.0","MIMS.me.human-gut.6.0":"MIMS.me.human-gut.6.0","MIMS.me.human-oral.6.0":"MIMS.me.human-oral.6.0","MIMS.me.human-skin.6.0":"MIMS.me.human-skin.6.0","MIMS.me.human-vaginal.6.0":"MIMS.me.human-vaginal.6.0","MIMS.me.hydrocarbon-cores.6.0":"MIMS.me.hydrocarbon-cores.6.0","MIMS.me.hydrocarbon-fluids.swabs.6.0":"MIMS.me.hydrocarbon-fluids.swabs.6.0","MIMS.me.microbial.6.0":"MIMS.me.microbial.6.0","MIMS.me.miscellaneous.6.0":"MIMS.me.miscellaneous.6.0","MIMS.me.plant-associated.6.0":"MIMS.me.plant-associated.6.0","MIMS.me.sediment.6.0":"MIMS.me.sediment.6.0","MIMS.me.soil.6.0":"MIMS.me.soil.6.0","MIMS.me.symbiont-associated.6.0":"MIMS.me.symbiont-associated.6.0","MIMS.me.wastewater.6.0":"MIMS.me.wastewater.6.0","MIMS.me.water.6.0":"MIMS.me.water.6.0"},
                        "MISAG": {"MISAG.6.0":"MISAG.6.0","MISAG.agriculture.6.0":"MISAG.agriculture.6.0","MISAG.air.6.0":"MISAG.air.6.0","MISAG.built.6.0":"MISAG.built.6.0","MISAG.food-animal.6.0":"MISAG.food-animal.6.0","MISAG.food-farm.env.6.0":"MISAG.food-farm.env.6.0","MISAG.food-human.foods.6.0":"MISAG.food-human.foods.6.0","MISAG.food-prod.facility.6.0":"MISAG.food-prod.facility.6.0","MISAG.host-associated.6.0":"MISAG.host-associated.6.0","MISAG.human-associated.6.0":"MISAG.human-associated.6.0","MISAG.human-gut.6.0":"MISAG.human-gut.6.0","MISAG.human-oral.6.0":"MISAG.human-oral.6.0","MISAG.human-skin.6.0":"MISAG.human-skin.6.0","MISAG.human-vaginal.6.0":"MISAG.human-vaginal.6.0","MISAG.hydrocarbon-cores.6.0":"MISAG.hydrocarbon-cores.6.0","MISAG.hydrocarbon-fluids.swabs.6.0":"MISAG.hydrocarbon-fluids.swabs.6.0","MISAG.microbial.6.0":"MISAG.microbial.6.0","MISAG.miscellaneous.6.0":"MISAG.miscellaneous.6.0","MISAG.plant-associated.6.0":"MISAG.plant-associated.6.0","MISAG.sediment.6.0":"MISAG.sediment.6.0","MISAG.soil.6.0":"MISAG.soil.6.0","MISAG.symbiont-associated.6.0":"MISAG.symbiont-associated.6.0","MISAG.wastewater.6.0":"MISAG.wastewater.6.0","MISAG.water.6.0":"MISAG.water.6.0"},
                        "MIUVIG": {"MIUVIG.6.0":"MIUVIG.6.0","MIUVIG.agriculture.6.0":"MIUVIG.agriculture.6.0","MIUVIG.air.6.0":"MIUVIG.air.6.0","MIUVIG.built.6.0":"MIUVIG.built.6.0","MIUVIG.food-animal.6.0":"MIUVIG.food-animal.6.0","MIUVIG.food-farm.env.6.0":"MIUVIG.food-farm.env.6.0","MIUVIG.food-human.foods.6.0":"MIUVIG.food-human.foods.6.0","MIUVIG.food-prod.facility.6.0":"MIUVIG.food-prod.facility.6.0","MIUVIG.host-associated.6.0":"MIUVIG.host-associated.6.0","MIUVIG.human-associated.6.0":"MIUVIG.human-associated.6.0","MIUVIG.human-gut.6.0":"MIUVIG.human-gut.6.0","MIUVIG.human-oral.6.0":"MIUVIG.human-oral.6.0","MIUVIG.human-skin.6.0":"MIUVIG.human-skin.6.0","MIUVIG.human-vaginal.6.0":"MIUVIG.human-vaginal.6.0","MIUVIG.hydrocarbon-cores.6.0":"MIUVIG.hydrocarbon-cores.6.0","MIUVIG.hydrocarbon-fluids.swabs.6.0":"MIUVIG.hydrocarbon-fluids.swabs.6.0","MIUVIG.microbial.6.0":"MIUVIG.microbial.6.0","MIUVIG.miscellaneous.6.0":"MIUVIG.miscellaneous.6.0","MIUVIG.plant-associated.6.0":"MIUVIG.plant-associated.6.0","MIUVIG.sediment.6.0":"MIUVIG.sediment.6.0","MIUVIG.soil.6.0":"MIUVIG.soil.6.0","MIUVIG.symbiont-associated.6.0":"MIUVIG.symbiont-associated.6.0","MIUVIG.wastewater.6.0":"MIUVIG.wastewater.6.0","MIUVIG.water.6.0":"MIUVIG.water.6.0"}
                        },
                    selected="Pathogen.cl.1.0",
                ),
            ),
            # If GenBank checkbox checked then load Genbank Schema Options
            ui.panel_conditional(
                "input.GenBank_checkbox",
                ui.input_select(
                    "GenBank_schemas",
                    label="Select SeqSender GenBank Schema:",
                    choices=["COV", "FLU", "OTHER"],
                    selected="OTHER",
                ),
            ),
            # If GISAID checkbox checked then load GISAID Database Options
            ui.panel_conditional(
                "input.GISAID_checkbox",
                ui.input_select(
                    "GISAID_databases",
                    label="Select GISAID Database:",
                    choices=["FLU", "COV", "POX", "ARBO", "RSV"],
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
                shiny_tools.create_help_tooltip(
                    id = "config_tooltip",
                    description = "Based on previous selections, required fields for the config file can change. The config file is used to generate the necessary files for submission.",
                    position="right",
                ),
            ),
            ui.card(
                ui.page_fluid(
                    div(ui.HTML("Submission:"), style="margin-top:-20px"),
                    ui.panel_conditional(
                        "input.BioSample_checkbox || input.SRA_checkbox || input.GenBank_checkbox",
                        config_indent(1, "NCBI:", custom_style = ""),
                        # Username
                        config_indent(2, "Username:"),
                        config_text_input("ncbi_config_username", placeholder = "NCBI FTP Username", help_msg = ui.p("Username for your NCBI FTP account.", ui.strong("Not your NCBI account username."))),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        # Password
                        config_indent(2, "Password:"),
                        config_text_input("ncbi_config_password", placeholder = "NCBI FTP Password", help_msg = ui.p("Password for your NCBI FTP account.", ui.strong("Not your NCBI account password."))),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        # Spuid Namespace
                        config_indent(2, "Spuid_Namespace:"),
                        config_text_input("ncbi_config_spuid_namespace", placeholder = "NCBI SPUID Namespace", help_msg = ui.p("NCBI SPUID namespace. This should have been provided when you set up your submisison group for NCBI FTP submissions.")),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        # BioSample Package
                        ui.panel_conditional(
                            "input.BioSample_checkbox",
                            config_indent(2, "BioSample_Package:"),
                            div(
                                ui.output_text("BioSample_Package_Name",inline=True,),
                                shiny_tools.create_help_tooltip("config_bs_package", description = ("This is predetermined based on database options selected."), position = "none"),
                                style="display:inline-block;height:5px;"),
                            div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        ),
                        ui.panel_conditional(
                            "input.GenBank_checkbox",
                            # Automatically remove GenBank samples if failed to upload
                            config_indent(2, "GenBank_Auto_Remove_Failed_Samples:"),
                            div(ui.input_radio_buttons(
                                    "ncbi_config_auto_remove_genbank",
                                    label=None,
                                    choices=["True", "False"],
                                    inline=True,
                                ),
                                style="display:inline-block;height:5px;font-size:medium;",
                            ),
                            shiny_tools.create_help_tooltip("ncbi_auto_remove", description = ("If you're submitting to GenBank via FTP, this will automatically remove sequences that fail to pass GenBank's validation and annotation, only if there are samples that did pass. A report file (", ui.code("seq-edit-report.html"), ") will also be generated explaining the reason why the samples were rejected."), position = "none"),
                            div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        ),
                        # Publication title
                        config_indent(2, "Publication_Title:"),
                        config_text_input("ncbi_config_publication_title", placeholder = "Publication Title", help_msg = ui.p("Publication title for your NCBI submissions. If your samples are part of a paper/article then that title should be used here, otherwise, use a descriptive title for the submissions you plan to make with this config file.")),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        # Publication status
                        config_indent(2, "Publication_Status:"),
                        div(ui.input_radio_buttons(
                                "ncbi_config_publication_status",
                                label=None,
                                choices=["Unpublished", "In-press", "Published"],
                                inline=True,
                            ),
                            style="display:inline-block;height:5px;font-size:medium;",
                        ),
                        shiny_tools.create_help_tooltip("ncbi_sub_status", description = ("If a paper is included with your submission then use \"In-press\" or \"Published\". If you haven't or aren't publishing a paper with your samples then use \"Unpublished\"."), position = "none"),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        ui.panel_conditional(
                            "input.GenBank_checkbox && input.GISAID_checkbox",
                            # NCBI submission position
                            config_indent(2, "Submission_Position:"),
                            div(ui.input_radio_buttons(
                                    "ncbi_submission_position",
                                    label=None,
                                    choices={1:"1", 2:"2", "":"None"},
                                    selected=1,
                                    inline=True,
                                ),
                                style="display:inline-block;height:5px;font-size:medium;",
                            ),
                            shiny_tools.create_help_tooltip("ncbi_sub_pos", description = ("If submitting to both GISAID and NCBI, then determine the order you want the accessions linked between databases. If you do not want to link accessions and want GenBank and GISAID submissions to be made at the same time then select ", ui.strong("None"), "."), position = "none"),
                            div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        ),
                        # Specified release date
                        config_indent(2, "Specified_Release_Date:"),
                        config_text_input("ncbi_config_release_date", placeholder = "Release Date", help_msg = "Release date for NCBI submissions. Allows you to specify when should NCBI make your submissions public, allowing you to submit your samples ahead of a planned release. If you want to release samples immediately leave blank or specify a date as \"YYYY-MM-DD\" to release at a specific date. This field also takes the option of a numeric value followed by one of the key words \"days\", \"weeks\", or \"months\", to have SeqSender create a specified release date based on when SeqSender is ran."),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        # Link sample between NCBI databases
                        config_indent(2, "Link_Sample_Between_NCBI_Databases:"),
                        div(ui.input_radio_buttons(
                                "ncbi_config_link_samples",
                                label=None,
                                choices=["True", "False"],
                                inline=True,
                            ),
                            style="display:inline-block;height:5px;font-size:medium;",
                        ),
                        shiny_tools.create_help_tooltip("ncbi_link_sample", description = ("If you're submitting your samples to multiple NCBI databases, do you want your sample's data linked between each database?"), position = "none"),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        # Description category
                        config_indent(2, "Description:", custom_style = ""),
                        # Organization category
                        config_indent(3, "Organization", custom_style = ""),
                        # Organization role
                        config_indent(4, "Role:"),
                        config_text_input("ncbi_config_role", placeholder="Organization Role", help_msg = "Provided by NCBI when setting up your FTP account."),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        # Organization type
                        config_indent(4, "Type:"),
                        config_text_input("ncbi_config_type", placeholder= "Organization Type", help_msg = "Provided by NCBI when setting up your FTP account."),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        # Organization name
                        config_indent(4, "Name:"),
                        config_text_input("ncbi_config_org_name", placeholder = "Organization Name", help_msg = "Submitting group's name."),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        # Address category
                        config_indent(5, "Address:", custom_style = ""),
                        # submitter affiliation
                        config_indent(6, "Affil:"),
                        config_text_input("ncbi_config_affil", placeholder = "Submitter Affiliation", help_msg = "Submitting organization's full name."),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        # submitter division
                        config_indent(6, "Div:"),
                        config_text_input("ncbi_config_div", placeholder = "Submitter Division", help_msg = "Submitting department/laboratory name."),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        # submitter street
                        config_indent(6, "Street:"),
                        config_text_input("ncbi_config_street", placeholder = "Submitting Organization Street Name", help_msg = "Street name."),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        # submitter city
                        config_indent(6, "City:"),
                        config_text_input("ncbi_config_city", placeholder = "Submitting Organization City", help_msg = "City."),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        # submitter state
                        config_indent(6, "Sub:"),
                        config_text_input("ncbi_config_state", placeholder = "Submitting Organization State/Providence", help_msg = "State/Providence."),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        # submitter zipcode
                        config_indent(6, "Postal_Code:"),
                        config_text_input("ncbi_config_postal", placeholder = "Submitting Organization Postal Code", help_msg = "Postal code."),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        # submitter country
                        config_indent(6, "Country:"),
                        config_text_input("ncbi_config_country", placeholder = "Submitting Organization Country", help_msg = "Country"),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        # submitter contact email
                        config_indent(6, "Email:"),
                        config_text_input("ncbi_config_email", placeholder = "Submitting Organization Email", help_msg = "Email for someone to contact you at if they have questions about your samples."),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        # submitter contact phone
                        config_indent(6, "Phone:"),
                        config_text_input("ncbi_config_phone", placeholder = "Submitting Organization Phone", help_msg = "Phone number for someone to contact you at if they have questions about your samples."),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        # submitter point of contact
                        config_indent(5, "Submitter", custom_style = ""),
                        # submitter poc email
                        config_indent(6, "Email:"),
                        config_text_input("ncbi_sub_email_one", "Submitter Contact Email", help_msg = "Email for someone at NCBI to contact about issues during the submission process."),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        # submitter poc email
                        config_indent(6, "Alt_Email:"),
                        config_text_input("ncbi_sub_email_two", "Submitter Contact Email", help_msg = "Alternate email for someone at NCBI to contact about issues during the submission process."),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        # poc name
                        config_indent(6, "Name:", custom_style = ""),
                        # poc name first
                        config_indent(7, "First:"),
                        config_text_input("ncbi_config_first_name", placeholder = "Submitter First Name", help_msg = "Who is making the submission: First Name."),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        config_indent(7, "Last:"),
                        config_text_input("ncbi_config_last_name", placeholder = "Submitter Last Name", help_msg = "Who is making the submission: Last Name."),
                    ),
                    ui.panel_conditional(
                        "input.GISAID_checkbox",
                        config_indent(1, "GISAID", custom_style = ""),
                        config_indent(2, "Client-Id:"),
                        config_text_input("gisaid_config_client", placeholder = "GISAID Client ID", help_msg = "Provided by GISAID when authorized to make FTP submissions."),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        config_indent(2, "Username:"),
                        config_text_input("gisaid_config_username", placeholder = "GISAID Username", help_msg = "Username to login to GISAID."),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        config_indent(2, "Password:"),
                        config_text_input("gisaid_config_password", placeholder = "GISAID Password", help_msg = "Password to login to GISAID."),
                        div(ui.HTML("<br>"), style="margin-top:-25px;"),
                        ui.panel_conditional(
                            "input.GISAID_checkbox && input.GenBank_checkbox",
                            config_indent(2, "Submission_Position:"),
                            div(ui.input_radio_buttons(
                                    "gisaid_submission_position",
                                    label=None,
                                    choices={1:"1", 2:"2", "":"None"},
                                    inline=True,
                                    selected=2,
                                ),
                                style="display:inline-block;height:5px;font-size:medium;",
                            ),
                            shiny_tools.create_help_tooltip("gisaid_sub_pos", description = ("If submitting to both GISAID and NCBI, then determine the order you want the accessions linked between databases. If you do not want to link accessions and want GenBank and GISAID submissions to be made at the same time then select ", ui.strong("None"), "."), position = "none"),
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
