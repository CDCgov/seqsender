
# Python Libraries
import os
import sys
import importlib
import pathlib
import pandas as pd
# Local imports
sys.path.insert(0, str(pathlib.Path(__file__).parent))

# Get program directory
PROG_DIR: str = os.path.dirname(os.path.abspath(__file__))

# When working with all schemas exclude these specific unique ones
SCHEMA_EXCLUSIONS = ["config.seqsender_upload_log_schema","config.seqsender_upload_log_schema","config.config_file.gisaid_schema","config.config_file.ncbi_gisaid_schema", "config.config_file.ncbi_schema"]

# Collect all schema files from config subdirectory
# Convert file name to importable string
def get_all_schema_files():
    # Load all schema files
    files = os.walk(os.path.join(PROG_DIR, "config"))
    schemas = []
    for root, dir_list, file_list in files:
        # Get seqsender subdirectory structure
        subdirectory = root.split("/seqsender/config")[-1].replace("/", ".")
        root_end = "config."
        if subdirectory and subdirectory != ".":
            root_end = "config" + subdirectory + "."
        # Collect schemas
        schemas += [(root_end + file.replace(".py", "")) for file in file_list if file.endswith(".py")]
    return schemas

# Import pandera schema as module
def load_schema(schema_name: str):
    return importlib.import_module(schema_name).schema

# Process schema file to retrieve field information
def process_schema(schema):
    schema_contents = []
    for column_name, context in schema.columns.items():
        required_field = context.required
        description_field = context.description
        # Convert required boolean to string
        if required_field:
            required_field = "Required"
        else:
            required_field = "Optional"
        # Update required field if required group of columns
        if description_field and "At least one required: Group" in description_field:
            required_field = "At least one field required. Group: " + description_field.split("Group: \"")[-1].split("\".")[0]
        schema_contents.append({"column_name": column_name, "required_column": required_field, "description": description_field})
    return pd.DataFrame(schema_contents)

# Update all shiny templates based on all schema files in config folder
def update_all_schema_templates():
    schemas_list = get_all_schema_files()
    for schema_name in schemas_list:
        # Skip schemas in exclusion list
        if schema_name in SCHEMA_EXCLUSIONS:
            continue
        try:
            schema = load_schema(schema_name)
        except Exception as e:
            print("Warning: Unable to load schema \"" + schema_name + "\".", file=sys.stderr)
            print(e, file=sys.stderr)
            continue
        template = dict()
        try:
            metadata_template = process_schema(schema)
        except:
            print("Warning: Unable to process schema into metadata template.", file=sys.stderr)
            print(e, file=sys.stderr)
            continue
        metadata_template.to_csv(os.path.join(PROG_DIR, "shiny", "templates", schema_name.replace("_", ".") + "_template.csv"), header = True, index = False)

if __name__ == "__main__":
	update_all_schema_templates()
