# tests/unit/test_biosample_sra_handler.py
from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path
from typing import Any
import os
import pandas as pd
import pytest
from lxml import etree

def _source_root() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        # During mutmut, prefer the mutated source tree.
        mutant_src = parent / "src"
        if parent.name == "mutants" and (mutant_src / "biosample_sra_handler.py").exists():
            return mutant_src

        # Normal pytest run.
        normal_src = parent / "src"
        if (normal_src / "biosample_sra_handler.py").exists() and parent.name != "mutants":
            return normal_src

    raise RuntimeError("Could not find src/biosample_sra_handler.py")


SOURCE_DIR = _source_root()
MODULE_PATH = SOURCE_DIR / "biosample_sra_handler.py"

if str(SOURCE_DIR) not in sys.path:
    sys.path.insert(0, str(SOURCE_DIR))

"""Unit tests initially generated with ChatGPT v5.5 "Deep Research" and "Thinking"
then modified for use and clarity with validation and coverage testing via mutmut."""

#*******************************************************************************
#              create test biosample_sra_handler.py connections
#*******************************************************************************

@pytest.fixture()
def handler(monkeypatch: pytest.MonkeyPatch):
    src_pkg: Any = types.ModuleType("src")
    src_pkg.__path__ = [str(SOURCE_DIR)]
    monkeypatch.setitem(sys.modules, "src", src_pkg)

    def alias_src_module(name, module):
        monkeypatch.setitem(sys.modules, f"src.{name}", module)
        monkeypatch.setitem(sys.modules, name, module)
        setattr(src_pkg, name, module)

    settings_stub: Any = types.ModuleType("settings")
    settings_stub.BIOSAMPLE_REGEX = "^bs-|^bioproject$|^organism$|^collection_date$"
    settings_stub.SRA_REGEX = "^sra-|^bioproject$|bs-sample_name|^organism$|^collection_date$"

    ncbi_stub: Any = types.ModuleType("ncbi_handler")

    def process_report_header(report_file):
        return({"SubmissionStatus": {}}, "SUBMITTED", "SUB123")

    ncbi_stub.process_report_header = process_report_header

    file_handler_stub: Any = types.ModuleType("file_handler")

    def save_csv(**kwargs):
        return None

    def save_xml(*args, **kwargs):
        return None

    file_handler_stub.save_csv = save_csv
    file_handler_stub.save_xml = save_xml

    upload_log_stub: Any = types.ModuleType("upload_log")

    def update_submission_status_csv(**kwargs):
        return None

    upload_log_stub.update_submission_status_csv = update_submission_status_csv

    alias_src_module("ncbi_handler", ncbi_stub)
    alias_src_module("file_handler", file_handler_stub)
    alias_src_module("upload_log", upload_log_stub)
    alias_src_module("settings", settings_stub)

    sys.modules.pop("biosample_sra_handler", None)
    sys.modules.pop("src.biosample_sra_handler", None)

    spec = importlib.util.spec_from_file_location("biosample_sra_handler", MODULE_PATH)
    assert spec and spec.loader

    module: Any = importlib.util.module_from_spec(spec)
    sys.modules["biosample_sra_handler"] = module
    setattr(src_pkg, "biosample_sra_handler", module)
    spec.loader.exec_module(module)
    return module

#*******************************************************************************
#                       create fake ncbi config file
#*******************************************************************************

@pytest.fixture()
def config_dict() -> dict[str, Any]:
    return {
        "Spuid_Namespace": "CDC-NS",
        "BioSample_Package": "Pathogen.cl.1.0",
        "Specified_Release_Date": "2027-01-15",
        "Description": {
            "Organization": {
                "Type": "center",
                "Role": "owner",
                "Name": "Example Lab",
                "Submitter": {
                    "Email": "submitter@example.org",
                    "Name": {"First": "Ada", "Last": "Lovelace"},
                },
            }
        },
    }

#*******************************************************************************
#                      create fake metadata file
#*******************************************************************************

@pytest.fixture()
def biosample_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "organism": ["Influenza A virus"],
            "collection_date": ["2025-01-02"],
            "bioproject": ["PRJNA000001"],
            "bs-sample_name": ["BS1"],
            "bs-sample_title": ["BioSample title"],
            "bs-sample_description": ["BioSample description"],
            "bs-strain": ["strain-1"],
            "bs-host": ["Homo sapiens"],
            "bs-title": ["Portal BioSample title"],
            "bs-comment": ["Portal BioSample comment"],
        }
    )


@pytest.fixture()
def sra_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "organism": ["Influenza A virus"],
            "collection_date": ["2025-01-02"],
            "bioproject": ["PRJNA000001"],
            "bs-sample_name": ["BS1"],
            "sra-sample_name": ["SRA1"],
            "sra-title": ["Portal SRA title"],
            "sra-comment": ["Portal SRA comment"],
            "sra-file_location": ["local"],
            "sra-file_1": ["reads_R1.fastq.gz"],
            "sra-file_2": [""],
            "sra-library_name": ["lib-1"],
            "sra-platform": ["ILLUMINA"],
            "sra-loader": ["ignored-loader"],
            "sra-file": ["legacy-file-column"],
        }
    )


def parse_xml(xml_bytes: bytes) -> etree._Element:
    return etree.fromstring(xml_bytes)

#*******************************************************************************
#                          check_raw_read_files
#*******************************************************************************

def test_check_raw_read_files__returns_relative_and_absolute_local_files(handler, tmp_path):
    raw_reads = tmp_path / "raw_reads"
    raw_reads.mkdir()
    rel_file_1 = raw_reads / "sample_R1.fastq.gz"
    rel_file_2 = raw_reads / "sample_R2.fastq.gz"
    rel_file_1.write_text("reads")
    rel_file_2.write_text("reads")

    metadata = pd.DataFrame(
        {
            "sra-sample_name": ["SRA1", "SRA2", "SRA3"],
            "sra-file_location": ["local", "local", "cloud"],
            "sra-file_1": ["sample_R1.fastq.gz", str(rel_file_1), "s3://bucket/ignored.fastq.gz"],
            "sra-file_2": ["sample_R2.fastq.gz", str(rel_file_2), "s3://bucket/ignored.fastq.gz"],
        }
    )

    observed = handler.check_raw_read_files("submission", str(tmp_path), metadata)
    assert observed == {
        str(raw_reads / "sample_R1.fastq.gz"),
        str(raw_reads / "sample_R2.fastq.gz"),
        str(rel_file_1),
        str(rel_file_2),
    }
    assert os.path.basename(raw_reads) == "raw_reads"
    assert (raw_reads / "sample_R1.fastq.gz").exists()
    assert (raw_reads / "sample_R2.fastq.gz").exists()

def test_check_raw_read_files__continues_after_blank_extra_file_to_validate_later_extra_file(handler, tmp_path):
    raw_reads = tmp_path / "raw_reads"
    raw_reads.mkdir()
    (raw_reads / "sample_R1.fastq.gz").write_text("reads")
    (raw_reads / "sample_R3.fastq.gz").write_text("reads")
    metadata = pd.DataFrame(
        {
            "sra-sample_name": ["SRA1"],
            "sra-file_location": ["local"],
            "sra-file_1": ["sample_R1.fastq.gz"],
            "sra-file_2": ["   "],
            "sra-file_3": ["sample_R3.fastq.gz"],
        }
    )
    assert handler.check_raw_read_files("submission", str(tmp_path), metadata) == {str(raw_reads / "sample_R1.fastq.gz"), str(raw_reads / "sample_R3.fastq.gz")}

def test_check_raw_read_files__relative_path_must_join_raw_reads_directory(handler, tmp_path):
    raw_reads = tmp_path / "raw_reads"
    raw_reads.mkdir()
    nested = raw_reads / "nested"
    nested.mkdir()
    expected_file = nested / "reads.fastq.gz"
    expected_file.write_text("reads")

    metadata = pd.DataFrame(
        {
            "sra-sample_name": ["SRA1"],
            "sra-file_location": ["local"],
            "sra-file_1": ["nested/reads.fastq.gz"],
        }
    )

    observed = handler.check_raw_read_files("submission", str(tmp_path), metadata)
    assert observed == {str(expected_file)}

def test_check_raw_read_files__only_lowercase_local_rows_are_checked(handler, tmp_path):
    raw_reads = tmp_path / "raw_reads"
    raw_reads.mkdir()
    (raw_reads / "lowercase.fastq.gz").write_text("reads")

    metadata = pd.DataFrame(
        {
            "sra-sample_name": ["LOWER", "UPPER", "TITLE", "CLOUD"],
            "sra-file_location": ["local", "LOCAL", "Local", "cloud"],
            "sra-file_1": [
                "lowercase.fastq.gz",
                "missing_upper.fastq.gz",
                "missing_title.fastq.gz",
                "missing_cloud.fastq.gz",
            ],
        }
    )

    observed = handler.check_raw_read_files("submission", str(tmp_path), metadata)
    assert observed == {str(raw_reads / "lowercase.fastq.gz")}

def test_check_raw_read_files__blank_required_first_file_exits_exactly(handler, tmp_path, capsys):
    metadata = pd.DataFrame(
        {
            "sra-sample_name": ["SRA1"],
            "sra-file_location": ["local"],
            "sra-file_1": ["   "],
            "sra-file_2": ["also_ignored.fastq.gz"],
        }
    )

    with pytest.raises(SystemExit) as exc:
        handler.check_raw_read_files("submission", str(tmp_path), metadata)

    assert exc.value.code == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == (
        "Error: Raw read files for SRA1 could not be found. Field 'sra-file_1' must either be the full file path, or if just the file name it must be stored at '<submission_dir>/<submission_name>/raw_reads/<sra-file>'.\n"
        "Error: Raw read files for SRA1 could not be found. Field 'sra-file_2' must either be the full file path, or if just the file name it must be stored at '<submission_dir>/<submission_name>/raw_reads/<sra-file>'.\n"
    )

def test_check_raw_read_files__exits_when_required_local_file_missing(handler, tmp_path, capsys):
    metadata = pd.DataFrame(
        {
            "sra-sample_name": ["SRA1"],
            "sra-file_location": ["local"],
            "sra-file_1": ["missing.fastq.gz"],
        }
    )

    with pytest.raises(SystemExit) as exc:
        handler.check_raw_read_files("submission", str(tmp_path), metadata)

    assert exc.value.code == 1
    assert "Error: Raw read files for SRA1 could not be found. Field 'sra-file_1' must either be the full file path, or if just the file name it must be stored at '<submission_dir>/<submission_name>/raw_reads/<sra-file>'.\n" == capsys.readouterr().err


def test_check_raw_read_files__ignores_blank_string_extra_file_columns(handler, tmp_path):
    raw_reads = tmp_path / "raw_reads"
    raw_reads.mkdir()
    (raw_reads / "sample_R1.fastq.gz").write_text("reads")
    metadata = pd.DataFrame(
        {
            "sra-sample_name": ["SRA1"],
            "sra-file_location": ["local"],
            "sra-file_1": ["sample_R1.fastq.gz"],
            "sra-file_2": ["   "],
            "sra-file_3": [None],
        }
    )

    observed = handler.check_raw_read_files("submission", str(tmp_path), metadata)

    assert observed == {str(raw_reads / "sample_R1.fastq.gz")}

#*******************************************************************************
#                      create_manual_submission_files
#*******************************************************************************

def test_create_manual_submission_files__sra_filters_renames_drops_prefixes_and_orders(handler, monkeypatch, sra_metadata, tmp_path):
    calls: list[dict[str, Any]] = []

    def capture_save_csv(**kwargs):
        calls.append(kwargs)

    monkeypatch.setattr(handler.file_handler, "save_csv", capture_save_csv)

    handler.create_manual_submission_files("SRA", str(tmp_path), sra_metadata, config_dict={})

    assert len(calls) == 1
    saved_df = calls[0]["df"]
    assert calls[0]["file_path"] == str(tmp_path)
    assert calls[0]["file_name"] == "metadata.tsv"
    assert calls[0]["sep"] == "\t"
    # Current implementation keeps portal-only sra-title/sra-comment in the
    # manual SRA TSV; this assertion documents that behavior.
    assert saved_df.columns.tolist() == [
        "sample_name",
        "library_ID",
        "organism",
        "collection_date",
        "title",
        "comment",
        "filename",
        "filename2",
        "platform",
    ]
    assert saved_df.to_dict("records") == [
        {
            "sample_name": "SRA1",
            "library_ID": "lib-1",
            "organism": "Influenza A virus",
            "collection_date": "2025-01-02",
            "title": "Portal SRA title",
            "comment": "Portal SRA comment",
            "filename": "reads_R1.fastq.gz",
            "filename2": "",
            "platform": "ILLUMINA",
        }
    ]
    assert "file_location" not in saved_df.columns
    assert "loader" not in saved_df.columns


def test_create_manual_submission_files__biosample_filters_renames_drops_prefixes_and_orders(handler, monkeypatch, biosample_metadata, tmp_path):
    calls: list[dict[str, Any]] = []
    monkeypatch.setattr(handler.file_handler, "save_csv", lambda **kwargs: calls.append(kwargs))

    handler.create_manual_submission_files("BIOSAMPLE", str(tmp_path), biosample_metadata, config_dict={})

    saved_df = calls[0]["df"]
    assert saved_df.columns.tolist() == [
        "sample_name",
        "organism",
        "collection_date",
        "strain",
        "host",
    ]
    assert saved_df.to_dict("records") == [
        {
            "sample_name": "BS1",
            "organism": "Influenza A virus",
            "collection_date": "2025-01-02",
            "strain": "strain-1",
            "host": "Homo sapiens",
        }
    ]
    assert "title" not in saved_df.columns
    assert "comment" not in saved_df.columns


def test_create_manual_submission_files__invalid_database_exits(handler, tmp_path, capsys):
    with pytest.raises(SystemExit) as exc:
        handler.create_manual_submission_files("GENBANK", str(tmp_path), pd.DataFrame(), config_dict={})

    assert exc.value.code == 1
    assert "Error: create_manual_submission_files function only for databases SRA/BioSample. Not '{database}'.\n" == capsys.readouterr().err

#*******************************************************************************
#                      create_submission_xml BioSample
#*******************************************************************************

def test_create_submission_xml__biosample_custom_title_comment_hold_and_attributes(handler, config_dict, biosample_metadata):
    root = parse_xml(handler.create_submission_xml(organism="FLU", database="BIOSAMPLE", submission_name="SUB1", config_dict=config_dict, metadata=biosample_metadata))

    assert root.findtext("Description/Title") == "Portal BioSample title"
    assert root.findtext("Description/Comment") == "Portal BioSample comment"
    assert root.find("Description/Hold").get("release_date") == "2027-01-15"
    assert root.find("Description/Organization").get("type") == "center"
    assert root.findtext("Description/Organization/Contact/Name/First") == "Ada"

    biosample = root.find(".//BioSample")
    assert biosample is not None
    assert biosample.findtext("SampleId/SPUID") == "BS1"
    assert biosample.find("SampleId/SPUID").get("spuid_namespace") == "CDC-NS"
    assert biosample.findtext("Descriptor/Title") == "BioSample title"
    assert biosample.findtext("Descriptor/Description") == "BioSample description"
    assert biosample.findtext("Organism/OrganismName") == "Influenza A virus"
    assert biosample.findtext("BioProject/PrimaryId") == "PRJNA000001"
    assert biosample.findtext("Package") == "Pathogen.cl.1.0"

    attrs = {node.get("attribute_name"): node.text for node in biosample.findall("Attributes/Attribute")}
    assert attrs["strain"] == "strain-1"
    assert attrs["host"] == "Homo sapiens"
    assert attrs["collection_date"] == "2025-01-02"

    identifier = root.find(".//AddData/Identifier/SPUID")
    assert identifier.text == "BS1"
    assert identifier.get("spuid_namespace") == "CDC-NS"

def test_create_submission_xml__biosample_defaults_title_comment_when_portal_fields_absent(handler, config_dict, biosample_metadata):
    config_dict["Specified_Release_Date"] = ""
    metadata = biosample_metadata.drop(columns=["bs-title", "bs-comment"])
    root = parse_xml(handler.create_submission_xml("FLU", "BIOSAMPLE", "SUB1", config_dict, metadata))
    assert root.findtext("Description/Title") == "SUB1-BS"
    assert root.findtext("Description/Comment") == "BioSample Submission"
    assert root.find("Description/Hold") is None

def test_create_submission_xml__biosample_blank_portal_title_comment_are_preserved_current_behavior(handler, config_dict, biosample_metadata):
    metadata = biosample_metadata.copy()
    metadata["bs-title"] = "   "
    metadata["bs-comment"] = "   "
    root = parse_xml(handler.create_submission_xml(organism="FLU", database="BIOSAMPLE", submission_name="SUB1", config_dict=config_dict, metadata=metadata))
    assert root.findtext("Description/Title") == "   "
    assert root.findtext("Description/Comment") == "   "

def test_create_submission_xml__serializes_with_exact_tostring_options(handler, monkeypatch, config_dict, biosample_metadata):
    original_tostring = handler.etree.tostring
    observed: list[dict[str, Any]] = []

    def fake_tostring(element, *args, **kwargs):
        observed.append({"tag": element.tag, "args": args, "kwargs": kwargs})
        return original_tostring(element, *args, **kwargs)

    monkeypatch.setattr(handler.etree, "tostring", fake_tostring)
    handler.create_submission_xml(organism="FLU", database="BIOSAMPLE", submission_name="SUB1", config_dict=config_dict, metadata=biosample_metadata)
    assert observed == [
        {
            "tag": "Submission",
            "args": (),
            "kwargs": {
                "encoding": "utf-8",
                "pretty_print": True,
                "xml_declaration": True,
            },
        }
    ]

def test_create_submission_xml__biosample_has_exact_required_xml_structure(handler, config_dict, biosample_metadata):
    xml = handler.create_submission_xml(organism="FLU", database="BIOSAMPLE", submission_name="SUB1", config_dict=config_dict, metadata=biosample_metadata)
    assert xml.startswith(b"<?xml version='1.0' encoding='utf-8'?>\n")
    root = parse_xml(xml)
    assert root.tag == "Submission"

    description = root.find("Description")
    assert description is not None
    assert description.findtext("Title") == "Portal BioSample title"
    assert description.findtext("Comment") == "Portal BioSample comment"
    assert description.find("Hold").attrib == {"release_date": "2027-01-15"}

    organization = description.find("Organization")
    assert organization is not None
    assert organization.attrib == {"type": "center", "role": "owner"}
    assert organization.findtext("Name") == "Example Lab"

    contact = organization.find("Contact")
    assert contact is not None
    assert contact.attrib == {"email": "submitter@example.org"}
    assert contact.findtext("Name/First") == "Ada"
    assert contact.findtext("Name/Last") == "Lovelace"

    add_data = root.find("Action/AddData")
    assert add_data is not None
    assert add_data.attrib == {"target_db": "BioSample"}

    data = add_data.find("Data")
    assert data is not None
    assert data.attrib == {"content_type": "xml"}

    biosample = add_data.find("Data/XmlContent/BioSample")
    assert biosample is not None
    assert biosample.attrib == {"schema_version": "2.0"}
    assert biosample.find("SampleId/SPUID").attrib == {"spuid_namespace": "CDC-NS"}
    assert biosample.findtext("SampleId/SPUID") == "BS1"
    assert biosample.findtext("Descriptor/Title") == "BioSample title"
    assert biosample.findtext("Descriptor/Description") == "BioSample description"
    assert biosample.findtext("Organism/OrganismName") == "Influenza A virus"

    primary_id = biosample.find("BioProject/PrimaryId")
    assert primary_id is not None
    assert primary_id.attrib == {"db": "BioProject"}
    assert primary_id.text == "PRJNA000001"
    assert biosample.findtext("Package") == "Pathogen.cl.1.0"
    attributes = [
        (node.attrib, node.text)
        for node in biosample.findall("Attributes/Attribute")
    ]
    assert attributes == [
        ({"attribute_name": "strain"}, "strain-1"),
        ({"attribute_name": "host"}, "Homo sapiens"),
        ({"attribute_name": "collection_date"}, "2025-01-02"),
    ]
    identifier = add_data.find("Identifier/SPUID")
    assert identifier is not None
    assert identifier.attrib == {"spuid_namespace": "CDC-NS"}
    assert identifier.text == "BS1"

def test_create_submission_xml__biosample_descriptor_with_blank_description_omits_description_node(handler, config_dict, biosample_metadata):
    metadata = biosample_metadata.copy()
    metadata["bs-sample_title"] = "Only title"
    metadata["bs-sample_description"] = ""
    root = parse_xml(handler.create_submission_xml(organism="FLU", database="BIOSAMPLE", submission_name="SUB1", config_dict=config_dict, metadata=metadata))
    descriptor = root.find(".//BioSample/Descriptor")
    assert descriptor is not None
    assert descriptor.findtext("Title") == "Only title"
    assert descriptor.find("Description") is None

def test_create_submission_xml__biosample_descriptor_with_blank_title_and_description_omits_descriptor(handler, config_dict, biosample_metadata):
    metadata = biosample_metadata.copy()
    metadata["bs-sample_title"] = ""
    metadata["bs-sample_description"] = ""
    root = parse_xml(handler.create_submission_xml(organism="FLU", database="BIOSAMPLE", submission_name="SUB1", config_dict=config_dict, metadata=metadata))
    assert root.find(".//BioSample/Descriptor") is None

def test_create_submission_xml__biosample_skips_blank_optional_attributes_but_keeps_collection_date(handler, config_dict, biosample_metadata):
    metadata = biosample_metadata.copy()
    metadata["bs-strain"] = ""
    metadata["bs-host"] = "Homo sapiens"
    metadata["bs-package"] = "must-not-be-an-attribute"
    root = parse_xml(handler.create_submission_xml(organism="FLU", database="BIOSAMPLE", submission_name="SUB1", config_dict=config_dict, metadata=metadata))
    attributes = [(node.attrib, node.text) for node in root.findall(".//BioSample/Attributes/Attribute")]
    assert attributes == [({"attribute_name": "host"}, "Homo sapiens"), ({"attribute_name": "collection_date"}, "2025-01-02")]

def test_create_submission_xml__biosample_blank_bioproject_omits_bioproject_node(handler, config_dict, biosample_metadata):
    metadata = biosample_metadata.copy()
    metadata["bioproject"] = ""
    root = parse_xml(handler.create_submission_xml(organism="FLU", database="BIOSAMPLE", submission_name="SUB1", config_dict=config_dict, metadata=metadata))
    assert root.find(".//BioSample/BioProject") is None

def test_create_submission_xml__biosample_without_bioproject_column_omits_bioproject_node(handler, config_dict, biosample_metadata):
    metadata = biosample_metadata.drop(columns=["bioproject"])
    root = parse_xml(handler.create_submission_xml(organism="FLU", database="BIOSAMPLE", submission_name="SUB1", config_dict=config_dict, metadata=metadata))
    assert root.find(".//BioSample/BioProject") is None

def test_create_submission_xml__biosample_descriptor_absent_when_title_and_description_columns_absent(handler, config_dict, biosample_metadata):
    metadata = biosample_metadata.drop(columns=["bs-sample_title", "bs-sample_description"])
    root = parse_xml(handler.create_submission_xml("FLU", "BIOSAMPLE", "SUB1", config_dict, metadata))
    assert root.find(".//BioSample/Descriptor") is None

def test_create_submission_xml__biosample_descriptor_created_when_only_title_has_value(handler, config_dict, biosample_metadata):
    metadata = biosample_metadata.copy()
    metadata["bs-sample_title"] = "Only title"
    metadata["bs-sample_description"] = pd.NA
    root = parse_xml(handler.create_submission_xml("FLU", "BIOSAMPLE", "SUB1", config_dict, metadata))
    descriptor = root.find(".//BioSample/Descriptor")
    assert descriptor is not None
    assert descriptor.findtext("Title") == "Only title"
    assert descriptor.find("Description") is None

def test_create_submission_xml__biosample_descriptor_created_when_only_description_has_value(handler, config_dict, biosample_metadata):
    metadata = biosample_metadata.copy()
    metadata["bs-sample_title"] = pd.NA
    metadata["bs-sample_description"] = "Only description"
    root = parse_xml(handler.create_submission_xml("FLU", "BIOSAMPLE", "SUB1", config_dict, metadata))
    descriptor = root.find(".//BioSample/Descriptor")
    assert descriptor is not None
    assert descriptor.find("Title") is None
    assert descriptor.findtext("Description") == "Only description"

#*******************************************************************************
#                        create_submission_xml SRA
#*******************************************************************************

def test_create_submission_xml__sra_has_exact_required_xml_structure(handler, config_dict, sra_metadata):
    xml = handler.create_submission_xml(organism="FLU", database="SRA", submission_name="SUB1", config_dict=config_dict, metadata=sra_metadata)

    assert xml.startswith(b"<?xml version='1.0' encoding='utf-8'?>\n")

    root = parse_xml(xml)
    assert root.tag == "Submission"

    description = root.find("Description")
    assert description is not None
    assert description.findtext("Title") == "Portal SRA title"
    assert description.findtext("Comment") == "Portal SRA comment"

    organization = description.find("Organization")
    assert organization is not None
    assert organization.attrib == {"type": "center", "role": "owner"}
    assert organization.findtext("Name") == "Example Lab"

    contact = organization.find("Contact")
    assert contact is not None
    assert contact.attrib == {"email": "submitter@example.org"}
    assert contact.findtext("Name/First") == "Ada"
    assert contact.findtext("Name/Last") == "Lovelace"

    addfiles = root.find("Action/AddFiles")
    assert addfiles is not None
    assert addfiles.attrib == {"target_db": "SRA"}

    file_node = addfiles.find("File")
    assert file_node is not None
    assert file_node.attrib == {"file_path": "reads_R1.fastq.gz"}
    assert file_node.findtext("DataType") == "generic-data"
    assert [
        (node.attrib, node.text)
        for node in addfiles.findall("Attribute")
    ] == [
        ({"name": "library_name"}, "lib-1"),
        ({"name": "platform"}, "ILLUMINA"),
        ({"name": "loader"}, "ignored-loader"),
        ({"name": "file"}, "legacy-file-column"),
    ]

    refids = addfiles.findall("AttributeRefId")
    assert [(node.attrib, node.findtext("RefId/PrimaryId"), node.findtext("RefId/SPUID")) for node in refids] == [
        ({"name": "BioProject"}, "PRJNA000001", None),
        ({"name": "BioSample"}, None, "BS1"),
    ]
    assert refids[1].find("RefId/SPUID").attrib == {"spuid_namespace": "CDC-NS"}

    identifier = addfiles.find("Identifier/SPUID")
    assert identifier is not None
    assert identifier.attrib == {"spuid_namespace": "CDC-NS"}
    assert identifier.text == "SRA1"

def test_create_submission_xml__sra_local_file_attributes_and_refs(handler, config_dict, sra_metadata):
    root = parse_xml(handler.create_submission_xml("FLU", "SRA", "SUB1", config_dict, sra_metadata))
    assert root.findtext("Description/Title") == "Portal SRA title"
    assert root.findtext("Description/Comment") == "Portal SRA comment"
    addfiles = root.find(".//AddFiles")
    assert addfiles.get("target_db") == "SRA"
    file_node = addfiles.find("File")
    assert file_node.get("file_path") == "reads_R1.fastq.gz"
    assert file_node.findtext("DataType") == "generic-data"

    attrs = {node.get("name"): node.text for node in addfiles.findall("Attribute")}
    assert attrs["library_name"] == "lib-1"
    assert attrs["platform"] == "ILLUMINA"

    refids = addfiles.findall("AttributeRefId")
    assert [node.get("name") for node in refids] == ["BioProject", "BioSample"]
    assert refids[0].findtext("RefId/PrimaryId") == "PRJNA000001"
    assert refids[1].findtext("RefId/SPUID") == "BS1"
    assert addfiles.findtext("Identifier/SPUID") == "SRA1"

def test_create_submission_xml__sra_cloud_file_has_exact_file_node(handler, config_dict, sra_metadata):
    metadata = sra_metadata.copy()
    metadata["sra-file_location"] = "cloud"
    metadata["sra-file_1"] = "s3://bucket/reads_R1.fastq.gz"
    root = parse_xml(handler.create_submission_xml("FLU", "SRA", "SUB1", config_dict, metadata))
    file_node = root.find(".//AddFiles/File")
    assert file_node is not None
    assert file_node.attrib == {"cloud_url": "s3://bucket/reads_R1.fastq.gz"}
    assert file_node.findtext("DataType") == "generic-data"

def test_create_submission_xml__sra_cloud_file_uses_exact_cloud_file_node(handler, config_dict, sra_metadata):
    metadata = sra_metadata.copy()
    metadata["sra-file_location"] = "cloud"
    metadata["sra-file_1"] = "s3://bucket/reads_R1.fastq.gz"
    root = parse_xml(handler.create_submission_xml("FLU", "SRA", "SUB1", config_dict, metadata))
    file_node = root.find(".//AddFiles/File")
    assert file_node is not None
    assert file_node.attrib == {"cloud_url": "s3://bucket/reads_R1.fastq.gz"}
    assert file_node.findtext("DataType") == "generic-data"

def test_create_submission_xml__sra_blank_portal_title_comment_are_preserved_current_behavior(handler, config_dict, sra_metadata):
    metadata = sra_metadata.copy()
    metadata["sra-title"] = "   "
    metadata["sra-comment"] = "   "
    root = parse_xml(handler.create_submission_xml(organism="FLU", database="SRA", submission_name="SUB1", config_dict=config_dict, metadata=metadata))
    assert root.findtext("Description/Title") == "   "
    assert root.findtext("Description/Comment") == "   "

def test_create_submission_xml__sra_cloud_multiple_files_use_exact_file_tags(handler, config_dict, sra_metadata):
    metadata = sra_metadata.copy()
    metadata["sra-file_location"] = "cloud"
    metadata["sra-file_1"] = "s3://bucket/reads_R1.fastq.gz"
    metadata["sra-file_2"] = "s3://bucket/reads_R2.fastq.gz"
    root = parse_xml(handler.create_submission_xml("FLU", "SRA", "SUB1", config_dict, metadata))
    file_nodes = root.findall(".//AddFiles/File")
    assert [node.tag for node in file_nodes] == ["File", "File"]
    assert [node.attrib for node in file_nodes] == [{"cloud_url": "s3://bucket/reads_R1.fastq.gz"}, {"cloud_url": "s3://bucket/reads_R2.fastq.gz"}]

def test_create_submission_xml__sra_skips_blank_second_file_and_keeps_third_file(handler, config_dict, sra_metadata):
    metadata = sra_metadata.copy()
    metadata["sra-file_1"] = "reads_R1.fastq.gz"
    metadata["sra-file_2"] = "   "
    metadata["sra-file_3"] = "reads_R3.fastq.gz"
    root = parse_xml(handler.create_submission_xml(organism="FLU", database="SRA", submission_name="SUB1", config_dict=config_dict, metadata=metadata))
    file_nodes = root.findall(".//AddFiles/File")
    assert [node.attrib for node in file_nodes] == [{"file_path": "reads_R1.fastq.gz"}, {"file_path": "reads_R3.fastq.gz"}]

def test_create_submission_xml__sra_defaults_title_comment_when_fields_absent(handler, config_dict, sra_metadata):
    metadata = sra_metadata.drop(columns=["sra-title", "sra-comment"])
    root = parse_xml(handler.create_submission_xml("FLU", "SRA", "SUB1", config_dict, metadata))
    assert root.findtext("Description/Title") == "SUB1-SRA"
    assert root.findtext("Description/Comment") == "SRA Submission"

@pytest.mark.parametrize(
    ("updates", "expected_message"),
    [
        (
            {"sra-file_1": ""},
            "Error: metadata must contain a file for SRA1 in column sra-file_1\n",
        ),
        (
            {"sra-file_location": "ftp"},
            "Error: Metadata field file_location must be either cloud or local. Field currently contains: ftp\n",
        ),
    ],
)
def test_create_submission_xml__sra_invalid_file_fields_exit(handler, config_dict, sra_metadata, updates, expected_message, capsys):
    metadata = sra_metadata.copy()
    for col, value in updates.items():
        metadata[col] = value

    with pytest.raises(SystemExit) as exc:
        handler.create_submission_xml("FLU", "SRA", "SUB1", config_dict, metadata)

    assert exc.value.code == 1
    assert expected_message == capsys.readouterr().err

def test_create_submission_xml__sra_keeps_second_file_and_skips_blank_third_file(handler, config_dict, sra_metadata):
    metadata = sra_metadata.copy()
    metadata["sra-file_1"] = "reads_R1.fastq.gz"
    metadata["sra-file_2"] = "reads_R2.fastq.gz"
    metadata["sra-file_3"] = "   "
    root = parse_xml(handler.create_submission_xml(organism="FLU", database="SRA", submission_name="SUB1", config_dict=config_dict, metadata=metadata))
    file_nodes = root.findall(".//AddFiles/File")
    assert [node.attrib for node in file_nodes] == [
        {"file_path": "reads_R1.fastq.gz"},
        {"file_path": "reads_R2.fastq.gz"},
    ]
    assert [node.findtext("DataType") for node in file_nodes] == [
        "generic-data",
        "generic-data",
    ]

def test_create_submission_xml__sra_skips_blank_attributes_and_blank_bioproject(handler, config_dict, sra_metadata):
    metadata = sra_metadata.copy()
    metadata["bioproject"] = ""
    metadata["sra-library_name"] = ""
    metadata["sra-platform"] = "ILLUMINA"
    metadata["sra-loader"] = "ignored-loader"
    metadata["sra-file"] = "legacy-file-column"
    root = parse_xml(handler.create_submission_xml(organism="FLU", database="SRA", submission_name="SUB1", config_dict=config_dict, metadata=metadata))
    attrs = [(node.attrib, node.text)for node in root.findall(".//AddFiles/Attribute")]
    assert attrs == [
        ({"name": "platform"}, "ILLUMINA"),
        ({"name": "loader"}, "ignored-loader"),
        ({"name": "file"}, "legacy-file-column"),
    ]
    assert root.find(".//AddFiles/AttributeRefId[@name='BioProject']") is None
    assert root.find(".//AddFiles/AttributeRefId[@name='BioSample']") is not None

#*******************************************************************************
#                         create_raw_reads_list
#*******************************************************************************

def test_create_raw_reads_list__opens_exact_file_with_write_plus_and_writes_one_path_per_line(handler, tmp_path, monkeypatch):
    original_open = open
    observed_open_calls: list[tuple[str, str]] = []

    def fake_open(file, mode="r", *args, **kwargs):
        observed_open_calls.append((str(file), mode))
        return original_open(file, mode, *args, **kwargs)

    monkeypatch.setattr("builtins.open", fake_open)

    handler.create_raw_reads_list(str(tmp_path), {"/data/a.fastq.gz", "/data/b.fastq.gz"})

    assert observed_open_calls == [(str(tmp_path / "raw_reads_location.txt"), "w+")]

    observed = (tmp_path / "raw_reads_location.txt").read_text().splitlines()
    assert sorted(observed) == ["/data/a.fastq.gz", "/data/b.fastq.gz"]

#*******************************************************************************
#                     create_biosample_sra_submission
#*******************************************************************************

def test_create_biosample_sra_submission__sra_runs_raw_reads_manual_xml_and_save(handler, monkeypatch, config_dict, sra_metadata, tmp_path):
    calls: list[tuple[str, Any]] = []

    def fake_check_raw_read_files(**kwargs):
        calls.append(("check", kwargs))
        return {"/raw/a.fastq.gz"}

    def fake_create_raw_reads_list(**kwargs):
        calls.append(("raw_list", kwargs))

    def fake_create_manual_submission_files(**kwargs):
        calls.append(("manual", kwargs))

    def fake_create_submission_xml(**kwargs):
        calls.append(("xml", kwargs))
        return b"<Submission/>"

    def fake_save_xml(*args, **kwargs):
        calls.append(("save_xml", {"args": args, "kwargs": kwargs}))

    monkeypatch.setattr(handler, "check_raw_read_files", fake_check_raw_read_files)
    monkeypatch.setattr(handler, "create_raw_reads_list", fake_create_raw_reads_list)
    monkeypatch.setattr(handler, "create_manual_submission_files", fake_create_manual_submission_files)
    monkeypatch.setattr(handler, "create_submission_xml", fake_create_submission_xml)
    monkeypatch.setattr(handler.file_handler, "save_xml", fake_save_xml)

    handler.create_biosample_sra_submission("FLU", "SRA", "SUB1", str(tmp_path), str(tmp_path / "SRA"), config_dict, sra_metadata)
    assert [name for name, _ in calls] == ["check", "raw_list", "manual", "xml", "save_xml"]

    assert calls[0][1]["submission_name"] == "SUB1"
    assert calls[0][1]["submission_dir"] == str(tmp_path)
    assert calls[0][1]["metadata"] is sra_metadata

    assert calls[1][1] == {
        "submission_dir": str(tmp_path / "SRA"),
        "raw_files_list": {"/raw/a.fastq.gz"},
    }

    assert calls[2][1]["database"] == "SRA"
    assert calls[2][1]["submission_dir"] == str(tmp_path / "SRA")
    pd.testing.assert_frame_equal(calls[2][1]["metadata"], sra_metadata)
    assert calls[2][1]["metadata"] is not sra_metadata
    assert calls[2][1]["config_dict"] is config_dict

    assert calls[3][1]["organism"] == "FLU"
    assert calls[3][1]["database"] == "SRA"
    assert calls[3][1]["submission_name"] == "SUB1"
    assert calls[3][1]["metadata"] is sra_metadata
    assert calls[3][1]["config_dict"] is config_dict

    assert calls[4][1] == {
        "args": (b"<Submission/>", str(tmp_path / "SRA")),
        "kwargs": {},
    }

def test_create_biosample_sra_submission__biosample_skips_raw_read_steps(handler, monkeypatch, config_dict, biosample_metadata, tmp_path):
    calls: list[tuple[str, Any]] = []

    def fake_check_raw_read_files(**kwargs):
        calls.append(("check", kwargs))
        return {"/raw/a.fastq.gz"}

    def fake_create_raw_reads_list(**kwargs):
        calls.append(("raw_list", kwargs))

    def fake_create_manual_submission_files(**kwargs):
        calls.append(("manual", kwargs))

    def fake_create_submission_xml(**kwargs):
        calls.append(("xml", kwargs))
        return b"<Submission/>"

    def fake_save_xml(*args, **kwargs):
        calls.append(("save_xml", {"args": args, "kwargs": kwargs}))

    monkeypatch.setattr(handler, "check_raw_read_files", fake_check_raw_read_files)
    monkeypatch.setattr(handler, "create_raw_reads_list", fake_create_raw_reads_list)
    monkeypatch.setattr(handler, "create_manual_submission_files", fake_create_manual_submission_files)
    monkeypatch.setattr(handler, "create_submission_xml", fake_create_submission_xml)
    monkeypatch.setattr(handler.file_handler, "save_xml", fake_save_xml)

    handler.create_biosample_sra_submission("FLU", "BIOSAMPLE", "SUB1", str(tmp_path), str(tmp_path / "BIOSAMPLE"), config_dict, biosample_metadata)
    assert [name for name, _ in calls] == ["manual", "xml", "save_xml"]

    assert calls[0][1]["database"] == "BIOSAMPLE"
    assert calls[0][1]["submission_dir"] == str(tmp_path / "BIOSAMPLE")
    pd.testing.assert_frame_equal(calls[0][1]["metadata"], biosample_metadata)
    assert calls[0][1]["metadata"] is not biosample_metadata
    assert calls[0][1]["config_dict"] is config_dict

    assert calls[1][1]["organism"] == "FLU"
    assert calls[1][1]["database"] == "BIOSAMPLE"
    assert calls[1][1]["submission_name"] == "SUB1"
    assert calls[1][1]["metadata"] is biosample_metadata
    assert calls[1][1]["config_dict"] is config_dict

    assert calls[2][1] == {
        "args": (b"<Submission/>", str(tmp_path / "BIOSAMPLE")),
        "kwargs": {},
    }

#*******************************************************************************
#                       process_biosample_sra_report
#*******************************************************************************

def test_process_biosample_sra_report__no_action_returns_header_status_without_update(handler, monkeypatch):

    def fake_process_report_header(report_file):
        return({"SubmissionStatus": {"@status": "submitted"}}, "SUBMITTED", "SUB1")

    monkeypatch.setattr(handler.ncbi_handler, "process_report_header", fake_process_report_header)
    update_calls: list[dict[str, Any]] = []

    def capture_update_submission_status_csv(**kwargs):
        update_calls.append(kwargs)

    monkeypatch.setattr(handler.upload_log, "update_submission_status_csv", capture_update_submission_status_csv)

    assert handler.process_biosample_sra_report("report.xml", "BIOSAMPLE", "/sub") == ("SUBMITTED", "SUB1")
    assert update_calls == []

def test_process_biosample_sra_report__dict_response_requires_both_accession_and_spuid(handler, monkeypatch, capsys):
    report_dict = {
        "SubmissionStatus": {
            "Action": {
                "@target_db": "BioSample",
                "@status": "processed-ok",
                "Response": {
                    "Object": {
                        "@spuid": "BS1",
                    }
                },
            }
        }
    }

    def fake_process_report_header(report_file):
        return report_dict, "PROCESSED", "SUB1"

    monkeypatch.setattr(handler.ncbi_handler, "process_report_header", fake_process_report_header)
    update_calls: list[dict[str, Any]] = []

    def capture_update_submission_status_csv(**kwargs):
        update_calls.append(kwargs)

    monkeypatch.setattr(handler.upload_log, "update_submission_status_csv", capture_update_submission_status_csv)
    assert handler.process_biosample_sra_report("report.xml", "BIOSAMPLE", "/sub") == ("PROCESSED", "SUB1")
    assert update_calls == []
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == "Error: Unable to process BIOSAMPLE report.xml to retrieve accessions at: report.xml\n"

def test_process_biosample_sra_report__list_action_skips_wrong_db_and_missing_response_then_uses_valid_sra(handler, monkeypatch):
    report_dict = {
        "SubmissionStatus": {
            "Action": [
                {
                    "@target_db": "BioSample",
                    "@status": "processed-ok",
                    "Response": {
                        "Object": {
                            "@spuid": "BS1",
                            "@accession": "SAMN000001",
                        }
                    },
                },
                {
                    "@target_db": "SRA",
                    "@status": "processed-ok",
                },
                {
                    "@target_db": "SRA",
                    "@status": "processed-ok",
                    "Response": [
                        {"Warning": "skip me"},
                        {
                            "@status": "processed-ok",
                            "Object": {
                                "@spuid": "SRA1",
                                "@accession": "SRR000001",
                            },
                        },
                    ],
                },
            ]
        }
    }

    def fake_process_report_header(report_file):
        return report_dict, "PROCESSED", "SUB1"

    monkeypatch.setattr(handler.ncbi_handler, "process_report_header", fake_process_report_header)
    update_calls: list[dict[str, Any]] = []

    def capture_update_submission_status_csv(**kwargs):
        update_calls.append(kwargs)

    monkeypatch.setattr(handler.upload_log, "update_submission_status_csv", capture_update_submission_status_csv)
    assert handler.process_biosample_sra_report("report.xml", "SRA", "/sub") == ("PROCESSED", "SUB1")
    assert update_calls[0]["submission_dir"] == "/sub"
    assert update_calls[0]["update_database"] == "SRA"
    assert update_calls[0]["update_df"].to_dict("records") == [
        {
            "sra-sample_name": "SRA1",
            "sra_status": "processed-ok",
            "sra_accession": "SRR000001",
            "sra_message": "",
        }
    ]

def test_process_biosample_sra_report__dict_action_updates_sra_status(handler, monkeypatch, capsys):
    report_dict = {
        "SubmissionStatus": {
            "Action": {
                "@target_db": "SRA",
                "@status": "processed-ok",
                "Response": {
                    "Object": {
                        "@spuid": "SRA1",
                        "@accession": "SRR000001",
                    }
                },
            }
        }
    }

    def fake_process_report_header(report_file):
        return report_dict, "PROCESSED", "SUB1"

    monkeypatch.setattr(handler.ncbi_handler, "process_report_header", fake_process_report_header)
    update_calls: list[dict[str, Any]] = []

    def capture_update_submission_status_csv(**kwargs):
        update_calls.append(kwargs)

    monkeypatch.setattr(handler.upload_log, "update_submission_status_csv", capture_update_submission_status_csv)
    result = handler.process_biosample_sra_report("report.xml", "SRA", "/sub")
    assert result == ("PROCESSED", "SUB1")
    assert update_calls[0]["submission_dir"] == "/sub"
    assert update_calls[0]["update_database"] == "SRA"
    assert update_calls[0]["update_df"].to_dict("records") == [
        {
            "sra-sample_name": "SRA1",
            "sra_status": "processed-ok",
            "sra_accession": "SRR000001",
            "sra_message": "",
        }
    ]
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""

def test_process_biosample_sra_report__dict_action_updates_biosample_status(handler, monkeypatch, capsys):
    report_dict = {
        "SubmissionStatus": {
            "Action": {
                "@target_db": "BioSample",
                "@status": "processed-ok",
                "Response": {"Object": {"@spuid": "BS1", "@accession": "SAMN000001"}},
            }
        }
    }

    def fake_process_report_header(report_file):
        return(report_dict, "PROCESSED", "SUB1")

    monkeypatch.setattr(handler.ncbi_handler, "process_report_header", fake_process_report_header)
    update_calls: list[dict[str, Any]] = []

    def capture_update_submission_status_csv(**kwargs):
        update_calls.append(kwargs)

    monkeypatch.setattr(handler.upload_log, "update_submission_status_csv", capture_update_submission_status_csv)

    status, submission_id = handler.process_biosample_sra_report("report.xml", "BIOSAMPLE", "/sub")

    assert (status, submission_id) == ("PROCESSED", "SUB1")
    assert update_calls[0]["submission_dir"] == "/sub"
    assert update_calls[0]["update_database"] == "BIOSAMPLE"
    update_df = update_calls[0]["update_df"]
    assert update_df.to_dict("records") == [
        {
            "bs-sample_name": "BS1",
            "biosample_status": "processed-ok",
            "biosample_accession": "SAMN000001",
            "biosample_message": "",
        }
    ]
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""

def test_process_biosample_sra_report__list_action_updates_sra_status_and_skips_other_db(handler, monkeypatch):
    report_dict = {
        "SubmissionStatus": {
            "Action": [
                {"@target_db": "BioSample", "@status": "processed-ok", "Response": {"Object": {"@spuid": "BS1", "@accession": "SAMN000001"}}},
                {
                    "@target_db": "SRA",
                    "@status": "processed-ok",
                    "Response": [
                        {"Warning": "not sample info"},
                        {"@status": "ignored-response-status", "Object": {"@spuid": "SRA1", "@accession": "SRR000001"}},
                    ],
                },
            ]
        }
    }

    def fake_process_report_header(report_file):
        return(report_dict, "PROCESSED", "SUB1")

    monkeypatch.setattr(handler.ncbi_handler, "process_report_header", fake_process_report_header)
    update_calls: list[dict[str, Any]] = []

    def capture_update_submission_status_csv(**kwargs):
        update_calls.append(kwargs)

    monkeypatch.setattr(handler.upload_log, "update_submission_status_csv", capture_update_submission_status_csv)

    status, submission_id = handler.process_biosample_sra_report("report.xml", "SRA", "/sub")

    assert (status, submission_id) == ("PROCESSED", "SUB1")
    update_df = update_calls[0]["update_df"]
    assert update_df.to_dict("records") == [
        {
            "sra-sample_name": "SRA1",
            "sra_status": "ignored-response-status",
            "sra_accession": "SRR000001",
            "sra_message": "",
        }
    ]

def test_process_biosample_sra_report__invalid_action_type_returns_without_update(handler, monkeypatch, capsys):
    report_dict = {"SubmissionStatus": {"Action": "bad-action"}}

    def fake_process_report_header(report_file):
        return(report_dict, "SUBMITTED", "SUB1")

    monkeypatch.setattr(handler.ncbi_handler, "process_report_header", fake_process_report_header)
    update_calls: list[dict[str, Any]] = []

    def capture_update_submission_status_csv(**kwargs):
        update_calls.append(kwargs)

    monkeypatch.setattr(handler.upload_log, "update_submission_status_csv", capture_update_submission_status_csv)

    assert handler.process_biosample_sra_report("report.xml", "BIOSAMPLE", "/sub") == ("SUBMITTED", "SUB1")
    assert update_calls == []
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == "Error: Unable to correctly process BioSample report at: report.xml\n"

def test_process_biosample_sra_report__processed_without_accessions_warns_but_does_not_update(handler, monkeypatch, capsys):
    report_dict = {"SubmissionStatus": {"Action": {"@target_db": "BioSample", "@status": "processed-ok", "Response": {}}}}

    def fake_process_report_header(report_file):
        return(report_dict, "PROCESSED", "SUB1")

    monkeypatch.setattr(handler.ncbi_handler, "process_report_header", fake_process_report_header)
    update_calls: list[dict[str, Any]] = []

    def capture_update_submission_status_csv(**kwargs):
        update_calls.append(kwargs)

    monkeypatch.setattr(handler.upload_log, "update_submission_status_csv", capture_update_submission_status_csv)

    assert handler.process_biosample_sra_report("report.xml", "BIOSAMPLE", "/sub") == ("PROCESSED", "SUB1")
    assert update_calls == []
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == "Error: Unable to process BIOSAMPLE report.xml to retrieve accessions at: report.xml\n"

def test_process_biosample_sra_report__malformed_action_is_swallowed_and_returns_header_status(handler, monkeypatch):
    report_dict = {"SubmissionStatus": {"Action": [{"@target_db": "SRA", "@status": "processed-ok", "Response": [None]}]}}

    def fake_process_report_header(report_file):
        return(report_dict, "PROCESSING", "SUB1")

    monkeypatch.setattr(handler.ncbi_handler, "process_report_header", fake_process_report_header)
    update_calls: list[dict[str, Any]] = []

    def capture_update_submission_status_csv(**kwargs):
        update_calls.append(kwargs)

    monkeypatch.setattr(handler.upload_log, "update_submission_status_csv", capture_update_submission_status_csv)
    assert handler.process_biosample_sra_report("report.xml", "SRA", "/sub") == ("PROCESSING", "SUB1")
    assert update_calls == []

def test_process_biosample_sra_report__list_response_ignores_list_without_object_and_does_not_update(handler, monkeypatch, capsys):
    report_dict = {
        "SubmissionStatus": {
            "Action": {
                "@target_db": "SRA",
                "@status": "processed-ok",
                "Response": [
                    {"Warning": "not sample info"},
                    {"Message": "still not sample info"},
                ],
            }
        }
    }

    def fake_process_report_header(report_file):
        return report_dict, "PROCESSED", "SUB1"

    monkeypatch.setattr(handler.ncbi_handler, "process_report_header", fake_process_report_header)
    update_calls: list[dict[str, Any]] = []

    def capture_update_submission_status_csv(**kwargs):
        update_calls.append(kwargs)

    monkeypatch.setattr(handler.upload_log, "update_submission_status_csv", capture_update_submission_status_csv)
    assert handler.process_biosample_sra_report("report.xml", "SRA", "/sub") == ("PROCESSED", "SUB1")
    assert update_calls == []
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == "Error: Unable to process SRA report.xml to retrieve accessions at: report.xml\n"

def test_process_biosample_sra_report__dict_response_without_object_does_not_update(handler, monkeypatch, capsys):
    report_dict = {
        "SubmissionStatus": {
            "Action": {
                "@target_db": "SRA",
                "@status": "processed-ok",
                "Response": {"Message": "no object"},
            }
        }
    }

    def fake_process_report_header(report_file):
        return report_dict, "PROCESSED", "SUB1"

    monkeypatch.setattr(handler.ncbi_handler, "process_report_header", fake_process_report_header)
    update_calls: list[dict[str, Any]] = []

    def capture_update_submission_status_csv(**kwargs):
        update_calls.append(kwargs)

    monkeypatch.setattr(handler.upload_log, "update_submission_status_csv", capture_update_submission_status_csv)
    assert handler.process_biosample_sra_report("report.xml", "SRA", "/sub") == ("PROCESSED", "SUB1")
    assert update_calls == []
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == "Error: Unable to process SRA report.xml to retrieve accessions at: report.xml\n"

def test_process_biosample_sra_report__dict_object_requires_accession_and_spuid(handler, monkeypatch, capsys):
    report_dict = {
        "SubmissionStatus": {
            "Action": {
                "@target_db": "SRA",
                "@status": "processed-ok",
                "Response": {"Object": {"@accession": "SRR000001"}},
            }
        }
    }

    def fake_process_report_header(report_file):
        return report_dict, "PROCESSED", "SUB1"

    monkeypatch.setattr(handler.ncbi_handler, "process_report_header", fake_process_report_header)
    update_calls: list[dict[str, Any]] = []

    def capture_update_submission_status_csv(**kwargs):
        update_calls.append(kwargs)

    monkeypatch.setattr(handler.upload_log, "update_submission_status_csv", capture_update_submission_status_csv)
    assert handler.process_biosample_sra_report("report.xml", "SRA", "/sub") == ("PROCESSED", "SUB1")
    assert update_calls == []
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == "Error: Unable to process SRA report.xml to retrieve accessions at: report.xml\n"
