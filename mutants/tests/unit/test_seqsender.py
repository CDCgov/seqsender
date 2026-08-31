from __future__ import annotations

import argparse
import importlib.util
import os
import sys
import types
from pathlib import Path

import pandas as pd
import pytest
from typing import Any

SOURCE_DIR = Path(__file__).resolve().parents[2]
MODULE_PATH = os.path.join(SOURCE_DIR, "seqsender.py")

"""Unit tests initially generated with ChatGPT v5.5 "Deep Research" and "Thinking"
then modified for use and clarity with validation and coverage testing via mutmut."""

#*******************************************************************************
#                    create test seqsender.py connections
#*******************************************************************************

@pytest.fixture()
def seqsender_module(monkeypatch: pytest.MonkeyPatch):
    calls: list[tuple[str, dict]] = []

    setup: Any = types.ModuleType("setup")

    def create_test_data(**kwargs):
        calls.append(("setup.create_test_data", kwargs))
        return None

    def download_biosample_xml_list():
        calls.append(("setup.download_biosample_xml_list", {}))
        return None

    def test_internet_connection(**kwargs):
        calls.append(("setup.test_internet_connection", kwargs))
        return None

    setup.create_test_data = create_test_data
    setup.download_biosample_xml_list = download_biosample_xml_list
    setup.test_internet_connection = test_internet_connection

    file_handler: Any = types.ModuleType("file_handler")

    def validate_directory(**kwargs):
        calls.append(("file_handler.validate_directory", kwargs))

    def validate_file(**kwargs):
        calls.append(("file_handler.validate_file", kwargs))

    def create_directory(path):
        calls.append(("file_handler.create_directory", {"path": path}))

    def validate_gisaid_installer(*args, **kwargs):
        calls.append(
            (
                "file_handler.validate_gisaid_installer",
                {"args": args, "kwargs": kwargs},
            )
        )
    file_handler.validate_directory = validate_directory
    file_handler.validate_file = validate_file
    file_handler.create_directory = create_directory
    file_handler.validate_gisaid_installer = validate_gisaid_installer

    argument_handler: Any = types.ModuleType("argument_handler")

    def args_parser() -> argparse.ArgumentParser:
        return argparse.ArgumentParser()

    argument_handler.args_parser = args_parser
    ncbi_handler: Any = types.ModuleType("ncbi_handler")

    def submit_ncbi(**kwargs):
        calls.append(("ncbi_handler.submit_ncbi", kwargs))

    def email_table2asn(**kwargs):
        calls.append(("ncbi_handler.email_table2asn", kwargs))
        return "EMAILED"

    ncbi_handler.submit_ncbi = submit_ncbi
    ncbi_handler.email_table2asn = email_table2asn

    genbank_handler: Any = types.ModuleType("genbank_handler")

    def create_genbank_submission(**kwargs):
        calls.append(("genbank_handler.create_genbank_submission", kwargs))

    genbank_handler.create_genbank_submission = create_genbank_submission

    biosample_sra_handler: Any = types.ModuleType("biosample_sra_handler")

    def create_biosample_sra_submission(**kwargs):
        calls.append(("biosample_sra_handler.create_biosample_sra_submission", kwargs))

    biosample_sra_handler.create_biosample_sra_submission = create_biosample_sra_submission

    gisaid_handler: Any = types.ModuleType("gisaid_handler")

    def create_gisaid_files(**kwargs):
        calls.append(("gisaid_handler.create_gisaid_files", kwargs))

    def submit_gisaid(**kwargs):
        calls.append(("gisaid_handler.submit_gisaid", kwargs))

    gisaid_handler.create_gisaid_files = create_gisaid_files
    gisaid_handler.submit_gisaid = submit_gisaid

    upload_log: Any = types.ModuleType("upload_log")

    def create_submission_status_csv(**kwargs):
        calls.append(("upload_log.create_submission_status_csv", kwargs))

    def create_submission_log(**kwargs):
        calls.append(("upload_log.create_submission_log", kwargs))

    def update_submission_status(**kwargs):
        calls.append(("upload_log.update_submission_status", kwargs))

    upload_log.create_submission_status_csv = create_submission_status_csv
    upload_log.create_submission_log = create_submission_log
    upload_log.update_submission_status = update_submission_status

    tools: Any = types.ModuleType("tools")

    def get_config(**kwargs):
        return {
            "NCBI": {
                "Link_Sample_Between_NCBI_Databases": False,
                "Username": "ncbi-user",
            },
            "GISAID": {"Username": "gisaid-user"},
        }
    def get_metadata(**kwargs):
        return pd.DataFrame(
            {
                "sequence_name": ["seq1"],
                "gb-sample_name": ["gb1"],
                "gs-sample_name": ["gs1"],
            }
        )
    def process_fasta_samples(**kwargs):
        return kwargs["metadata"].assign(fasta_sequence_orig=["ATGC"])

    def get_submission_type(test):
        return "TEST" if test else "PRODUCTION"

    def get_submission_position(config_dict, database):
        return None

    tools.get_config = get_config
    tools.get_metadata = get_metadata
    tools.process_fasta_samples = process_fasta_samples
    tools.get_submission_type = get_submission_type
    tools.get_submission_position = get_submission_position

    settings: Any = types.ModuleType("settings")
    settings.VERSION = "9.9.9-test"
    settings.GENBANK_FTP_ORGANISMS = ["FLU", "COV"]

    # Clear cached real modules/stubs so seqsender.py imports only these test doubles.
    for module_name in [
        "src",
        "src.setup",
        "src.file_handler",
        "src.argument_handler",
        "src.ncbi_handler",
        "src.genbank_handler",
        "src.biosample_sra_handler",
        "src.gisaid_handler",
        "src.upload_log",
        "src.tools",
        "src.settings",
        "setup",
        "file_handler",
        "argument_handler",
        "ncbi_handler",
        "genbank_handler",
        "biosample_sra_handler",
        "gisaid_handler",
        "upload_log",
        "tools",
        "settings",
        "seqsender_under_test",
    ]:
        sys.modules.pop(module_name, None)

    # seqsender.py may use either old top-level imports (`import setup`) or
    # package-style imports (`from src import setup`). Register both names.
    src_pkg: Any = types.ModuleType("src")
    src_pkg.__path__ = []  # Mark as a package for `from src import ...`.
    monkeypatch.setitem(sys.modules, "src", src_pkg)

    def alias_src_module(name, module):
        monkeypatch.setitem(sys.modules, name, module)
        monkeypatch.setitem(sys.modules, f"src.{name}", module)
        setattr(src_pkg, name, module)

    module_map: dict[str, Any] = {
        "setup": setup,
        "file_handler": file_handler,
        "argument_handler": argument_handler,
        "ncbi_handler": ncbi_handler,
        "genbank_handler": genbank_handler,
        "biosample_sra_handler": biosample_sra_handler,
        "gisaid_handler": gisaid_handler,
        "upload_log": upload_log,
        "tools": tools,
        "settings": settings,
    }

    for name, stub in module_map.items():
        alias_src_module(name, stub)

    spec = importlib.util.spec_from_file_location("seqsender_under_test", MODULE_PATH)
    assert spec and spec.loader
    module: Any = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    module._test_calls = calls
    return module


def called(module, name: str):
    return [kwargs for call_name, kwargs in module._test_calls if call_name == name]

#*******************************************************************************
#                            get_execution_time
#*******************************************************************************

def test_get_execution_time__prints_runtime(seqsender_module, capsys):
    seqsender_module.get_execution_time()
    assert "Total runtime (HRS:MIN:SECS):" in capsys.readouterr().out

#*******************************************************************************
#                                    prep
#*******************************************************************************

def test_prep__biosample_and_sra_creates_expected_submissions(seqsender_module, tmp_path, capsys):
    result = seqsender_module.prep(
        database=["BIOSAMPLE", "SRA"],
        organism="COV",
        submission_dir=str(tmp_path),
        submission_name="sub1",
        config_file=str(tmp_path / "config.yaml"),
        metadata_file=str(tmp_path / "metadata.csv"),
        fasta_file=None,
        gff_file=None,
        table2asn=False,
        publication_title=None,
        publication_status=None,
        decrypt_key="test-key",
    )

    assert result[0] == str(tmp_path / "config.yaml")
    assert len(called(seqsender_module, "biosample_sra_handler.create_biosample_sra_submission")) == 2
    assert len(called(seqsender_module, "file_handler.create_directory")) == 2
    assert "turned off" in capsys.readouterr().out

def test_prep__sra_without_biosample_warns(seqsender_module, tmp_path, capsys):
    seqsender_module.prep(
        database=["SRA"],
        organism="OTHER",
        submission_dir=str(tmp_path),
        submission_name="sub1",
        config_file=str(tmp_path / "config.yaml"),
        metadata_file=str(tmp_path / "metadata.csv"),
        fasta_file=None,
        gff_file=None,
        table2asn=False,
        publication_title=None,
        publication_status=None,
        decrypt_key="test-key",
    )
    assert "SRA requires a BioSample submission" in capsys.readouterr().out

@pytest.mark.parametrize("database", [["GENBANK"], ["GISAID"]])
def test_prep_requires_fasta_for_genbank_or_gisaid(seqsender_module, tmp_path, database):
    with pytest.raises(SystemExit):
        seqsender_module.prep(
            database=database,
            organism="FLU",
            submission_dir=str(tmp_path),
            submission_name="sub1",
            config_file="config.yaml",
            metadata_file="metadata.csv",
            fasta_file=None,
            gff_file=None,
            table2asn=False,
            publication_title=None,
            publication_status=None,
            decrypt_key="test-key",
        )

def test_prep__genbank_and_gisaid_processes_fasta_and_routes_handlers(seqsender_module, tmp_path):
    seqsender_module.prep(
        database=["GENBANK", "GISAID"],
        organism="FLU",
        submission_dir=str(tmp_path),
        submission_name="sub1",
        config_file=str(tmp_path / "config.yaml"),
        metadata_file=str(tmp_path / "metadata.csv"),
        fasta_file=str(tmp_path / "seqs.fasta"),
        gff_file=str(tmp_path / "ann.gff"),
        table2asn=True,
        publication_title="Title",
        publication_status="Published",
        decrypt_key="test-key",
    )

    assert len(called(seqsender_module, "genbank_handler.create_genbank_submission")) == 1
    assert len(called(seqsender_module, "gisaid_handler.create_gisaid_files")) == 1
    genbank_call = called(seqsender_module, "genbank_handler.create_genbank_submission")[0]
    assert genbank_call["gff_file"] == str(tmp_path / "ann.gff")
    assert genbank_call["table2asn"] is True
    assert "fasta_sequence_orig" in genbank_call["metadata"].columns

def test_prep__relative_file_names_are_resolved_under_submission_folder(seqsender_module, tmp_path):
    seqsender_module.prep(
        database=["BIOSAMPLE"],
        organism="COV",
        submission_dir=str(tmp_path),
        submission_name="sub1",
        config_file="config.yaml",
        metadata_file="metadata.csv",
        fasta_file=None,
        gff_file=None,
        table2asn=False,
        publication_title=None,
        publication_status=None,
        decrypt_key="test-key",
    )
    validate_file_calls = called(seqsender_module, "file_handler.validate_file")
    assert {Path(call["file_path"]).name for call in validate_file_calls} == {"config.yaml", "metadata.csv"}
    assert all(str(tmp_path / "sub1") in call["file_path"] for call in validate_file_calls)

def test_prep__invalid_database_exits(seqsender_module, tmp_path):
    with pytest.raises(SystemExit):
        seqsender_module.prep(
            database=["BAD"],
            organism="COV",
            submission_dir=str(tmp_path),
            submission_name="sub1",
            config_file=str(tmp_path / "config.yaml"),
            metadata_file=str(tmp_path / "metadata.csv"),
            fasta_file=None,
            gff_file=None,
            table2asn=False,
            publication_title=None,
            publication_status=None,
            decrypt_key="test-key",
        )

#*******************************************************************************
#                                    submit
#*******************************************************************************

def test_submit__biosample_and_sra_submit_to_ncbi_and_log(seqsender_module, tmp_path, monkeypatch):
    def fake_prep(**kwargs):
        return (
            "config.yaml",
            {"NCBI": {"Link_Sample_Between_NCBI_Databases": False}, "GISAID": {}},
            pd.DataFrame({"sample": ["x"]}),
        )
    monkeypatch.setattr(seqsender_module, "prep", fake_prep)
    seqsender_module.submit(
        database=["BIOSAMPLE", "SRA"],
        organism="COV",
        submission_dir=str(tmp_path),
        submission_name="sub1",
        config_file="config.yaml",
        metadata_file="metadata.csv",
        fasta_file=None,
        gff_file=None,
        publication_title=None,
        publication_status=None,
        test=True,
        decrypt_key="test-key",
    )
    assert len(called(seqsender_module, "ncbi_handler.submit_ncbi")) == 2
    assert len(called(seqsender_module, "upload_log.create_submission_log")) == 2
    assert all(call["submission_status"] == "SUBMITTED" for call in called(seqsender_module, "upload_log.create_submission_log"))
    assert all(call["submission_type"] == "TEST" for call in called(seqsender_module, "upload_log.create_submission_log"))

def test_submit__genbank_ftp(seqsender_module, tmp_path, monkeypatch):
    def fake_prep(**kwargs):
        return (
            "config.yaml",
            {"NCBI": {"Link_Sample_Between_NCBI_Databases": False}, "GISAID": {}},
            pd.DataFrame({"sample": ["x"]}),
        )
    monkeypatch.setattr(seqsender_module, "prep", fake_prep)
    seqsender_module.submit(
        database=["GENBANK"],
        organism="FLU",
        submission_dir=str(tmp_path),
        submission_name="sub1",
        config_file="config.yaml",
        metadata_file="metadata.csv",
        fasta_file="seqs.fasta",
        gff_file=None,
        publication_title=None,
        publication_status=None,
        decrypt_key="test-key",
    )
    ncbi_calls = called(seqsender_module, "ncbi_handler.submit_ncbi")
    assert ncbi_calls[0]["database"] == "GENBANK"
    log = called(seqsender_module, "upload_log.create_submission_log")[0]
    assert log["database"] == "GENBANK-FTP"
    assert log["submission_status"] == "SUBMITTED"

def test_submit__genbank_table2asn_forced_by_flag(seqsender_module, tmp_path, monkeypatch):
    def fake_prep(**kwargs):
        return (
            "config.yaml",
            {"NCBI": {"Link_Sample_Between_NCBI_Databases": False}, "GISAID": {}},
            pd.DataFrame({"sample": ["x"]}),
        )
    monkeypatch.setattr(seqsender_module, "prep", fake_prep)
    seqsender_module.submit(
        database=["GENBANK"],
        organism="FLU",
        submission_dir=str(tmp_path),
        submission_name="sub1",
        config_file="config.yaml",
        metadata_file="metadata.csv",
        fasta_file="seqs.fasta",
        gff_file=None,
        publication_title=None,
        publication_status=None,
        table2asn=True,
        decrypt_key="test-key",
    )
    assert len(called(seqsender_module, "ncbi_handler.email_table2asn")) == 1
    log = called(seqsender_module, "upload_log.create_submission_log")[0]
    assert log["database"] == "GENBANK-TBL2ASN"
    assert log["submission_status"] == "EMAILED"

def test_submit__genbank_table2asn_for_non_ftp_organism(seqsender_module, tmp_path, monkeypatch):
    def fake_prep(**kwargs):
        return (
            "config.yaml",
            {"NCBI": {"Link_Sample_Between_NCBI_Databases": False}, "GISAID": {}},
            pd.DataFrame({"sample": ["x"]}),
        )
    monkeypatch.setattr(seqsender_module, "prep", fake_prep)
    seqsender_module.submit(
        database=["GENBANK"],
        organism="OTHER",
        submission_dir=str(tmp_path),
        submission_name="sub1",
        config_file="config.yaml",
        metadata_file="metadata.csv",
        fasta_file="seqs.fasta",
        gff_file=None,
        publication_title=None,
        publication_status=None,
        decrypt_key="test-key",
    )
    assert len(called(seqsender_module, "ncbi_handler.email_table2asn")) == 1
    assert called(seqsender_module, "upload_log.create_submission_log")[0]["database"] == "GENBANK-TBL2ASN"


def test_submit__genbank_waits_when_linked_to_ncbi_first(seqsender_module, tmp_path, monkeypatch):
    def fake_prep(**kwargs):
        return (
            "config.yaml",
            {"NCBI": {"Link_Sample_Between_NCBI_Databases": True}, "GISAID": {}},
            pd.DataFrame({"sample": ["x"]}),
        )
    monkeypatch.setattr(seqsender_module, "prep", fake_prep)
    seqsender_module.submit(
        database=["BIOSAMPLE", "GENBANK"],
        organism="FLU",
        submission_dir=str(tmp_path),
        submission_name="sub1",
        config_file="config.yaml",
        metadata_file="metadata.csv",
        fasta_file="seqs.fasta",
        gff_file=None,
        publication_title=None,
        publication_status=None,
        decrypt_key="test-key",
    )
    assert len(called(seqsender_module, "ncbi_handler.submit_ncbi")) == 1  # only BioSample now
    genbank_log = called(seqsender_module, "upload_log.create_submission_log")[1]
    assert genbank_log["database"] == "GENBANK-FTP"
    assert genbank_log["submission_status"] == "WAITING"


def test_submit__gisaid_validates_cli_and_submits(seqsender_module, tmp_path, monkeypatch):
    def fake_prep(**kwargs):
        return (
            "config.yaml",
            {"NCBI": {"Link_Sample_Between_NCBI_Databases": False}, "GISAID": {"CLI_Path": "/bin/cli"}},
            pd.DataFrame({"sample": ["x"]}),
        )
    monkeypatch.setattr(seqsender_module, "prep", fake_prep)
    seqsender_module.submit(
        database=["GISAID"],
        organism="COV",
        submission_dir=str(tmp_path),
        submission_name="sub1",
        config_file="config.yaml",
        metadata_file="metadata.csv",
        fasta_file="seqs.fasta",
        gff_file=None,
        publication_title=None,
        publication_status=None,
        decrypt_key="test-key",
    )
    assert len(called(seqsender_module, "file_handler.validate_gisaid_installer")) == 1
    assert len(called(seqsender_module, "gisaid_handler.submit_gisaid")) == 1
    assert called(seqsender_module, "upload_log.create_submission_log")[0]["submission_status"] is None


def test_submit__gisaid_waits_if_genbank_first(seqsender_module, tmp_path, monkeypatch):
    def fake_prep(**kwargs):
        return (
            "config.yaml",
            {"NCBI": {"Link_Sample_Between_NCBI_Databases": False}, "GISAID": {}},
            pd.DataFrame({"sample": ["x"]}),
        )

    def fake_get_submission_position(config_dict, database):
        return 2 if database == "GISAID" else 1

    monkeypatch.setattr(seqsender_module, "prep", fake_prep)
    monkeypatch.setattr(seqsender_module.tools, "get_submission_position", fake_get_submission_position)
    seqsender_module.submit(
        database=["GENBANK", "GISAID"],
        organism="FLU",
        submission_dir=str(tmp_path),
        submission_name="sub1",
        config_file="config.yaml",
        metadata_file="metadata.csv",
        fasta_file="seqs.fasta",
        gff_file=None,
        publication_title=None,
        publication_status=None,
        decrypt_key="test-key",
    )
    assert len(called(seqsender_module, "gisaid_handler.submit_gisaid")) == 0
    gisaid_log = called(seqsender_module, "upload_log.create_submission_log")[1]
    assert gisaid_log["database"] == "GISAID"
    assert gisaid_log["submission_status"] == "WAITING"


def test_submit__invalid_database_exits(seqsender_module, tmp_path, monkeypatch):
    def fake_prep(**kwargs):
        return (
            "config.yaml",
            {"NCBI": {}, "GISAID": {}},
            pd.DataFrame({"sample": ["x"]}),
        )
    monkeypatch.setattr(seqsender_module, "prep", fake_prep)
    with pytest.raises(SystemExit):
        seqsender_module.submit(
            database=["BAD"],
            organism="COV",
            submission_dir=str(tmp_path),
            submission_name="sub1",
            config_file="config.yaml",
            metadata_file="metadata.csv",
            fasta_file=None,
            gff_file=None,
            publication_title=None,
            publication_status=None,
            decrypt_key="test-key",
        )

#*******************************************************************************
#                            main argument parser
#*******************************************************************************

class FakeParser:
    def __init__(self, args):
        self.args = args
        self.help_printed = False

    def parse_args(self):
        return self.args

    def print_help(self):
        self.help_printed = True

def run_main_with_args(seqsender_module, monkeypatch, args):
    parser = FakeParser(args)

    def fake_args_parser():
        return parser

    monkeypatch.setattr(seqsender_module.argument_handler, "args_parser", fake_args_parser)
    return parser

def test_main__dispatches_prep(seqsender_module, tmp_path, monkeypatch):
    captured = {}

    def fake_prep(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(seqsender_module, "prep", fake_prep)
    args = argparse.Namespace(
        command="prep",
        biosample="BIOSAMPLE",
        sra="",
        genbank="",
        gisaid="",
        organism="COV",
        submission_name="sub1",
        submission_dir=str(tmp_path),
        config_file="config.yaml",
        metadata_file="metadata.csv",
        fasta_file=None,
        gff_file=None,
        table2asn=False,
        skip_validation=False,
        publication_title=None,
        publication_status=None,
    )
    run_main_with_args(seqsender_module, monkeypatch, args)
    seqsender_module.main()
    assert captured["database"] == ["BIOSAMPLE"]
    assert captured["submission_dir"] == os.path.abspath(str(tmp_path))

def test_main__dispatches_submit(seqsender_module, tmp_path, monkeypatch):
    captured = {}

    def fake_submit(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(seqsender_module, "submit", fake_submit)
    args = argparse.Namespace(
        command="submit",
        biosample="",
        sra="SRA",
        genbank="",
        gisaid="",
        organism="OTHER",
        submission_name="sub1",
        submission_dir=str(tmp_path),
        config_file="config.yaml",
        metadata_file="metadata.csv",
        fasta_file=None,
        gff_file=None,
        table2asn=False,
        test=True,
        skip_validation=True,
        publication_title=None,
        publication_status=None,
        key="test-key",
    )
    run_main_with_args(seqsender_module, monkeypatch, args)
    seqsender_module.main()
    assert captured["database"] == ["SRA"]
    assert captured["test"] is True
    assert captured["skip_validation"] is True

def test_main__missing_database_prints_help_and_exits(seqsender_module, tmp_path, monkeypatch):
    args = argparse.Namespace(
        command="prep",
        biosample="",
        sra="",
        genbank="",
        gisaid="",
        organism="COV",
        submission_name="sub1",
        submission_dir=str(tmp_path),
        config_file="config.yaml",
        metadata_file="metadata.csv",
        fasta_file=None,
        gff_file=None,
        table2asn=False,
        skip_validation=False,
        publication_title=None,
        publication_status=None,
    )
    parser = run_main_with_args(seqsender_module, monkeypatch, args)
    with pytest.raises(SystemExit):
        seqsender_module.main()
    assert parser.help_printed is True

def test_main__submission_status_dispatch(seqsender_module, tmp_path, monkeypatch):
    args = argparse.Namespace(command="submission_status", submission_dir=str(tmp_path), submission_name="sub1", key="test-key")
    run_main_with_args(seqsender_module, monkeypatch, args)
    seqsender_module.main()
    call = called(seqsender_module, "upload_log.update_submission_status")[0]
    assert call["submission_dir"] == os.path.abspath(str(tmp_path))
    assert call["submission_name"] == "sub1"

def test_main__test_data_dispatch(seqsender_module, tmp_path, monkeypatch):
    args = argparse.Namespace(
        command="test_data",
        biosample="BIOSAMPLE",
        sra="SRA",
        genbank="",
        gisaid="",
        organism="FLU",
        submission_dir=str(tmp_path),
    )
    run_main_with_args(seqsender_module, monkeypatch, args)
    seqsender_module.main()
    call = called(seqsender_module, "setup.create_test_data")[0]
    assert call["database"] == ["BIOSAMPLE", "SRA"]
    assert call["submission_dir"] == os.path.abspath(str(tmp_path))

def test_main__version_prints_version(seqsender_module, monkeypatch, capsys):
    args = argparse.Namespace(command="version")
    run_main_with_args(seqsender_module, monkeypatch, args)
    seqsender_module.main()
    assert "Version: 9.9.9-test" in capsys.readouterr().out

def test_main__update_biosample_dispatch(seqsender_module, monkeypatch):
    args = argparse.Namespace(command="update_biosample")
    run_main_with_args(seqsender_module, monkeypatch, args)
    seqsender_module.main()
    assert len(called(seqsender_module, "setup.download_biosample_xml_list")) == 1

def test_main__test_network_connection_dispatch(seqsender_module, monkeypatch):
    args = argparse.Namespace(command="test_network_connection")
    run_main_with_args(seqsender_module, monkeypatch, args)
    seqsender_module.main()
    call = called(seqsender_module, "setup.test_internet_connection")[0]
    assert call["databases"] == ["GENERAL", "NCBI", "GISAID"]

def test_main__unknown_command_prints_help_and_exits(seqsender_module, monkeypatch):
    parser = run_main_with_args(seqsender_module, monkeypatch, argparse.Namespace(command=None))
    with pytest.raises(SystemExit):
        seqsender_module.main()
    assert parser.help_printed is True
