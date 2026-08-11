from __future__ import annotations

import argparse
import importlib.util
import sys
import types
from typing import Any
import os
from pathlib import Path

import pytest

def _source_root() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        # During mutmut, prefer the mutated source tree.
        mutant_src = parent / "src"
        if parent.name == "mutants" and (mutant_src / "argument_handler.py").exists():
            return mutant_src

        # Normal pytest run.
        normal_src = parent / "src"
        if (normal_src / "argument_handler.py").exists() and parent.name != "mutants":
            return normal_src

    raise RuntimeError("Could not find src/argument_handler.py")


SOURCE_DIR = _source_root()

if str(SOURCE_DIR) not in sys.path:
    sys.path.insert(0, str(SOURCE_DIR))

MODULE_PATH = SOURCE_DIR / "argument_handler.py"

settings_stub: Any = types.ModuleType("settings")
settings_stub.ORGANISM_CHOICES = ["FLU", "COV", "POX", "ARBO", "RSV", "OTHER"]
src_pkg: Any = types.ModuleType("src")
src_pkg.__path__ = [str(SOURCE_DIR)]

sys.modules.pop("argument_handler", None)
sys.modules.pop("src.argument_handler", None)
sys.modules.pop("settings", None)
sys.modules.pop("src.settings", None)

sys.modules["src"] = src_pkg
sys.modules["settings"] = settings_stub
sys.modules["src.settings"] = settings_stub
setattr(src_pkg, "settings", settings_stub)

spec = importlib.util.spec_from_file_location("argument_handler", MODULE_PATH)
assert spec and spec.loader

argument_handler: Any = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = argument_handler
setattr(src_pkg, "argument_handler", argument_handler)

spec.loader.exec_module(argument_handler)

def parse(argv: list[str]) -> argparse.Namespace:
    return argument_handler.args_parser().parse_args(argv)

def _option(parser: argparse.ArgumentParser, dest: str) -> argparse.Action:
    matches = [action for action in parser._actions if action.dest == dest]
    assert len(matches) == 1
    return matches[0]

def _subparsers(parser: argparse.ArgumentParser) -> argparse._SubParsersAction:
    matches = [action for action in parser._actions if isinstance(action, argparse._SubParsersAction)]
    assert len(matches) == 1
    return matches[0]

def test_prep_requires_config_file(capsys: pytest.CaptureFixture[str]):
    with pytest.raises(SystemExit) as exc:
        parse(["prep", "--biosample", "--organism", "FLU", "--submission_name", "sub1", "--submission_dir", "/tmp/out", "--metadata_file", "metadata.csv"])
    assert exc.value.code == 2
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "the following arguments are required: --config_file" in captured.err

def test_args_parser__top_level_parser_metadata_is_exact():
    parser = argument_handler.args_parser()
    assert parser.description == "Genomic tool to simplify/automate the process of submitting organism samples to public repositories. With built-in tools to create/submit/link/log organism samples for the databases: BioSample, SRA, GenBank, and GISAID."
    assert parser.formatter_class is argparse.ArgumentDefaultsHelpFormatter



def test_args_parser__prep_shared_options_have_exact_metadata():
    prep = _subparsers(argument_handler.args_parser()).choices["prep"]

    expected = {
        "biosample": {
            "option_strings": ["--biosample", "-b"],
            "help": "Create/Submit BioSample data.",
            "const": "BIOSAMPLE",
            "default": "",
            "required": False,
        },
        "sra": {
            "option_strings": ["--sra", "-s"],
            "help": "Create/Submit SRA data.",
            "const": "SRA",
            "default": "",
            "required": False,
        },
        "genbank": {
            "option_strings": ["--genbank", "-n"],
            "help": "Create/Submit GenBank data. (requires --fasta_file)",
            "const": "GENBANK",
            "default": "",
            "required": False,
        },
        "gisaid": {
            "option_strings": ["--gisaid", "-g"],
            "help": "Create/Submit GISAID data. (requires --fasta_file)",
            "const": "GISAID",
            "default": "",
            "required": False,
        },
        "organism": {
            "option_strings": ["--organism"],
            "help": "Type of organism data. Listed organism options have unique submissions options/processes, if your specific organism is not listed, use 'OTHER' for options available to all organisms.",
            "default": "",
            "required": True,
            "choices": ["FLU", "COV", "POX", "ARBO", "RSV", "OTHER"],
        },
        "submission_name": {
            "option_strings": ["--submission_name"],
            "help": "Unique name for the submission of your data. Reusing the same name can cause issues during the submission process. A folder will be created at: 'submission_dir/submission_name'.",
            "required": True,
        },
        "submission_dir": {
            "option_strings": ["--submission_dir"],
            "help": (
                "Output directory where all files for your submission will be stored. A folder will be created at '<submission_dir>/<submission_name>'; this is the location where: all of the submission files will be created, SeqSender "
                "will stage each step of the submission process automatically, and where SeqSender will generate all the output from your submission."
            ),
            "required": True,
        },
        "config_file": {
            "option_strings": ["--config_file"],
            "help": (
                "Config file to be used in the creation/submission of your samples. SeqSender will store this file location in your 'submission_log.csv' where it will use it to manage your submission, be careful when modifying and "
                "ensure SeqSender maintains access to this file. Input either full file path or if just file name it must be stored at '<submission_dir>/<submission_name>/<config_file>'."
            ),
            "required": True,
        },
        "metadata_file": {
            "option_strings": ["--metadata_file"],
            "help": "Metadata file to be used in the creation/submission of your samples. Input either full file path or if just file name it must be stored at '<submission_dir>/<submission_name>/<metadata_file>'.",
            "required": True,
        },
    }

    for dest, expected_values in expected.items():
        action = _option(prep, dest)
        for attr, value in expected_values.items():
            assert getattr(action, attr) == value

def test_args_parser__prep_optional_flags_have_exact_metadata():
    prep = _subparsers(argument_handler.args_parser()).choices["prep"]

    expected = {
        "skip_validation": {
            "option_strings": ["--skip_validation"],
            "help": (
                "Skip initial validation for metadata file. Validation will still occur for the 'config_file' and for any subsequent submissions made via "
                "'submission_status'. Warning, this can cause unexpected errors using SeqSender if required columns are missing."
            ),
            "required": False,
            "const": True,
            "default": False,
        },
        "fasta_file": {
            "option_strings": ["--fasta_file"],
            "help": (
                "Fasta file used to generate submission files; fasta header should match the column 'sequence_name' stored in your metadata. Input either full file "
                "path or if just file name it must be stored at '<submission_dir>/<submission_name>/<fasta_file>'."
            ),
            "default": None,
            "required": False,
        },
        "table2asn": {
            "option_strings": ["--table2asn"],
            "help": "Perform a table2asn submission instead of GenBank FTP submission for organism choices 'FLU' or 'COV'.",
            "required": False,
            "const": True,
            "default": False,
        },
        "gff_file": {
            "option_strings": ["--gff_file"],
            "help": "Annotation file only available for table2asn submissions. (requires '--table2asn' for organism choices 'FLU', or 'COV').",
            "default": None,
            "required": False,
        },
        "publication_title": {
            "option_strings": ["--publication_title"],
            "help": "Publication Title associated with sample submission. For GenBank only, overwrites value given via config file.",
            "required": False,
            "default": None,
        },
        "publication_status": {
            "option_strings": ["--publication_status"],
            "help": "Status of publication associated with sample submission. For GenBank only, overwrites value given via config file.",
            "required": False,
            "default": None,
            "choices": ["Unpublished", "In-press", "Published"],
        },
    }
    for dest, expected_values in expected.items():
        action = _option(prep, dest)
        for attr, value in expected_values.items():
            assert getattr(action, attr) == value

def test_args_parser__submit_test_flag_has_exact_metadata():
    submit = _subparsers(argument_handler.args_parser()).choices["submit"]
    action = _option(submit, "test")
    assert action.option_strings == ["--test"]
    assert action.help == "Perform a test submission."
    assert action.const is True
    assert action.default is False
    assert action.required is False

def test_args_parser__submission_status_submission_name_is_optional_with_exact_help():
    status = _subparsers(argument_handler.args_parser()).choices["submission_status"]
    action = _option(status, "submission_name")
    assert action.option_strings == ["--submission_name"]
    assert action.required is False
    assert action.default is None
    assert action.help == "Unique name for the submission of your data. This is an optional field if you want Seqsender to only update the specified submission in the 'submission_log.csv'."

# Test main commands
@pytest.mark.parametrize("command", ["prep", "submit"])
# Test organism flags
@pytest.mark.parametrize("organism", ["FLU", "COV", "POX", "ARBO", "RSV", "OTHER"])
# Test all databases flags
@pytest.mark.parametrize(("bs_flag", "bs_value"), [("", ""), ("--biosample", "BIOSAMPLE"), ("-b", "BIOSAMPLE")])
@pytest.mark.parametrize(("sra_flag", "sra_value"), [("", ""), ("--sra", "SRA"), ("-s", "SRA")])
@pytest.mark.parametrize(("gb_flag", "gb_value"), [("", ""), ("--genbank", "GENBANK"), ("-n", "GENBANK")])
@pytest.mark.parametrize(("gs_flag", "gs_value"), [("", ""), ("--gisaid", "GISAID"), ("-g", "GISAID")])
def test_prep_and_submit_database_flags(command, organism,
    bs_flag, bs_value, sra_flag, sra_value, gb_flag, gb_value, gs_flag, gs_value):
    args = parse(
        [
            arg
            for arg in [
                command,
                bs_flag,
                sra_flag,
                gb_flag,
                gs_flag,
                "--organism",
                organism,
                "--submission_name",
                "sub1",
                "--submission_dir",
                "/tmp/out",
                "--config_file",
                "config.yaml",
                "--metadata_file",
                "metadata.csv",
            ]
            if arg
        ],
    )
    if bs_value:
        assert args.biosample == "BIOSAMPLE"
    if sra_value:
        assert args.sra == "SRA"
    if gb_value:
        assert args.genbank == "GENBANK"
    if gs_value:
        assert args.gisaid == "GISAID"
    assert args.command == command
    assert args.organism == organism


def test_prep_defaults():
    args = parse(
        [
            "prep",
            "--biosample",
            "--organism",
            "COV",
            "--submission_name",
            "sub1",
            "--submission_dir",
            "/tmp/out",
            "--config_file",
            "config.yaml",
            "--metadata_file",
            "metadata.csv",
        ],
    )
    assert args.skip_validation is False
    assert args.table2asn is False
    assert args.fasta_file is None
    assert args.gff_file is None
    assert args.publication_title is None
    assert args.publication_status is None


def test_prep_optional_file_and_publication_flags():
    args = parse(
        [
            "prep",
            "--genbank",
            "--organism",
            "FLU",
            "--submission_name",
            "sub1",
            "--submission_dir",
            "/tmp/out",
            "--config_file",
            "config.yaml",
            "--metadata_file",
            "metadata.csv",
            "--fasta_file",
            "seqs.fasta",
            "--gff_file",
            "ann.gff",
            "--table2asn",
            "--skip_validation",
            "--publication_title",
            "A title",
            "--publication_status",
            "Published",
        ],
    )
    assert args.fasta_file == "seqs.fasta"
    assert args.gff_file == "ann.gff"
    assert args.table2asn is True
    assert args.skip_validation is True
    assert args.publication_title == "A title"
    assert args.publication_status == "Published"


def test_submit_accepts_test_flag():
    args = parse(
        [
            "submit",
            "--sra",
            "--organism",
            "OTHER",
            "--submission_name",
            "sub1",
            "--submission_dir",
            "/tmp/out",
            "--config_file",
            "config.yaml",
            "--metadata_file",
            "metadata.csv",
            "--test",
        ],
    )
    assert args.command == "submit"
    assert args.test is True

def test_submission_status_command_has_optional_submission_name():
    args = parse(["submission_status", "--submission_dir", "/tmp/out"])
    assert args.command == "submission_status"
    assert args.submission_name is None

    args = parse(["submission_status", "--submission_dir", "/tmp/out", "--submission_name", "sub1"])
    assert args.submission_name == "sub1"

def test_generate_test_data_command():
    args = parse(["test_data", "--gisaid", "--organism", "FLU", "--submission_dir", "/tmp/out"])
    assert args.command == "test_data"
    assert args.gisaid == "GISAID"

@pytest.mark.parametrize("command", ["test_network_connection"])
def test_miscellaneous_commands(command):
    args = parse([command])
    assert args.command == command

@pytest.mark.parametrize("command", ["prep", "submit", "test_data"])
def test_organism_is_required_for_database_submission_commands(command):
    argv = [command, "--biosample", "--submission_dir", "/tmp/out"]
    if command != "test_data":
        argv += ["--submission_name", "sub1", "--config_file", "config.yaml", "--metadata_file", "metadata.csv"]
    with pytest.raises(SystemExit):
        parse(argv)

def test_organism_choices_are_enforced():
    with pytest.raises(SystemExit):
        parse(
            [
                "prep",
                "--biosample",
                "--organism",
                "BAD",
                "--submission_name",
                "sub1",
                "--submission_dir",
                "/tmp/out",
                "--config_file",
                "config.yaml",
                "--metadata_file",
                "metadata.csv",
            ],
        )

def test_publication_status_choices_are_enforced():
    with pytest.raises(SystemExit):
        parse(
            [
                "prep",
                "--genbank",
                "--organism",
                "FLU",
                "--submission_name",
                "sub1",
                "--submission_dir",
                "/tmp/out",
                "--config_file",
                "config.yaml",
                "--metadata_file",
                "metadata.csv",
                "--fasta_file",
                "seqs.fasta",
                "--publication_status",
                "Draft",
            ],
        )


def test_subparser_defaults():
    args = parse(
        [
            "prep",
            "--biosample",
            "--sra",
            "--organism",
            "FLU",
            "--submission_name",
            "sub1",
            "--submission_dir",
            "/tmp/out",
            "--config_file",
            "config.yaml",
            "--metadata_file",
            "metadata.csv",
        ],
    )
    assert args.biosample == "BIOSAMPLE"
    assert args.sra == "SRA"
    assert args.gisaid == ""
    assert args.genbank == ""
    assert args.organism == "FLU"
    assert args.skip_validation is False
    assert args.config_file == "config.yaml"
    assert args.metadata_file == "metadata.csv"
    assert args.fasta_file is None
    assert args.table2asn is False
    assert args.gff_file is None
    assert args.publication_title is None
    assert args.publication_status is None

def test_subparser_defaults_fasta():
    args = parse(
        [
            "prep",
            "--genbank",
            "--gisaid",
            "--organism",
            "COV",
            "--submission_name",
            "sub1",
            "--submission_dir",
            "/tmp/out",
            "--config_file",
            "config.yaml",
            "--metadata_file",
            "metadata.csv",
            "--fasta_file",
            "seqs.fasta",
        ],
    )
    assert args.biosample == ""
    assert args.sra == ""
    assert args.gisaid == "GISAID"
    assert args.genbank == "GENBANK"
    assert args.organism == "COV"
    assert args.skip_validation is False
    assert args.fasta_file == "seqs.fasta"
    assert args.config_file == "config.yaml"
    assert args.metadata_file == "metadata.csv"
    assert args.table2asn is False
    assert args.gff_file is None
    assert args.publication_title is None
    assert args.publication_status is None

def test_version_command_parses_exactly():
    args = parse(["version"])
    assert args.command == "version"
