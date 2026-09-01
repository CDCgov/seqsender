# tests/unit/test_file_handler.py
import pandas as pd
import pytest
import builtins
import importlib.util
import sys
import types
import os
from Bio.Seq import Seq
from typing import Any
from pathlib import Path

def _source_root() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        # During mutmut, prefer the mutated source tree.
        mutant_src = parent / "src"
        if parent.name == "mutants" and (mutant_src / "file_handler.py").exists():
            return mutant_src

        # Normal pytest run.
        normal_src = parent / "src"
        if (normal_src / "file_handler.py").exists() and parent.name != "mutants":
            return normal_src

    raise RuntimeError("Could not find src/file_handler.py")

SOURCE_DIR = _source_root()
MODULE_PATH = SOURCE_DIR / "file_handler.py"
if str(SOURCE_DIR) not in sys.path:
    sys.path.insert(0, str(SOURCE_DIR))

"""Unit tests initially generated with ChatGPT v5.5 "Deep Research" and "Thinking"
then modified for use and clarity with validation and coverage testing via mutmut."""

#*******************************************************************************
#                 create test file_handler.py connections
#*******************************************************************************

@pytest.fixture()
def file_handler_module(monkeypatch):
    monkeypatch.syspath_prepend(str(SOURCE_DIR))
    src_pkg: Any = types.ModuleType("src")
    src_pkg.__path__ = [str(SOURCE_DIR / "src")]
    monkeypatch.setitem(sys.modules, "src", src_pkg)

    def alias_src_module(name, module):
        monkeypatch.setitem(sys.modules, f"src.{name}", module)
        monkeypatch.setitem(sys.modules, name, module)
        setattr(src_pkg, name, module)

    settings_stub: Any = types.ModuleType("settings")
    settings_stub.SAMPLE_NAME_DATABASE_PREFIX = {
        "BIOSAMPLE": "bs-",
        "SRA": "sra-",
        "GENBANK": "gb-",
    }
    settings_stub.PROG_DIR = str(SOURCE_DIR)

    for module_name in [
        "file_handler",
        "src.file_handler",
        "settings",
        "src.settings",
    ]:
        sys.modules.pop(module_name, None)

    monkeypatch.setitem(sys.modules, "src", src_pkg)
    alias_src_module("src.settings", settings_stub)
    spec = importlib.util.spec_from_file_location("file_handler", MODULE_PATH)
    assert spec and spec.loader
    module: Any = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    setattr(src_pkg, "file_handler", module)
    spec.loader.exec_module(module)
    return module

#*******************************************************************************
#                               copy_file
#*******************************************************************************

def test_copy_file__copies_contents(file_handler_module, tmp_path):
    file_handler = file_handler_module
    source = tmp_path / "source.txt"
    destination = tmp_path / "destination.txt"
    source.write_text("hello", encoding="utf-8")
    assert source.exists()
    file_handler.copy_file(str(source), str(destination))
    assert destination.exists()
    assert destination.read_text(encoding="utf-8") == "hello"

#*******************************************************************************
#                               validate_file
#*******************************************************************************

def test_validate_file__passes_for_existing_file(file_handler_module, tmp_path):
    file_handler = file_handler_module
    existing = tmp_path / "metadata.csv"
    existing.write_text("a\n1\n", encoding="utf-8")

    file_handler.validate_file("metadata_file", str(existing))

def test_validate_file__exits_for_missing_file(file_handler_module, tmp_path, capsys):
    file_handler = file_handler_module
    missing = tmp_path / "missing.csv"

    with pytest.raises(SystemExit) as exc:
        file_handler.validate_file("metadata_file", str(missing))

    assert exc.value.code == 1
    assert f"Error: Input metadata file does not exist at: {missing}\n" == capsys.readouterr().err
#*******************************************************************************
#                             validate_directory
#*******************************************************************************

def test_validate_directory__passes_for_existing_directory(file_handler_module, tmp_path):
    file_handler = file_handler_module
    try:
        file_handler.validate_directory("submission directory", str(tmp_path))
    except SystemExit:
        pytest.fail("Function unexpectedly called sys.exit(1)")

def test_validate_directory__exits_for_missing_directory(file_handler_module, tmp_path, capsys):
    file_handler = file_handler_module
    missing = tmp_path / "missing-dir"

    with pytest.raises(SystemExit) as exc:
        file_handler.validate_directory("submission directory", str(missing))

    assert f"There is no submission directory at: {missing}\n" == capsys.readouterr().err
    assert exc.value.code == 1

#*******************************************************************************
#                            create_directory
#*******************************************************************************

def test_create_directory__creates_nested_directory_and_is_idempotent(file_handler_module, tmp_path):
    file_handler = file_handler_module
    directory = tmp_path / "nested" / "output"

    file_handler.create_directory(str(directory))
    file_handler.create_directory(str(directory))

    assert directory.is_dir()

#*******************************************************************************
#                                load_yaml
#*******************************************************************************

def test_load_yaml__returns_parsed_yaml(file_handler_module, tmp_path):
    file_handler = file_handler_module
    yaml_file = tmp_path / "config.yaml"
    yaml_file.write_text("Submission:\n  NCBI:\n    Username: user\n", encoding="utf-8")

    result = file_handler.load_yaml("Config file", str(yaml_file))

    assert result == {"Submission": {"NCBI": {"Username": "user"}}}


def test_load_yaml__exits_for_invalid_yaml(file_handler_module, tmp_path, capsys):
    file_handler = file_handler_module
    yaml_file = tmp_path / "bad.yaml"
    yaml_file.write_text("Submission: [unterminated\n", encoding="utf-8")

    with pytest.raises(SystemExit) as exc:
        file_handler.load_yaml("Config file", str(yaml_file))

    assert exc.value.code == 1
    assert "Error: Config file is incorrect. File must be a valid yaml format.\n" == capsys.readouterr().err

def test_load_yaml__opens_file_in_explicit_read_mode(file_handler_module, tmp_path, monkeypatch):
    yaml_file = tmp_path / "config.yaml"
    yaml_file.write_text("a: 1\n", encoding="utf-8")
    original_open = open
    observed_open_calls = []

    def fake_open(file, *args, **kwargs):
        observed_open_calls.append((str(file), args, kwargs))
        return original_open(file, *args, **kwargs)

    monkeypatch.setattr("builtins.open", fake_open)
    assert file_handler_module.load_yaml("Config file", str(yaml_file)) == {"a": 1}
    assert observed_open_calls == [(str(yaml_file), ("r",), {})]

#*******************************************************************************
#                                is_row_empty
#*******************************************************************************

@pytest.mark.parametrize(
    ("row", "expected"),
    [
        (pd.Series([None, "", "   "]), True),
        (pd.Series([None, "  x  ", ""]), False),
        (pd.Series([0, "", None]), False),
    ],
)
def test_is_row_empty__classifies_empty_rows(file_handler_module, row, expected):
    file_handler = file_handler_module
    assert file_handler.is_row_empty(row) is expected

#*******************************************************************************
#                                load_csv
#*******************************************************************************

def test_load_csv__strips_column_names_and_drops_wholly_empty_rows(file_handler_module, tmp_path):
    file_handler = file_handler_module
    csv_file = tmp_path / "metadata.csv"
    csv_file.write_text(" col_a , col_b \n1,x\n   ,   \n2,y\n", encoding="utf-8")

    df = file_handler.load_csv(str(csv_file))

    assert list(df.columns) == ["col_a", "col_b"]
    assert df.shape == (2, 2)
    assert df["col_a"].tolist() == ["1", "2"]

def test_load_csv__default_separator_is_comma(file_handler_module, tmp_path):
    csv_file = tmp_path / "metadata.csv"
    csv_file.write_text("a,b\n1,2\n", encoding="utf-8")

    df = file_handler_module.load_csv(str(csv_file))

    assert list(df.columns) == ["a", "b"]
    assert df.to_dict(orient="records") == [{"a": "1", "b": "2"}]

def test_load_csv__supports_custom_separator(file_handler_module, tmp_path):
    file_handler = file_handler_module
    tsv_file = tmp_path / "metadata.tsv"
    tsv_file.write_text("a\tb\n1\t2\n", encoding="utf-8")

    df = file_handler.load_csv(str(tsv_file), sep="\t")

    assert df.to_dict(orient="records") == [{"a": "1", "b": "2"}]

def test_load_csv__passes_exact_read_csv_settings(file_handler_module, monkeypatch, tmp_path):
    file_handler = file_handler_module
    csv_file = tmp_path / "metadata.csv"
    csv_file.write_text(" a , b \n1,2\n", encoding="utf-8")

    observed_calls: list[dict[str, Any]] = []

    def fake_read_csv(*args, **kwargs):
        observed_calls.append({"args": args, "kwargs": kwargs})
        return pd.DataFrame([{" a ": "1", " b ": "2"},{" a ": "   ", " b ": "   "}])

    monkeypatch.setattr(file_handler.pd, "read_csv", fake_read_csv)

    df = file_handler.load_csv(str(csv_file), sep="|")
    # pd.read_csv specific settings
    assert observed_calls == [
        {
            "args": (str(csv_file),),
            "kwargs": {
                "header": 0,
                "dtype": str,
                "sep": "|",
                "engine": "python",
                "encoding": "utf-8",
                "index_col": False,
                "na_filter": False,
            },
        }
    ]

    assert list(df.columns) == ["a", "b"]
    assert df.to_dict("records") == [{"a": "1", "b": "2"}]

#*******************************************************************************
#                                load_fasta_file
#*******************************************************************************

def test_load_fasta_file__returns_expected_columns(file_handler_module, tmp_path):
    file_handler = file_handler_module
    fasta_file = tmp_path / "sequence.fasta"
    fasta_file.write_text(">seq1 description one\nACTG\n>seq2\nTTAA\n", encoding="utf-8")

    df = file_handler.load_fasta_file(str(fasta_file))

    assert list(df.columns) == ["fasta_name_orig", "fasta_sequence_orig", "fasta_description_orig"]
    assert df["fasta_name_orig"].tolist() == ["seq1", "seq2"]
    assert [str(seq) for seq in df["fasta_sequence_orig"]] == ["ACTG", "TTAA"]
    assert df.loc[0, "fasta_description_orig"] == "seq1 description one"

def test_load_fasta_file__passes_read_mode_explicitly(file_handler_module, tmp_path, monkeypatch):
    fasta_file = tmp_path / "sequence.fasta"
    fasta_file.write_text(">seq1 description one\nACTG\n", encoding="utf-8")

    original_open = open
    observed_open_calls = []

    def fake_open(file, *args, **kwargs):
        observed_open_calls.append((str(file), args, kwargs))
        return original_open(file, *args, **kwargs)

    monkeypatch.setattr("builtins.open", fake_open)

    df = file_handler_module.load_fasta_file(str(fasta_file))

    assert observed_open_calls == [(str(fasta_file), ("r",), {})]
    assert list(df.columns) == [
        "fasta_name_orig",
        "fasta_sequence_orig",
        "fasta_description_orig",
    ]
    assert df["fasta_name_orig"].tolist() == ["seq1"]
    assert [str(seq) for seq in df["fasta_sequence_orig"]] == ["ACTG"]
    assert df.loc[0, "fasta_description_orig"] == "seq1 description one"
#*******************************************************************************
#                                save_xml
#*******************************************************************************

def test_save_xml__writes_submission_xml(file_handler_module, tmp_path):
    file_handler = file_handler_module
    file_handler.save_xml(b"<Submission />", str(tmp_path))

    assert (tmp_path / "submission.xml").read_bytes() == b"<Submission />"


def test_save_xml__exits_on_permission_error(file_handler_module, tmp_path, monkeypatch, capsys):
    file_handler = file_handler_module

    def fake_open(*_args, **_kwargs):
        raise PermissionError("denied")

    monkeypatch.setattr(builtins, "open", fake_open)

    with pytest.raises(SystemExit) as exc:
        file_handler.save_xml(b"<Submission />", str(tmp_path))

    assert exc.value.code == 1
    assert f"Error: Permission error when trying to save 'submission.xml' to path: {tmp_path}\ndenied\n" == capsys.readouterr().err

def test_save_xml__exits_on_unexpected_error(file_handler_module, tmp_path, monkeypatch, capsys):
    file_handler = file_handler_module
    def fake_open(*_args, **_kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(builtins, "open", fake_open)

    with pytest.raises(SystemExit) as exc:
        file_handler.save_xml(b"<Submission />", str(tmp_path))

    assert exc.value.code == 1
    assert f"Error: An unexpected error occurred when trying to save 'submission.xml' to path: {tmp_path}\ndisk full\n" == capsys.readouterr().err

def test_save_xml__checks_expected_submission_xml_path(file_handler_module, tmp_path, monkeypatch):
    expected_path = str(tmp_path / "submission.xml")
    observed_exists_paths = []
    original_exists = os.path.exists

    def fake_exists(path):
        observed_exists_paths.append(path)
        if path == expected_path:
            return True
        return False

    def fake_sleep(seconds):
        pytest.fail(f"save_xml unexpectedly slept for {seconds} seconds")

    monkeypatch.setattr(file_handler_module.os.path, "exists", fake_exists)
    monkeypatch.setattr(file_handler_module.time, "sleep", fake_sleep)
    file_handler_module.save_xml(b"<Submission />", str(tmp_path))
    assert observed_exists_paths == [expected_path]
    assert (tmp_path / "submission.xml").read_bytes() == b"<Submission />"

def test_save_xml__wait_loop_sleeps_for_10_seconds_until_file_exists(file_handler_module, tmp_path, monkeypatch):
    expected_path = str(tmp_path / "submission.xml")
    exists_calls = 0
    sleep_calls = []

    def fake_exists(path):
        nonlocal exists_calls
        assert path == expected_path
        exists_calls += 1
        return exists_calls >= 2

    def fake_sleep(seconds):
        sleep_calls.append(seconds)

    monkeypatch.setattr(file_handler_module.os.path, "exists", fake_exists)
    monkeypatch.setattr(file_handler_module.time, "sleep", fake_sleep)
    file_handler_module.save_xml(b"<Submission />", str(tmp_path))
    assert sleep_calls == [10]

#*******************************************************************************
#                                save_csv
#*******************************************************************************

def test_save_csv__writes_csv_to_explicit_file_path(file_handler_module, tmp_path):
    file_handler = file_handler_module
    output = tmp_path / "metadata.csv"
    df = pd.DataFrame([{"a": "1", "b": "2"}])
    file_handler.save_csv(df, str(output))
    assert output.read_text(encoding="utf-8") == "a,b\n1,2\n"

def test_save_csv__passes_exact_to_csv_settings_by_default(file_handler_module, tmp_path, monkeypatch):
    observed_calls = []

    def fake_to_csv(self, *args, **kwargs):
        observed_calls.append({"args": args, "kwargs": kwargs})

    monkeypatch.setattr(pd.DataFrame, "to_csv", fake_to_csv)
    df = pd.DataFrame([{"a": "1", "b": "2"}])
    output = tmp_path / "metadata.csv"
    file_handler_module.save_csv(df, str(output))
    assert observed_calls == [{"args": (str(output),), "kwargs": {"header": True, "index": False, "sep": ","}}]

def test_save_csv__writes_csv_with_file_name_and_separator(file_handler_module, tmp_path):
    file_handler = file_handler_module
    df = pd.DataFrame([{"a": "1", "b": "2"}])
    file_handler.save_csv(df, str(tmp_path), file_name="metadata.tsv", sep="\t")
    assert (tmp_path / "metadata.tsv").read_text(encoding="utf-8") == "a\tb\n1\t2\n"


def test_save_csv__exits_on_permission_error(file_handler_module, monkeypatch, tmp_path, capsys):
    file_handler = file_handler_module
    def fake_to_csv(self, *_args, **_kwargs):
        raise PermissionError("denied")

    monkeypatch.setattr(pd.DataFrame, "to_csv", fake_to_csv)

    with pytest.raises(SystemExit) as exc:
        file_handler.save_csv(pd.DataFrame([{"a": 1}]), str(tmp_path), file_name="out.csv")

    assert exc.value.code == 1
    assert f"Error: Permission error when trying to save 'out.csv' to path: {tmp_path / 'out.csv'}\ndenied\n" == capsys.readouterr().err

def test_save_csv__exits_on_unexpected_error(file_handler_module, monkeypatch, tmp_path, capsys):
    file_handler = file_handler_module
    def fake_to_csv(self, *_args, **_kwargs):
        raise RuntimeError("boom")

    monkeypatch.setattr(pd.DataFrame, "to_csv", fake_to_csv)

    with pytest.raises(SystemExit) as exc:
        file_handler.save_csv(pd.DataFrame([{"a": 1}]), str(tmp_path), file_name="out.csv")

    assert exc.value.code == 1
    assert f"Error: An unexpected error occurred when trying to save 'out.csv' to path: {tmp_path / 'out.csv'}\nboom\n" == capsys.readouterr().err

#*******************************************************************************
#                                create_fasta
#*******************************************************************************

def test_create_fasta__writes_genbank_headers_with_bioproject_and_modifiers(file_handler_module, tmp_path):
    file_handler = file_handler_module
    metadata = pd.DataFrame(
        [
            {
                "gb-sample_name": "GB1",
                "bioproject": "PRJNA123",
                "gb-fasta_definition_line_modifiers": "[country=USA]",
                "fasta_sequence_orig": Seq("ACTG"),
            }
        ]
    )

    file_handler.create_fasta(
        database="GENBANK",
        metadata=metadata,
        submission_dir=str(tmp_path),
        config_dict={"Add_Definition_Line_Accessions": True},
    )

    content = (tmp_path / "sequence.fsa").read_text(encoding="utf-8")
    assert ">GB1 [BioProject=PRJNA123] [country=USA]" in content
    assert "ACTG" in content

def test_create_fasta__omits_blank_bioproject_even_when_flag_true(file_handler_module, tmp_path):
    metadata = pd.DataFrame(
        [
            {
                "gb-sample_name": "GB1",
                "bioproject": "   ",
                "gb-fasta_definition_line_modifiers": "",
                "fasta_sequence_orig": Seq("ACTG"),
            }
        ]
    )

    file_handler_module.create_fasta(
        database="GENBANK",
        metadata=metadata,
        submission_dir=str(tmp_path),
        config_dict={"Add_Definition_Line_Accessions": True},
    )

    content = (tmp_path / "sequence.fsa").read_text(encoding="utf-8")
    assert ">GB1" in content
    assert "BioProject" not in content
    assert "XXXX" not in content

def test_create_fasta__omits_blank_genbank_definition_line_modifiers(file_handler_module, tmp_path):
    metadata = pd.DataFrame(
        [
            {
                "gb-sample_name": "GB1",
                "bioproject": "",
                "gb-fasta_definition_line_modifiers": "   ",
                "fasta_sequence_orig": Seq("ACTG"),
            }
        ]
    )

    file_handler_module.create_fasta(
        database="GENBANK",
        metadata=metadata,
        submission_dir=str(tmp_path),
        config_dict={"Add_Definition_Line_Accessions": True},
    )

    content = (tmp_path / "sequence.fsa").read_text(encoding="utf-8")
    assert ">GB1" in content
    assert "[BioProject=" not in content
    assert "XXXX" not in content
    assert "[country=" not in content

def test_create_fasta__creates_seqrecords_with_empty_description(file_handler_module, tmp_path, monkeypatch):
    observed_records = []

    class FakeSeqRecord:
        def __init__(self, seq, id, description="DEFAULT"):
            self.seq = seq
            self.id = id
            self.description = description
            observed_records.append(
                {"seq": seq, "id": id, "description": description}
            )

    def fake_write(records, file_obj, fmt):
        file_obj.write("written\n")
        return len(records)

    monkeypatch.setattr(file_handler_module, "SeqRecord", FakeSeqRecord)
    monkeypatch.setattr(file_handler_module.SeqIO, "write", fake_write)
    metadata = pd.DataFrame(
        [
            {
                "gb-sample_name": "GB1",
                "bioproject": "",
                "gb-fasta_definition_line_modifiers": "",
                "fasta_sequence_orig": Seq("ACTG"),
            }
        ]
    )
    file_handler_module.create_fasta("GENBANK", metadata, str(tmp_path), {})
    assert observed_records == [{"seq": Seq("ACTG"), "id": "GB1", "description": ""}]

def test_create_fasta__exits_on_permission_error(file_handler_module, tmp_path, monkeypatch, capsys):
    file_handler = file_handler_module
    metadata = pd.DataFrame([{"gb-sample_name": "GB1", "fasta_sequence_orig": Seq("ACTG")}])

    def fake_open(*_args, **_kwargs):
        raise PermissionError("denied")

    monkeypatch.setattr(builtins, "open", fake_open)

    with pytest.raises(SystemExit) as exc:
        file_handler.create_fasta("GENBANK", metadata, str(tmp_path), {})

    assert exc.value.code == 1
    assert f"Error: Permission error when trying to save 'sequence.fsa' to path: {tmp_path}\ndenied\n" == capsys.readouterr().err

def test_create_fasta__exits_on_unexpected_error(file_handler_module, tmp_path, monkeypatch, capsys):
    file_handler = file_handler_module
    metadata = pd.DataFrame([{"gb-sample_name": "GB1", "fasta_sequence_orig": Seq("ACTG")}])

    def fake_open(*_args, **_kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(builtins, "open", fake_open)

    with pytest.raises(SystemExit) as exc:
        file_handler.create_fasta("GENBANK", metadata, str(tmp_path), {})

    assert exc.value.code == 1
    assert f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {tmp_path}\ndisk full\n" == capsys.readouterr().err
