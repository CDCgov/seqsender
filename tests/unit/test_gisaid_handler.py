from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path
import os
import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import Any

def _source_root() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        # During mutmut, prefer the mutated source tree.
        mutant_src = parent / "src"
        if parent.name == "mutants" and (mutant_src / "gisaid_handler.py").exists():
            return mutant_src

        # Normal pytest run.
        normal_src = parent / "src"
        if (normal_src / "gisaid_handler.py").exists() and parent.name != "mutants":
            return normal_src

    raise RuntimeError("Could not find src/gisaid_handler.py")


SOURCE_DIR = _source_root()
MODULE_PATH = SOURCE_DIR / "gisaid_handler.py"
if str(SOURCE_DIR) not in sys.path:
    sys.path.insert(0, str(SOURCE_DIR))

"""Unit tests initially generated with ChatGPT v5.5 "Deep Research" and "Thinking"
then modified for use and clarity with validation and coverage testing via mutmut."""

#*******************************************************************************
#                     create test gisaid.py connections
#*******************************************************************************

@pytest.fixture()
def gisaid_handler_module(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.syspath_prepend(str(SOURCE_DIR))
    src_pkg: Any = types.ModuleType("src")
    src_pkg.__path__ = [str(SOURCE_DIR)]
    monkeypatch.setitem(sys.modules, "src", src_pkg)

    def alias_src_module(name, module):
        monkeypatch.setitem(sys.modules, f"src.{name}", module)
        monkeypatch.setitem(sys.modules, name, module)
        setattr(src_pkg, name, module)

    upload_log_stub: Any = types.ModuleType("upload_log")
    upload_log_stub.update_calls = []

    def update_submission_status_csv(**kwargs):
        upload_log_stub.update_calls.append(kwargs)

    upload_log_stub.update_submission_status_csv = update_submission_status_csv

    tools_stub: Any = types.ModuleType("tools")
    tools_stub.check_credentials_calls = []

    def check_credentials(**kwargs):
        tools_stub.check_credentials_calls.append(kwargs)

    tools_stub.check_credentials = check_credentials

    file_handler_stub: Any = types.ModuleType("file_handler")
    file_handler_stub.saved_csvs = []
    file_handler_stub.created_fastas = []
    file_handler_stub.validate_file_calls = []
    file_handler_stub.validate_gisaid_installer_calls = []

    def save_csv(**kwargs):
        file_handler_stub.saved_csvs.append(kwargs)
        df = kwargs["df"]
        path = Path(kwargs["file_path"]) / kwargs.get("file_name", "")
        if kwargs.get("file_name"):
            df.to_csv(path, index=False, sep=kwargs.get("sep", ","))

    def create_fasta(**kwargs):
        file_handler_stub.created_fastas.append(kwargs)
        submission_dir = Path(kwargs["submission_dir"])
        submission_dir.mkdir(parents=True, exist_ok=True)
        records = []
        metadata = kwargs["metadata"]
        for _, row in metadata.iterrows():
            records.append(SeqRecord(Seq("ACGT"), id=row.get("gs-sample_name", "sample"), description=""))
        with open(submission_dir / "sequence.fsa", "w") as handle:
            SeqIO.write(records, handle, "fasta")

    def validate_file(file_type: str, file_path: str):
        file_handler_stub.validate_file_calls.append({"file_type": file_type, "file_path": file_path})
        if not Path(file_path).is_file():
            raise SystemExit(1)

    def validate_gisaid_installer(submission_dir: str, organism: str, config_dict: dict):
        file_handler_stub.validate_gisaid_installer_calls.append({"submission_dir": submission_dir, "organism": organism, "config_dict": config_dict})
        return str(Path(submission_dir) / f"{organism.lower()}CLI")

    file_handler_stub.save_csv = save_csv
    file_handler_stub.create_fasta = create_fasta
    file_handler_stub.validate_file = validate_file
    file_handler_stub.validate_gisaid_installer = validate_gisaid_installer

    settings_stub: Any = types.ModuleType("settings")
    settings_stub.GISAID_REGEX = "^gs-|^collection_date$|^authors$"
    for module_name in [
        "gisaid_handler",
        "src.gisaid_handler",
        "upload_log",
        "src.upload_log",
        "tools",
        "src.tools",
        "file_handler",
        "src.file_handler",
        "settings",
        "src.settings",
    ]:
        sys.modules.pop(module_name, None)

    monkeypatch.setitem(sys.modules, "src", src_pkg)
    alias_src_module("upload_log", upload_log_stub)
    alias_src_module("tools", tools_stub)
    alias_src_module("file_handler", file_handler_stub)
    alias_src_module("settings", settings_stub)

    sys.modules.pop("gisaid_handler", None)
    sys.modules.pop("src.gisaid_handler", None)
    spec = importlib.util.spec_from_file_location("gisaid_handler", MODULE_PATH)
    assert spec and spec.loader
    module: Any = importlib.util.module_from_spec(spec)
    sys.modules["giasid_handler"] = module
    setattr(src_pkg, "gisaid_handler", module)
    spec.loader.exec_module(module)
    return module


def assert_update_call(call: dict[str, Any], submission_dir: Path, expected_records: list[dict[str, str]]):
    assert set(call) == {"submission_dir", "update_database", "update_df"}
    assert call["submission_dir"] == str(submission_dir)
    assert call["update_database"] == "GISAID"
    assert call["update_df"] is not None
    assert call["update_df"].to_dict("records") == expected_records

def assert_no_none_strings(df: pd.DataFrame):
    assert all(column is not None for column in df.columns)
    for row in df.to_dict("records"):
        for key, value in row.items():
            assert key is not None
            assert value is not None
#*******************************************************************************
#                              process_flu_dates
#*******************************************************************************

def test_process_flu_dates__full_year_month_and_full_date(gisaid_handler_module):
    module = gisaid_handler_module

    year_only = module.process_flu_dates("2024").tolist()
    year_month = module.process_flu_dates("2024-05").tolist()
    full_date = module.process_flu_dates("2024-05-19").tolist()

    assert year_only == ["", "2024", ""]
    assert year_month == ["", "2024", "05"]
    assert full_date == ["2024-05-19", "", ""]

    assert all(value is not None for value in year_only)
    assert all(value is not None for value in year_month)
    assert all(value is not None for value in full_date)

def test_process_flu_dates__invalid_date_parts_exit(gisaid_handler_module, capsys):
    with pytest.raises(SystemExit) as exc:
        gisaid_handler_module.process_flu_dates("2024-05-19-extra")

    assert exc.value.code == 1
    assert "Error: Unable to process 'Collection_Date' column for FLU GISAID submission. The field should be in format 'YYYY-MM-DD'. Value unable to process: 2024-05-19-extra\n" == capsys.readouterr().err

#*******************************************************************************
#                         create_gisaid_files
#*******************************************************************************

def test_create_gisaid_files__cov_renames_required_fields_and_creates_files(
    tmp_path: Path,
    gisaid_handler_module,
    monkeypatch: pytest.MonkeyPatch,
):
    module = gisaid_handler_module
    copied = []

    def fake_copy(src, dst):
        copied.append((src, dst))

    monkeypatch.setattr(module.shutil, "copy", fake_copy)

    metadata = pd.DataFrame(
        {
            "gs-sample_name": ["hCoV-19/USA/1/2024"],
            "gs-covv_sex": ["Female"],
            "collection_date": ["2024-01-02"],
            "authors": ["Smith, Jane"],
            "organism": ["SARS-CoV-2"],
            "fasta_sequence_orig": ["ACGTACGT"],
        }
    )

    module.create_gisaid_files(
        organism="COV",
        database="GISAID",
        submission_name="sub1",
        submission_dir=str(tmp_path),
        config_dict={"Username": "submitter1"},
        metadata=metadata,
    )

    assert module.file_handler.saved_csvs[0]["file_name"] == "metadata.csv"
    assert module.file_handler.saved_csvs[0]["file_path"] == str(tmp_path)
    saved_df = module.file_handler.saved_csvs[0]["df"]
    assert_no_none_strings(saved_df)
    assert list(saved_df.columns[:3]) == ["submitter", "fn", "covv_virus_name"]
    assert saved_df.loc[0, "submitter"] == "submitter1"
    assert saved_df.loc[0, "fn"] == "sequence.fsa"
    assert saved_df.loc[0, "covv_gender"] == "Female"
    assert module.file_handler.created_fastas[0]["database"] == "GISAID"
    assert module.file_handler.created_fastas[0]["submission_dir"] == str(tmp_path)
    assert copied == [
        (str(tmp_path / "metadata.csv"), str(tmp_path / "orig_metadata.csv")),
        (str(tmp_path / "sequence.fsa"), str(tmp_path / "orig_sequence.fsa")),
    ]
    assert saved_df.to_dict("records") == [
        {
            "submitter": "submitter1",
            "fn": "sequence.fsa",
            "covv_virus_name": "hCoV-19/USA/1/2024",
            "covv_gender": "Female",
            "collection_date": "2024-01-02",
            "authors": "Smith, Jane",
        }
    ]

@pytest.mark.parametrize(
    ("organism", "sample_column", "sex_in", "sex_out"),
    [
        ("POX", "pox_virus_name", "pox_sex", "pox_gender"),
        ("ARBO", "arbo_virus_name", "arbo_sex", "arbo_gender"),
        ("RSV", "rsv_virus_name", "rsv_sex", "rsv_gender"),
    ],
)
def test_create_gisaid_files__non_cov_flu_organisms(
    tmp_path: Path,
    gisaid_handler_module,
    monkeypatch: pytest.MonkeyPatch,
    organism: str,
    sample_column: str,
    sex_in: str,
    sex_out: str,
):
    module = gisaid_handler_module

    def fake_copy(src, dst):
        return None

    monkeypatch.setattr(module.shutil, "copy", fake_copy)
    module.file_handler.saved_csvs.clear()

    metadata = pd.DataFrame(
        {
            "gs-sample_name": [f"{organism}/sample/1"],
            f"gs-{sex_in}": ["unknown"],
            "collection_date": ["2024"],
            "authors": ["Doe, John"],
        }
    )

    module.create_gisaid_files(
        organism=organism,
        database="GISAID",
        submission_name="sub1",
        submission_dir=str(tmp_path),
        config_dict={"Username": "submitter1"},
        metadata=metadata,
    )

    saved_df = module.file_handler.saved_csvs[0]["df"]
    assert_no_none_strings(saved_df)
    assert list(saved_df.columns[:3]) == ["submitter", "fn", sample_column]
    assert saved_df.loc[0, "submitter"] == "submitter1"
    assert saved_df.loc[0, "fn"] == "sequence.fsa"
    assert saved_df.loc[0, sample_column] == f"{organism}/sample/1"
    assert saved_df.loc[0, sex_out] == "unknown"
    assert sample_column not in ["", None]
    assert sex_out not in ["", None]

def test_create_gisaid_files__flu_pivots_segments_and_splits_dates(tmp_path: Path, gisaid_handler_module, monkeypatch: pytest.MonkeyPatch):
    module = gisaid_handler_module

    def fake_copy(src, dst):
        return None
    monkeypatch.setattr(module.shutil, "copy", fake_copy)

    metadata = pd.DataFrame(
        {
            "gs-sample_name": ["segHA", "segNA"],
            "gs-Isolate_Name": ["iso1", "iso1"],
            "gs-segment": ["HA", "NA"],
            "gs-Host_Sex": ["Female", "Female"],
            "collection_date": ["2024-05", "2024-05"],
            "authors": ["Smith, Jane", "Smith, Jane"],
        }
    )

    module.create_gisaid_files(
        organism="FLU",
        database="GISAID",
        submission_name="sub1",
        submission_dir=str(tmp_path),
        config_dict={"Username": "submitter1"},
        metadata=metadata,
    )

    saved_df = module.file_handler.saved_csvs[0]["df"]
    assert list(saved_df.columns[:3]) == ["Isolate_Id", "Segment_Ids", "Isolate_Name"]
    assert saved_df.loc[0, "Isolate_Id"] == ""
    assert saved_df.loc[0, "Segment_Ids"] == ""
    assert saved_df.loc[0, "Isolate_Name"] == "iso1"
    assert saved_df.loc[0, "Collection_Date"] == ""
    assert saved_df.loc[0, "Collection_Year"] == "2024"
    assert saved_df.loc[0, "Collection_Month"] == "05"
    assert saved_df.loc[0, "Authors"] == "Smith, Jane"
    assert saved_df.loc[0, "Host_Gender"] == "Female"
    assert saved_df.loc[0, "Seq_Id (HA)"] == "segHA"
    assert saved_df.loc[0, "Seq_Id (NA)"] == "segNA"

def test_create_gisaid_files__flu_uses_exact_segment_merge_arguments(tmp_path: Path, gisaid_handler_module, monkeypatch: pytest.MonkeyPatch):
    module = gisaid_handler_module
    monkeypatch.setattr(module.shutil, "copy", lambda src, dst: None)
    observed_merges: list[dict[str, Any]] = []
    original_merge = module.pd.DataFrame.merge

    def fake_merge(self, right, *args, **kwargs):
        observed_merges.append({"left_columns": self.columns.tolist(), "right_columns": right.columns.tolist(), "args": args, "kwargs": kwargs})
        return original_merge(self, right, *args, **kwargs)

    monkeypatch.setattr(module.pd.DataFrame, "merge", fake_merge)
    metadata = pd.DataFrame(
        {
            "gs-sample_name": ["segHA", "segNA"],
            "gs-Isolate_Name": ["iso1", "iso1"],
            "gs-segment": ["HA", "NA"],
            "gs-Host_Sex": ["Female", "Female"],
            "collection_date": ["2024-05", "2024-05"],
            "authors": ["Smith, Jane", "Smith, Jane"],
        }
    )

    module.create_gisaid_files(organism="FLU", database="GISAID", submission_name="sub1", submission_dir=str(tmp_path), config_dict={"Username": "submitter1"}, metadata=metadata)
    assert observed_merges == [
        {
            "left_columns": ["Isolate_Name", "Host_Gender", "Authors", "Collection_Date", "Collection_Year", "Collection_Month", "Isolate_Id", "Segment_Ids"],
            "right_columns": ["Isolate_Name", "Seq_Id (HA)", "Seq_Id (NA)"],
            "args": (),
            "kwargs": {
                "on": "Isolate_Name",
                "how": "inner",
                "validate": "1:1",
            },
        }
    ]
#*******************************************************************************
#                            process_gisaid_log
#*******************************************************************************
def test_process_gisaid_log__validates_exact_file_type_and_opens_read_mode(tmp_path: Path, gisaid_handler_module, monkeypatch: pytest.MonkeyPatch):
    module = gisaid_handler_module
    log_file = tmp_path / "gisaid.log"
    log_file.write_text('"msg": "hCoV-19/USA/CA-10/2024; EPI_ISL_101010"\n')

    original_open = open
    observed_open_calls: list[tuple[str, str]] = []

    def fake_open(file, mode="r", *args, **kwargs):
        observed_open_calls.append((str(file), mode))
        return original_open(file, mode, *args, **kwargs)

    monkeypatch.setattr("builtins.open", fake_open)

    returned = module.process_gisaid_log(log_file=str(log_file), submission_dir=str(tmp_path))

    assert returned.empty
    assert module.file_handler.validate_file_calls == [{"file_type": "GISAID log", "file_path": str(log_file)}]
    assert observed_open_calls == [(str(log_file), "r")]
def test_process_gisaid_log__standard_epi_id_regex_branch_uses_exact_uppercase_pattern(tmp_path: Path, gisaid_handler_module, monkeypatch: pytest.MonkeyPatch, capsys):
    module = gisaid_handler_module
    log_file = tmp_path / "gisaid.log"
    log_file.write_text("epi_id: hCoV-19/USA/CA-21/2024; EPI12345\n")

    def fake_findall(pattern, line):
        assert pattern == r'(?:[a-zA-Z0-9_-]+(?:/[a-zA-Z0-9_-]+)+|EPI_\w*)'
        return ["hCoV-19/USA/CA-21/2024", "EPI12345"]

    monkeypatch.setattr(module.re, "findall", fake_findall)

    with pytest.raises(KeyError, match="gisaid_accession_epi_isl_id"):
        module.process_gisaid_log(log_file=str(log_file), submission_dir=str(tmp_path))

    captured = capsys.readouterr()
    assert captured.out == "GISAID segments found.\n"
    assert len(module.upload_log.update_calls) == 1
    assert module.upload_log.update_calls[0]["update_database"] == "GISAID"
    assert module.upload_log.update_calls[0]["submission_dir"] == str(tmp_path)
    assert module.upload_log.update_calls[0]["update_df"].to_dict("records") == [
        {
            "gs-segment_name": "hCoV-19/USA/CA-21/2024",
            "gisaid_accession_epi_id": "EPI12345",
        }
    ]

def test_process_gisaid_log__updates_for_isolate_and_existing_segment(tmp_path: Path, gisaid_handler_module, capsys):
    module = gisaid_handler_module
    log_file = tmp_path / "gisaid.log"
    log_file.write_text(
        '\n'.join(
            [
                '"msg": "hCoV-19/USA/CA-1/2024; EPI_ISL_12345"',
                '"code": "validation_error" sample already exists; existing_virus_name: hCoV-19/USA/CA-2/2024; [\'EPI_67890\']',
            ]
        )
    )

    returned = module.process_gisaid_log(log_file=str(log_file), submission_dir=str(tmp_path))

    assert returned.empty
    assert capsys.readouterr().out == "GISAID isolates and GISAID segments found.\n"
    assert len(module.upload_log.update_calls) == 2
    assert module.upload_log.update_calls[0]["update_database"] == "GISAID"
    assert_update_call(
        module.upload_log.update_calls[0],
        tmp_path,
        [
            {
                "gs-sample_name": "hCoV-19/USA/CA-1/2024",
                "gisaid_accession_epi_isl_id": "EPI_ISL_12345",
            }
        ],
    )

    assert_update_call(
        module.upload_log.update_calls[1],
        tmp_path,
        [
            {
                "gs-segment_name": "hCoV-19/USA/CA-2/2024",
                "gisaid_accession_epi_id": "EPI_67890",
            }
        ],
    )

def test_process_gisaid_log__returns_only_failed_isolate_rows_and_filters_successes(tmp_path: Path, gisaid_handler_module, monkeypatch: pytest.MonkeyPatch):
    module = gisaid_handler_module
    log_file = tmp_path / "gisaid.log"
    log_file.write_text("placeholder\n")
    original_dataframe = module.pd.DataFrame

    def fake_validate_file(file_type: str, file_path: str):
        module.file_handler.validate_file_calls.append({"file_type": file_type, "file_path": file_path})

    monkeypatch.setattr(module.file_handler, "validate_file", fake_validate_file)

    class FakeFile:
        def __init__(self):
            self.lines = iter(["placeholder\n", ""])

        def readline(self):
            return next(self.lines)

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    monkeypatch.setattr("builtins.open", lambda *args, **kwargs: FakeFile())

    def fake_search(pattern, line):
        return False

    monkeypatch.setattr(module.re, "search", fake_search)
    call_count = {"value": 0}

    def fake_dataframe(rows=None, *args, **kwargs):
        call_count["value"] += 1
        if call_count["value"] == 1:
            return original_dataframe(
                [
                    {
                        "gs-sample_name": "success",
                        "gisaid_accession_epi_isl_id": "EPI_ISL_123",
                    },
                    {
                        "gs-sample_name": "failed",
                        "gisaid_accession_epi_isl_id": "",
                    },
                    {
                        "gs-sample_name": "missing",
                        "gisaid_accession_epi_isl_id": None,
                    },
                ]
            )
        return original_dataframe(rows, *args, **kwargs)

    monkeypatch.setattr(module.pd, "DataFrame", fake_dataframe)
    returned = module.process_gisaid_log(log_file=str(log_file), submission_dir=str(tmp_path))
    assert returned.to_dict("records") == [
        {"gs-sample_name": "failed"},
        {"gs-sample_name": "missing"},
    ]

def test_process_gisaid_log__returns_failed_isolate_rows_without_accessions(tmp_path: Path, gisaid_handler_module):
    module = gisaid_handler_module
    log_file = tmp_path / "gisaid.log"
    log_file.write_text(
        '"code": "validation_error" sample already exists; existing_virus_name: hCoV-19/USA/CA-3/2024; []\n'
    )

    # This line matches the duplicate-sample branch but has no accession. The code records
    # no isolate and currently falls through to the empty-dataframe bug covered below.
    with pytest.raises(KeyError):
        module.process_gisaid_log(log_file=str(log_file), submission_dir=str(tmp_path))

def test_process_gisaid_log__no_accessions_currently_raises_keyerror(tmp_path: Path, gisaid_handler_module, capsys):
    module = gisaid_handler_module
    log_file = tmp_path / "gisaid.log"
    log_file.write_text("no accession on this line\n")

    with pytest.raises(KeyError):
        module.process_gisaid_log(log_file=str(log_file), submission_dir=str(tmp_path))

    assert capsys.readouterr().out == (
        "Finished reading GISAID log. If workflow has failed here, "
        "it's likely no GISAID IDs were returned. Check results in GISAID upload log.\n"
        "Warning: no GISAID isolates or segments found\n"
    )

#*******************************************************************************
#                           submit_gisaid
#*******************************************************************************

def test_submit_gisaid__success_when_status_report_has_all_accessions(tmp_path: Path, gisaid_handler_module, monkeypatch: pytest.MonkeyPatch, capsys):
    module = gisaid_handler_module
    submission_dir = tmp_path / "GISAID"
    submission_dir.mkdir()
    (submission_dir / "metadata.csv").write_text("placeholder\n")
    (submission_dir / "sequence.fsa").write_text(">s1\nACGT\n")
    pd.DataFrame({"covv_virus_name": ["hCoV-19/USA/CA-1/2024"]}).to_csv(submission_dir / "orig_metadata.csv", index=False)
    (submission_dir / "orig_sequence.fsa").write_text(">hCoV-19/USA/CA-1/2024\nACGT\n")
    pd.DataFrame(
        {
            "gs-sample_name": ["hCoV-19/USA/CA-1/2024"],
            "gisaid_accession_epi_isl_id": ["EPI_ISL_12345"],
        }
    ).to_csv(tmp_path / "submission_status_report.csv", index=False)

    sleep_calls: list[int] = []

    def fake_sleep(seconds):
        sleep_calls.append(seconds)

    def fake_process_gisaid_log(**kwargs):
        return pd.DataFrame()

    monkeypatch.setattr(module.time, "sleep", fake_sleep)
    monkeypatch.setattr(module, "process_gisaid_log", fake_process_gisaid_log)

    def fake_run(command, cwd, stdout, stderr):
        assert command == [
            str(submission_dir / "covCLI"),
            "upload",
            "--username",
            "u",
            "--password",
            "p",
            "--clientid",
            "cid",
            "--metadata",
            str(submission_dir / "metadata.csv"),
            "--fasta",
            str(submission_dir / "sequence.fsa"),
            "--log",
            str(submission_dir / "gisaid_upload_log_1.txt"),
            "--debug",
        ]
        assert cwd == str(submission_dir)
        assert stdout == module.subprocess.PIPE
        assert stderr == module.subprocess.PIPE

        log_path = Path(command[command.index("--log") + 1])
        log_path.write_text("ok\n")
        return types.SimpleNamespace(returncode=0, stdout=b"", stderr=b"")

    monkeypatch.setattr(module.subprocess, "run", fake_run)

    status = module.submit_gisaid(
        organism="COV",
        submission_dir=str(submission_dir),
        submission_name="sub1",
        config_dict={"Username": "u", "Password": "p", "Client-Id": "cid"},
        submission_type="TEST",
    )

    assert status == "PROCESSED"
    assert sleep_calls == [5]
    assert capsys.readouterr().out == (
        "Uploading sample files to GISAID-COV, as a 'TEST' submission. If this is not intended, interrupt immediately.\n"
        "\n"
        "Submission attempt: 1\n"
        "Uploading successfully\n"
        f"Log file is stored at: {submission_dir}/gisaid_upload_log_attempt_1.txt\n"
    )
    assert module.tools.check_credentials_calls == [
        {"config_dict": {"Username": "u", "Password": "p", "Client-Id": "cid"}, "database": "GISAID"}
    ]

def test_submit_gisaid__subprocess_failure_exits(tmp_path: Path, gisaid_handler_module, monkeypatch: pytest.MonkeyPatch, capsys):
    module = gisaid_handler_module
    submission_dir = tmp_path / "GISAID"
    submission_dir.mkdir()
    (submission_dir / "metadata.csv").write_text("placeholder\n")
    (submission_dir / "orig_metadata.csv").write_text("covv_virus_name\ns1\n")
    (submission_dir / "sequence.fsa").write_text(">s1\nACGT\n")
    (submission_dir / "orig_sequence.fsa").write_text(">s1\nACGT\n")
    (tmp_path / "submission_status_report.csv").write_text("gs-sample_name,gisaid_accession_epi_isl_id\ns1,\n")
    sleep_calls: list[int] = []

    def fake_sleep(seconds):
        sleep_calls.append(seconds)

    def fake_run(*args, **kwargs):
        return types.SimpleNamespace(returncode=1, stdout=b"bad out", stderr=b"bad err")
    monkeypatch.setattr(module.time, "sleep", fake_sleep)
    monkeypatch.setattr(module.subprocess, "run", fake_run)

    with pytest.raises(SystemExit) as exc:
        module.submit_gisaid(
            organism="COV",
            submission_dir=str(submission_dir),
            submission_name="sub1",
            config_dict={"Username": "u", "Password": "p", "Client-Id": "cid"},
            submission_type="TEST",
        )
    assert sleep_calls == [5]
    captured = capsys.readouterr()
    assert exc.value.code == 1
    assert captured.out == (
        "Uploading sample files to GISAID-COV, as a 'TEST' submission. If this is not intended, interrupt immediately.\n"
        "\n"
        "Submission attempt: 1\n"
        "b'bad out'\n"
        "b'bad err'\n"
    )
    assert captured.err == "Error: upload command error\n"

#*******************************************************************************
#                           update_gisaid_files
#*******************************************************************************

@pytest.mark.parametrize(
    ("organism", "metadata_column", "fasta_name", "status_rows"),
    [
        ("COV", "covv_virus_name", "sample1", {"gs-sample_name": ["sample1"]}),
        ("POX", "pox_virus_name", "sample1", {"gs-sample_name": ["sample1"]}),
        (
            "FLU",
            "Isolate_Name",
            "segHA",
            {"gs-sample_name": ["iso1"], "gs-segment_name": ["segHA"]},
        ),
    ],
)
def test_update_gisaid_files__filters_to_processed_genbank_rows_and_rewrites_files(
    tmp_path: Path,
    gisaid_handler_module,
    organism: str,
    metadata_column: str,
    fasta_name: str,
    status_rows: dict,
):
    module = gisaid_handler_module
    submission_dir = tmp_path / "GISAID"
    submission_dir.mkdir()

    original_metadata_value = "iso1" if organism == "FLU" else "sample1"
    pd.DataFrame({metadata_column: [original_metadata_value, "dropme"]}).to_csv(
        submission_dir / "metadata.csv", index=False
    )
    SeqIO.write(
        [
            SeqRecord(Seq("AAAA"), id=fasta_name, description=""),
            SeqRecord(Seq("CCCC"), id="dropme", description=""),
        ],
        submission_dir / "orig_sequence.fsa",
        "fasta",
    )
    pd.DataFrame(
        {
            "genbank-status": ["processed-ok", "failed"],
            **{key: values + ["dropme"] for key, values in status_rows.items()},
        }
    ).to_csv(tmp_path / "submission_status_report.csv", index=False)

    module.update_gisaid_files(
        organism=organism,
        submission_dir=str(submission_dir),
        submission_status_file=str(tmp_path / "submission_status_report.csv"),
    )

    assert (submission_dir / "metadata.csv").exists()
    assert (submission_dir / "orig_metadata.csv").exists()
    assert (submission_dir / "sequence.fsa").exists()
    assert (submission_dir / "orig_sequence.fsa").exists()
    rewritten_metadata = pd.read_csv(submission_dir / "orig_metadata.csv", dtype=str)
    assert rewritten_metadata[metadata_column].tolist() == [original_metadata_value]
    rewritten_records = list(SeqIO.parse(submission_dir / "sequence.fsa", "fasta"))
    assert [record.id for record in rewritten_records] == [fasta_name]

def test_update_gisaid_files__genbank_status_filter_uses_na_false(tmp_path: Path, gisaid_handler_module):
    module = gisaid_handler_module
    submission_dir = tmp_path / "GISAID"
    submission_dir.mkdir()

    pd.DataFrame({"covv_virus_name": ["sample1", "sample2"]}).to_csv(
        submission_dir / "metadata.csv",
        index=False,
    )
    SeqIO.write(
        [
            SeqRecord(Seq("AAAA"), id="sample1", description=""),
            SeqRecord(Seq("CCCC"), id="sample2", description=""),
        ],
        submission_dir / "orig_sequence.fsa",
        "fasta",
    )

    pd.DataFrame(
        {
            "genbank-status": ["processed-ok", None],
            "gs-sample_name": ["sample1", "sample2"],
        }
    ).to_csv(tmp_path / "submission_status_report.csv", index=False)

    module.update_gisaid_files(
        organism="COV",
        submission_dir=str(submission_dir),
        submission_status_file=str(tmp_path / "submission_status_report.csv"),
    )

    rewritten_metadata = pd.read_csv(submission_dir / "orig_metadata.csv", dtype=str)
    assert rewritten_metadata["covv_virus_name"].tolist() == ["sample1"]

    rewritten_fasta_ids = [
        record.id for record in SeqIO.parse(submission_dir / "sequence.fsa", "fasta")
    ]
    assert rewritten_fasta_ids == ["sample1"]

def test_update_gisaid_files__uses_exact_read_merge_drop_to_csv_and_open_args(tmp_path: Path, gisaid_handler_module, monkeypatch: pytest.MonkeyPatch):
    module = gisaid_handler_module
    submission_dir = tmp_path / "GISAID"
    submission_dir.mkdir()

    metadata_path = submission_dir / "metadata.csv"
    orig_metadata_path = submission_dir / "orig_metadata.csv"
    fasta_path = submission_dir / "sequence.fsa"
    orig_fasta_path = submission_dir / "orig_sequence.fsa"
    status_path = tmp_path / "submission_status_report.csv"

    pd.DataFrame({"covv_virus_name": ["sample1", "sample2"]}).to_csv(metadata_path, index=False)
    SeqIO.write(
        [
            SeqRecord(Seq("AAAA"), id="sample1", description=""),
            SeqRecord(Seq("CCCC"), id="sample2", description=""),
        ],
        orig_fasta_path,
        "fasta",
    )
    pd.DataFrame(
        {
            "genbank-status": ["processed-ok", "failed"],
            "gs-sample_name": ["sample1", "sample2"],
        }
    ).to_csv(status_path, index=False)

    original_read_csv = module.pd.read_csv
    read_csv_calls: list[dict[str, Any]] = []

    def fake_read_csv(*args, **kwargs):
        read_csv_calls.append({"args": args, "kwargs": kwargs})
        return original_read_csv(*args, **kwargs)

    monkeypatch.setattr(module.pd, "read_csv", fake_read_csv)
    original_merge = module.pd.DataFrame.merge
    merge_calls: list[dict[str, Any]] = []

    def fake_merge(self, right, *args, **kwargs):
        merge_calls.append(
            {
                "left_columns": self.columns.tolist(),
                "right_columns": right.columns.tolist(),
                "args": args,
                "kwargs": kwargs,
            }
        )
        return original_merge(self, right, *args, **kwargs)

    monkeypatch.setattr(module.pd.DataFrame, "merge", fake_merge)
    original_drop = module.pd.DataFrame.drop
    drop_calls: list[dict[str, Any]] = []

    def fake_drop(self, *args, **kwargs):
        drop_calls.append({"args": args, "kwargs": kwargs})
        return original_drop(self, *args, **kwargs)

    monkeypatch.setattr(module.pd.DataFrame, "drop", fake_drop)
    original_to_csv = module.pd.DataFrame.to_csv
    to_csv_calls: list[dict[str, Any]] = []

    def fake_to_csv(self, *args, **kwargs):
        to_csv_calls.append({"args": args, "kwargs": kwargs})
        return original_to_csv(self, *args, **kwargs)

    monkeypatch.setattr(module.pd.DataFrame, "to_csv", fake_to_csv)
    original_open = open
    open_calls: list[tuple[str, str]] = []

    def fake_open(file, mode="r", *args, **kwargs):
        if str(file).endswith(("orig_sequence.fsa", "sequence.fsa")):
            open_calls.append((str(file), mode))
        return original_open(file, mode, *args, **kwargs)

    monkeypatch.setattr("builtins.open", fake_open)
    module.update_gisaid_files(organism="COV", submission_dir=str(submission_dir), submission_status_file=str(status_path))
    assert read_csv_calls == [
        {
            "args": (str(status_path),),
            "kwargs": {
                "header": 0,
                "dtype": str,
                "engine": "python",
                "encoding": "utf-8",
                "index_col": False,
            },
        },
        {
            "args": (str(metadata_path),),
            "kwargs": {
                "header": 0,
                "dtype": str,
                "engine": "python",
                "encoding": "utf-8",
                "index_col": False,
            },
        },
    ]

    assert merge_calls[0]["kwargs"] == {
        "how": "inner",
        "left_on": "covv_virus_name",
        "right_on": "gs-sample_name",
    }

    assert drop_calls[0]["kwargs"] == {
        "columns": ["gs-sample_name", "gs-segment_name"],
        "errors": "ignore",
    }

    assert to_csv_calls[-1] == {
        "args": (str(orig_metadata_path),),
        "kwargs": {
            "header": True,
            "index": False,
        },
    }

    assert open_calls == [
        (str(orig_fasta_path), "r"),
        (str(fasta_path), "w+"),
    ]

def test_update_gisaid_files__no_processed_rows_writes_empty_files(
    tmp_path: Path,
    gisaid_handler_module,
):
    module = gisaid_handler_module
    submission_dir = tmp_path / "GISAID"
    submission_dir.mkdir()
    pd.DataFrame({"covv_virus_name": ["sample1"]}).to_csv(submission_dir / "metadata.csv", index=False)
    SeqIO.write([SeqRecord(Seq("AAAA"), id="sample1", description="")], submission_dir / "orig_sequence.fsa", "fasta")
    pd.DataFrame(
        {
            "genbank-status": ["failed"],
            "gs-sample_name": ["sample1"],
        }
    ).to_csv(tmp_path / "submission_status_report.csv", index=False)

    module.update_gisaid_files(
        organism="COV",
        submission_dir=str(submission_dir),
        submission_status_file=str(tmp_path / "submission_status_report.csv"),
    )

    assert pd.read_csv(submission_dir / "orig_metadata.csv", dtype=str).empty
    assert list(SeqIO.parse(submission_dir / "sequence.fsa", "fasta")) == []

def test_process_flu_dates__trims_whitespace(gisaid_handler_module):
    assert gisaid_handler_module.process_flu_dates(" 2024-07-03 ").tolist() == ["2024-07-03", "", ""]


def test_create_gisaid_files__unsupported_other_organism_currently_raises_unboundlocalerror(
    tmp_path: Path,
    gisaid_handler_module,
    monkeypatch: pytest.MonkeyPatch,
):
    module = gisaid_handler_module

    def fake_copy(src, dst):
        return None

    monkeypatch.setattr(module.shutil, "copy", fake_copy)

    metadata = pd.DataFrame(
        {
            "gs-sample_name": ["OTHER/sample/1"],
            "collection_date": ["2024"],
            "authors": ["Doe, John"],
        }
    )

    # Current behavior: create_gisaid_files has branches for COV/POX/ARBO/RSV and FLU,
    # but not OTHER, leaving first_cols undefined.
    with pytest.raises(UnboundLocalError):
        module.create_gisaid_files(
            organism="OTHER",
            database="GISAID",
            submission_name="sub1",
            submission_dir=str(tmp_path),
            config_dict={"Username": "submitter1"},
            metadata=metadata,
        )


def test_create_gisaid_files__flu_full_date_populates_collection_date_column(
    tmp_path: Path,
    gisaid_handler_module,
    monkeypatch: pytest.MonkeyPatch,
):
    module = gisaid_handler_module

    def fake_copy(src, dst):
        return None
    monkeypatch.setattr(module.shutil, "copy", fake_copy)
    module.file_handler.saved_csvs.clear()

    metadata = pd.DataFrame(
        {
            "gs-sample_name": ["segPB2"],
            "gs-Isolate_Name": ["iso-full-date"],
            "gs-segment": ["PB2"],
            "gs-Host_Sex": ["Male"],
            "collection_date": ["2024-07-03"],
            "authors": ["Smith, Jane"],
        }
    )

    module.create_gisaid_files(
        organism="FLU",
        database="GISAID",
        submission_name="sub1",
        submission_dir=str(tmp_path),
        config_dict={"Username": "submitter1"},
        metadata=metadata,
    )

    saved_df = module.file_handler.saved_csvs[0]["df"]
    assert saved_df.loc[0, "Collection_Date"] == "2024-07-03"
    assert saved_df.loc[0, "Collection_Year"] == ""
    assert saved_df.loc[0, "Collection_Month"] == ""
    assert saved_df.loc[0, "Seq_Id (PB2)"] == "segPB2"


def test_process_gisaid_log__missing_log_file_exits(tmp_path: Path, gisaid_handler_module):
    with pytest.raises(SystemExit):
        gisaid_handler_module.process_gisaid_log(
            log_file=str(tmp_path / "missing.log"),
            submission_dir=str(tmp_path),
        )
def test_submit_gisaid__post_loop_empty_metadata_returns_exact_processed(tmp_path: Path, gisaid_handler_module, monkeypatch: pytest.MonkeyPatch):
    module = gisaid_handler_module
    submission_dir = tmp_path / "GISAID"
    submission_dir.mkdir()

    pd.DataFrame({"covv_virus_name": ["sample1"]}).to_csv(
        submission_dir / "orig_metadata.csv",
        index=False,
    )
    pd.DataFrame({"covv_virus_name": ["sample1"]}).to_csv(
        submission_dir / "metadata.csv",
        index=False,
    )
    SeqIO.write(
        [SeqRecord(Seq("AAAA"), id="sample1", description="")],
        submission_dir / "orig_sequence.fsa",
        "fasta",
    )
    SeqIO.write(
        [SeqRecord(Seq("AAAA"), id="sample1", description="")],
        submission_dir / "sequence.fsa",
        "fasta",
    )

    pd.DataFrame(
        {
            "gs-sample_name": ["sample1"],
            "gisaid_accession_epi_isl_id": [""],
        }
    ).to_csv(tmp_path / "submission_status_report.csv", index=False)

    monkeypatch.setattr(module.time, "sleep", lambda seconds: None)
    monkeypatch.setattr(module, "process_gisaid_log", lambda **kwargs: pd.DataFrame())

    def fake_run(command, cwd, stdout, stderr):
        Path(command[command.index("--log") + 1]).write_text("ok\n")

        # On attempt 4, mutate the status file so the loop's final metadata_df becomes empty.
        if command[command.index("--log") + 1].endswith("gisaid_upload_log_4.txt"):
            pd.DataFrame(
                {
                    "gs-sample_name": ["sample1"],
                    "gisaid_accession_epi_isl_id": ["EPI_ISL_FINAL"],
                }
            ).to_csv(tmp_path / "submission_status_report.csv", index=False)

        return types.SimpleNamespace(returncode=0, stdout=b"", stderr=b"")

    monkeypatch.setattr(module.subprocess, "run", fake_run)

    status = module.submit_gisaid(
        organism="COV",
        submission_dir=str(submission_dir),
        submission_name="sub1",
        config_dict={"Username": "u", "Password": "p", "Client-Id": "cid"},
        submission_type="TEST",
    )

    assert status == "PROCESSED"

def test_process_gisaid_log__existing_duplicate_isolate_accession_updates_status(tmp_path: Path, gisaid_handler_module, capsys):
    module = gisaid_handler_module
    log_file = tmp_path / "gisaid.log"
    log_file.write_text(
        '"code": "validation_error" sample already exists; '
        "existing_virus_name: hCoV-19/USA/CA-3/2024; ['EPI_ISL_99999']\n"
    )

    returned = module.process_gisaid_log(log_file=str(log_file), submission_dir=str(tmp_path))

    assert returned.empty
    assert capsys.readouterr().out == "GISAID isolates found.\n"
    assert len(module.upload_log.update_calls) == 1
    assert_update_call(
    module.upload_log.update_calls[0],
    tmp_path,
    [
        {
            "gs-sample_name": "hCoV-19/USA/CA-3/2024",
            "gisaid_accession_epi_isl_id": "EPI_ISL_99999",
        }
    ],
)

def test_process_gisaid_log__standard_epi_id_message_currently_not_returned_as_segment(
    tmp_path: Path,
    gisaid_handler_module,
):
    module = gisaid_handler_module
    log_file = tmp_path / "gisaid.log"
    log_file.write_text('epi_id: hCoV-19/USA/CA-4/2024; EPI_ID_123456\n')

    # Current behavior: the line matches the broad regex, but the extracted accession
    # EPI_ID_123456 does not satisfy re.match(r"EPI\d+", ...), so no segment row is
    # recorded and the empty-dataframe bug is reached.
    with pytest.raises(KeyError):
        module.process_gisaid_log(log_file=str(log_file), submission_dir=str(tmp_path))

def test_submit_gisaid__partial_failure_rewrites_remaining_samples_then_retries_successfully(
    tmp_path: Path,
    gisaid_handler_module,
    monkeypatch: pytest.MonkeyPatch,
):
    module = gisaid_handler_module
    submission_dir = tmp_path / "GISAID"
    submission_dir.mkdir()
    pd.DataFrame({"covv_virus_name": ["sample1", "sample2"]}).to_csv(
        submission_dir / "orig_metadata.csv", index=False
    )
    pd.DataFrame({"covv_virus_name": ["sample1", "sample2"]}).to_csv(
        submission_dir / "metadata.csv", index=False
    )
    SeqIO.write(
        [SeqRecord(Seq("AAAA"), id="sample1", description=""), SeqRecord(Seq("CCCC"), id="sample2", description="")],
        submission_dir / "orig_sequence.fsa",
        "fasta",
    )
    (submission_dir / "sequence.fsa").write_text(">sample1\nAAAA\n>sample2\nCCCC\n")
    status_file = tmp_path / "submission_status_report.csv"
    pd.DataFrame(
        {
            "gs-sample_name": ["sample1", "sample2"],
            "gisaid_accession_epi_isl_id": ["EPI_ISL_111", ""],
        }
    ).to_csv(status_file, index=False)
    sleep_calls: list[int] = []

    def fake_sleep(seconds):
        sleep_calls.append(seconds)

    monkeypatch.setattr(module.time, "sleep", fake_sleep)
    process_calls = {"count": 0}

    def fake_process_gisaid_log(**kwargs):
        process_calls["count"] += 1
        if process_calls["count"] == 2:
            pd.DataFrame(
                {
                    "gs-sample_name": ["sample1", "sample2"],
                    "gisaid_accession_epi_isl_id": ["EPI_ISL_111", "EPI_ISL_222"],
                }
            ).to_csv(status_file, index=False)
        return pd.DataFrame()

    def fake_run(command, cwd, stdout, stderr):
        log_path = Path(command[command.index("--log") + 1])
        log_path.write_text("ok\n")
        return types.SimpleNamespace(returncode=0, stdout=b"", stderr=b"")

    monkeypatch.setattr(module, "process_gisaid_log", fake_process_gisaid_log)
    monkeypatch.setattr(module.subprocess, "run", fake_run)

    status = module.submit_gisaid(
        organism="COV",
        submission_dir=str(submission_dir),
        submission_name="sub1",
        config_dict={"Username": "u", "Password": "p", "Client-Id": "cid"},
        submission_type="TEST",
    )

    assert status == "PROCESSED"
    assert process_calls["count"] == 2
    # The first retry should rewrite the active files to only the unsubmitted sample.
    active_metadata = pd.read_csv(submission_dir / "metadata.csv", dtype=str)
    assert active_metadata["covv_virus_name"].tolist() == ["sample2"]
    active_fasta_ids = [record.id for record in SeqIO.parse(submission_dir / "sequence.fsa", "fasta")]
    assert active_fasta_ids == ["sample2"]


def test_submit_gisaid__returns_error_after_all_retry_attempts_have_remaining_samples(tmp_path: Path, gisaid_handler_module, monkeypatch: pytest.MonkeyPatch, capsys):
    module = gisaid_handler_module
    submission_dir = tmp_path / "GISAID"
    submission_dir.mkdir()
    pd.DataFrame({"covv_virus_name": ["sample1"]}).to_csv(submission_dir / "orig_metadata.csv", index=False)
    pd.DataFrame({"covv_virus_name": ["sample1"]}).to_csv(submission_dir / "metadata.csv", index=False)
    SeqIO.write([SeqRecord(Seq("AAAA"), id="sample1", description="")], submission_dir / "orig_sequence.fsa", "fasta")
    (submission_dir / "sequence.fsa").write_text(">sample1\nAAAA\n")
    pd.DataFrame(
        {"gs-sample_name": ["sample1"], "gisaid_accession_epi_isl_id": [""]}
    ).to_csv(tmp_path / "submission_status_report.csv", index=False)
    sleep_calls: list[int] = []

    def fake_sleep(seconds):
        sleep_calls.append(seconds)

    def fake_process_gisaid_log(**kwargs):
        return pd.DataFrame()

    monkeypatch.setattr(module.time, "sleep", fake_sleep)
    monkeypatch.setattr(module, "process_gisaid_log", fake_process_gisaid_log)
    run_calls = []

    def fake_run(command, cwd, stdout, stderr):
        run_calls.append(command)
        Path(command[command.index("--log") + 1]).write_text("ok\n")
        return types.SimpleNamespace(returncode=0, stdout=b"", stderr=b"")

    monkeypatch.setattr(module.subprocess, "run", fake_run)

    status = module.submit_gisaid(
        organism="COV",
        submission_dir=str(submission_dir),
        submission_name="sub1",
        config_dict={"Username": "u", "Password": "p", "Client-Id": "cid"},
        submission_type="TEST",
    )

    captured = capsys.readouterr()

    assert status == "ERROR"
    assert sleep_calls == [5]
    assert len(run_calls) == 4
    assert captured.out == (
        "Uploading sample files to GISAID-COV, as a 'TEST' submission. "
        "If this is not intended, interrupt immediately.\n"
        "\n"
        "Submission attempt: 1\n"
        "\n"
        "Submission attempt: 2\n"
        "\n"
        "Submission attempt: 3\n"
        "\n"
        "Submission attempt: 4\n"
    )
    assert captured.err == (
        "Error: 1 sample(s) failed to upload to GISAID\n"
        f"Please check log file at: {submission_dir}/gisaid_upload_log_attempt_{{1,2,3}}.txt\n"
    )


def test_submit_gisaid__other_organism_uses_lowercase_virus_name_for_remaining_metadata(
    tmp_path: Path,
    gisaid_handler_module,
    monkeypatch: pytest.MonkeyPatch,
):
    module = gisaid_handler_module
    submission_dir = tmp_path / "GISAID"
    submission_dir.mkdir()
    pd.DataFrame({"rsv_virus_name": ["rsv1", "rsv2"]}).to_csv(submission_dir / "orig_metadata.csv", index=False)
    pd.DataFrame({"rsv_virus_name": ["rsv1", "rsv2"]}).to_csv(submission_dir / "metadata.csv", index=False)
    SeqIO.write(
        [SeqRecord(Seq("AAAA"), id="rsv1", description=""), SeqRecord(Seq("CCCC"), id="rsv2", description="")],
        submission_dir / "orig_sequence.fsa",
        "fasta",
    )
    (submission_dir / "sequence.fsa").write_text(">rsv1\nAAAA\n>rsv2\nCCCC\n")
    status_file = tmp_path / "submission_status_report.csv"
    pd.DataFrame(
        {"gs-sample_name": ["rsv1", "rsv2"], "gisaid_accession_epi_isl_id": ["EPI_ISL_1", ""]}
    ).to_csv(status_file, index=False)
    sleep_calls: list[int] = []

    def fake_sleep(seconds):
        sleep_calls.append(seconds)

    def fake_process_gisaid_log(**kwargs):
        return pd.DataFrame()

    monkeypatch.setattr(module.time, "sleep", fake_sleep)
    monkeypatch.setattr(module, "process_gisaid_log", fake_process_gisaid_log)

    def fake_run(command, cwd, stdout, stderr):
        Path(command[command.index("--log") + 1]).write_text("ok\n")
        return types.SimpleNamespace(returncode=0, stdout=b"", stderr=b"")

    monkeypatch.setattr(module.subprocess, "run", fake_run)

    status = module.submit_gisaid(
        organism="RSV",
        submission_dir=str(submission_dir),
        submission_name="sub1",
        config_dict={"Username": "u", "Password": "p", "Client-Id": "cid"},
        submission_type="TEST",
    )

    assert status == "ERROR"
    active_metadata = pd.read_csv(submission_dir / "metadata.csv", dtype=str)
    assert active_metadata["rsv_virus_name"].tolist() == ["rsv2"]


def test_update_gisaid_files__arbo_and_rsv_use_lowercase_virus_name_columns(
    tmp_path: Path,
    gisaid_handler_module,
):
    module = gisaid_handler_module

    for organism in ["ARBO", "RSV"]:
        submission_dir = tmp_path / organism / "GISAID"
        submission_dir.mkdir(parents=True)
        metadata_column = organism.lower() + "_virus_name"
        pd.DataFrame({metadata_column: [f"{organism.lower()}1", "dropme"]}).to_csv(
            submission_dir / "metadata.csv", index=False
        )
        SeqIO.write(
            [
                SeqRecord(Seq("AAAA"), id=f"{organism.lower()}1", description=""),
                SeqRecord(Seq("CCCC"), id="dropme", description=""),
            ],
            submission_dir / "orig_sequence.fsa",
            "fasta",
        )
        status_file = submission_dir.parent / "submission_status_report.csv"
        pd.DataFrame(
            {
                "genbank-status": ["processed-ok", "failed"],
                "gs-sample_name": [f"{organism.lower()}1", "dropme"],
            }
        ).to_csv(status_file, index=False)

        module.update_gisaid_files(
            organism=organism,
            submission_dir=str(submission_dir),
            submission_status_file=str(status_file),
        )

        rewritten_metadata = pd.read_csv(submission_dir / "orig_metadata.csv", dtype=str)
        assert rewritten_metadata[metadata_column].tolist() == [f"{organism.lower()}1"]
        assert [record.id for record in SeqIO.parse(submission_dir / "sequence.fsa", "fasta")] == [
            f"{organism.lower()}1"
        ]


def test_submit_gisaid__rewrite_block_updates_active_metadata_and_fasta_before_next_attempt(
    tmp_path: Path,
    gisaid_handler_module,
    monkeypatch: pytest.MonkeyPatch,
):
    """Regression coverage for the retry block that rewrites metadata.csv and sequence.fsa."""
    module = gisaid_handler_module
    submission_dir = tmp_path / "GISAID"
    submission_dir.mkdir()

    # sample1 has an accession after the first attempt; sample2 should remain active.
    pd.DataFrame({"covv_virus_name": ["sample1", "sample2"]}).to_csv(
        submission_dir / "orig_metadata.csv", index=False
    )
    pd.DataFrame({"covv_virus_name": ["sample1", "sample2"]}).to_csv(
        submission_dir / "metadata.csv", index=False
    )
    SeqIO.write(
        [
            SeqRecord(Seq("AAAA"), id="sample1", description=""),
            SeqRecord(Seq("CCCC"), id="sample2", description=""),
        ],
        submission_dir / "orig_sequence.fsa",
        "fasta",
    )
    SeqIO.write(
        [
            SeqRecord(Seq("AAAA"), id="sample1", description=""),
            SeqRecord(Seq("CCCC"), id="sample2", description=""),
        ],
        submission_dir / "sequence.fsa",
        "fasta",
    )
    pd.DataFrame(
        {
            "gs-sample_name": ["sample1", "sample2"],
            "gisaid_accession_epi_isl_id": ["EPI_ISL_111", ""],
        }
    ).to_csv(tmp_path / "submission_status_report.csv", index=False)
    sleep_calls: list[int] = []

    def fake_sleep(seconds):
        sleep_calls.append(seconds)

    def fake_process_gisaid_log(**kwargs):
        return pd.DataFrame()

    monkeypatch.setattr(module.time, "sleep", fake_sleep)
    monkeypatch.setattr(module, "process_gisaid_log", fake_process_gisaid_log)

    run_calls = {"count": 0}

    def fake_run(command, cwd, stdout, stderr):
        run_calls["count"] += 1
        Path(command[command.index("--log") + 1]).write_text("ok\n")
        if run_calls["count"] == 2:
            active_metadata = pd.read_csv(submission_dir / "metadata.csv", dtype=str)
            active_ids = [record.id for record in SeqIO.parse(submission_dir / "sequence.fsa", "fasta")]
            assert active_metadata["covv_virus_name"].tolist() == ["sample2"]
            assert active_ids == ["sample2"]
            raise RuntimeError("stop after proving retry inputs were rewritten")
        return types.SimpleNamespace(returncode=0, stdout=b"", stderr=b"")

    monkeypatch.setattr(module.subprocess, "run", fake_run)

    with pytest.raises(RuntimeError, match="retry inputs were rewritten"):
        module.submit_gisaid(
            organism="COV",
            submission_dir=str(submission_dir),
            submission_name="sub1",
            config_dict={"Username": "u", "Password": "p", "Client-Id": "cid"},
            submission_type="TEST",
        )

def test_submit_gisaid__uses_exact_read_csv_settings_for_status_and_orig_metadata(tmp_path: Path, gisaid_handler_module, monkeypatch: pytest.MonkeyPatch):
    module = gisaid_handler_module
    submission_dir = tmp_path / "GISAID"
    submission_dir.mkdir()

    metadata_path = submission_dir / "metadata.csv"
    orig_metadata_path = submission_dir / "orig_metadata.csv"
    fasta_path = submission_dir / "sequence.fsa"
    orig_fasta_path = submission_dir / "orig_sequence.fsa"
    status_path = tmp_path / "submission_status_report.csv"

    metadata_path.write_text("covv_virus_name\nsample1\n")
    orig_metadata_path.write_text("covv_virus_name\nsample1\n")
    fasta_path.write_text(">sample1\nAAAA\n")
    orig_fasta_path.write_text(">sample1\nAAAA\n")
    status_path.write_text("gs-sample_name,gisaid_accession_epi_isl_id\nsample1,EPI_ISL_1\n")

    sleep_calls: list[int] = []
    monkeypatch.setattr(module.time, "sleep", lambda seconds: sleep_calls.append(seconds))
    monkeypatch.setattr(module, "process_gisaid_log", lambda **kwargs: pd.DataFrame())

    def fake_run(command, cwd, stdout, stderr):
        Path(command[command.index("--log") + 1]).write_text("ok\n")
        return types.SimpleNamespace(returncode=0, stdout=b"", stderr=b"")

    monkeypatch.setattr(module.subprocess, "run", fake_run)

    original_read_csv = module.pd.read_csv
    observed_read_csv: list[dict[str, Any]] = []

    def fake_read_csv(*args, **kwargs):
        observed_read_csv.append({"args": args, "kwargs": kwargs})
        return original_read_csv(*args, **kwargs)

    monkeypatch.setattr(module.pd, "read_csv", fake_read_csv)

    status = module.submit_gisaid(
        organism="COV",
        submission_dir=str(submission_dir),
        submission_name="sub1",
        config_dict={"Username": "u", "Password": "p", "Client-Id": "cid"},
        submission_type="TEST",
    )

    assert status == "PROCESSED"
    assert observed_read_csv == [
        {
            "args": (str(status_path),),
            "kwargs": {
                "header": 0,
                "dtype": str,
                "engine": "python",
                "encoding": "utf-8",
                "index_col": False,
            },
        },
        {
            "args": (str(orig_metadata_path),),
            "kwargs": {
                "header": 0,
                "dtype": str,
                "engine": "python",
                "encoding": "utf-8",
                "index_col": False,
            },
        },
    ]

def test_submit_gisaid__retry_rewrite_uses_exact_merge_drop_to_csv_and_open_args(tmp_path: Path, gisaid_handler_module, monkeypatch: pytest.MonkeyPatch):
    module = gisaid_handler_module
    submission_dir = tmp_path / "GISAID"
    submission_dir.mkdir()
    pd.DataFrame({"covv_virus_name": ["sample1", "sample2"]}).to_csv(submission_dir / "orig_metadata.csv", index=False)
    pd.DataFrame({"covv_virus_name": ["sample1", "sample2"]}).to_csv(submission_dir / "metadata.csv", index=False)
    SeqIO.write(
        [
            SeqRecord(Seq("AAAA"), id="sample1", description=""),
            SeqRecord(Seq("CCCC"), id="sample2", description=""),
        ],
        submission_dir / "orig_sequence.fsa",
        "fasta",
    )
    SeqIO.write(
        [
            SeqRecord(Seq("AAAA"), id="sample1", description=""),
            SeqRecord(Seq("CCCC"), id="sample2", description=""),
        ],
        submission_dir / "sequence.fsa",
        "fasta",
    )
    status_path = tmp_path / "submission_status_report.csv"
    pd.DataFrame(
        {
            "gs-sample_name": ["sample1", "sample2"],
            "gisaid_accession_epi_isl_id": ["EPI_ISL_111", ""],
        }
    ).to_csv(status_path, index=False)

    sleep_calls: list[int] = []
    monkeypatch.setattr(module.time, "sleep", lambda seconds: sleep_calls.append(seconds))
    monkeypatch.setattr(module, "process_gisaid_log", lambda **kwargs: pd.DataFrame())
    run_calls = {"count": 0}

    def fake_run(command, cwd, stdout, stderr):
        run_calls["count"] += 1
        Path(command[command.index("--log") + 1]).write_text("ok\n")
        if run_calls["count"] == 2:
            raise RuntimeError("stop after first rewrite")
        return types.SimpleNamespace(returncode=0, stdout=b"", stderr=b"")

    monkeypatch.setattr(module.subprocess, "run", fake_run)
    original_merge = module.pd.DataFrame.merge
    merge_calls: list[dict[str, Any]] = []

    def fake_merge(self, right, *args, **kwargs):
        merge_calls.append(
            {
                "left_columns": self.columns.tolist(),
                "right_columns": right.columns.tolist(),
                "args": args,
                "kwargs": kwargs,
            }
        )
        return original_merge(self, right, *args, **kwargs)

    monkeypatch.setattr(module.pd.DataFrame, "merge", fake_merge)
    original_drop = module.pd.DataFrame.drop
    drop_calls: list[dict[str, Any]] = []

    def fake_drop(self, *args, **kwargs):
        drop_calls.append({"args": args, "kwargs": kwargs})
        return original_drop(self, *args, **kwargs)

    monkeypatch.setattr(module.pd.DataFrame, "drop", fake_drop)
    original_to_csv = module.pd.DataFrame.to_csv
    to_csv_calls: list[dict[str, Any]] = []

    def fake_to_csv(self, *args, **kwargs):
        to_csv_calls.append({"args": args, "kwargs": kwargs})
        return original_to_csv(self, *args, **kwargs)

    monkeypatch.setattr(module.pd.DataFrame, "to_csv", fake_to_csv)
    original_open = open
    open_calls: list[tuple[str, str]] = []

    def fake_open(file, mode="r", *args, **kwargs):
        if str(file).endswith(("orig_sequence.fsa", "sequence.fsa")):
            open_calls.append((str(file), mode))
        return original_open(file, mode, *args, **kwargs)

    monkeypatch.setattr("builtins.open", fake_open)
    with pytest.raises(RuntimeError, match="first rewrite"):
        module.submit_gisaid(
            organism="COV",
            submission_dir=str(submission_dir),
            submission_name="sub1",
            config_dict={"Username": "u", "Password": "p", "Client-Id": "cid"},
            submission_type="TEST",
        )

    assert merge_calls[0]["kwargs"] == {
        "how": "inner",
        "left_on": "covv_virus_name",
        "right_on": "gs-sample_name",
    }

    assert drop_calls[0]["kwargs"] == {
        "columns": ["gs-sample_name", "gs-segment_name"],
        "errors": "ignore",
    }

    assert to_csv_calls[-1] == {
        "args": (os.path.join(str(submission_dir), "metadata.csv"),),
        "kwargs": {
            "header": True,
            "index": False,
        },
    }

    assert open_calls == [
        (str(submission_dir / "orig_sequence.fsa"), "r"),
        (str(submission_dir / "sequence.fsa"), "w+"),
    ]

def test_submit_gisaid__after_retry_loop_returns_error_and_leaves_remaining_active_files(
    tmp_path: Path,
    gisaid_handler_module,
    monkeypatch: pytest.MonkeyPatch,
    capsys,
):
    """Covers the post-loop `if not metadata_df.empty: return "ERROR"` branch."""
    module = gisaid_handler_module
    submission_dir = tmp_path / "GISAID"
    submission_dir.mkdir()

    pd.DataFrame({"covv_virus_name": ["sample1", "sample2"]}).to_csv(
        submission_dir / "orig_metadata.csv", index=False
    )
    pd.DataFrame({"covv_virus_name": ["sample1", "sample2"]}).to_csv(
        submission_dir / "metadata.csv", index=False
    )
    SeqIO.write(
        [
            SeqRecord(Seq("AAAA"), id="sample1", description=""),
            SeqRecord(Seq("CCCC"), id="sample2", description=""),
        ],
        submission_dir / "orig_sequence.fsa",
        "fasta",
    )
    SeqIO.write(
        [
            SeqRecord(Seq("AAAA"), id="sample1", description=""),
            SeqRecord(Seq("CCCC"), id="sample2", description=""),
        ],
        submission_dir / "sequence.fsa",
        "fasta",
    )
    # sample1 is complete, sample2 remains incomplete for every attempt.
    pd.DataFrame(
        {
            "gs-sample_name": ["sample1", "sample2"],
            "gisaid_accession_epi_isl_id": ["EPI_ISL_111", ""],
        }
    ).to_csv(tmp_path / "submission_status_report.csv", index=False)
    sleep_calls: list[int] = []

    def fake_sleep(seconds):
        sleep_calls.append(seconds)

    def fake_process_gisaid_log(**kwargs):
        return pd.DataFrame()

    monkeypatch.setattr(module.time, "sleep", fake_sleep)
    monkeypatch.setattr(module, "process_gisaid_log", fake_process_gisaid_log)

    def fake_run(command, cwd, stdout, stderr):
        Path(command[command.index("--log") + 1]).write_text("ok\n")
        return types.SimpleNamespace(returncode=0, stdout=b"", stderr=b"")

    monkeypatch.setattr(module.subprocess, "run", fake_run)

    status = module.submit_gisaid(
        organism="COV",
        submission_dir=str(submission_dir),
        submission_name="sub1",
        config_dict={"Username": "u", "Password": "p", "Client-Id": "cid"},
        submission_type="TEST",
    )
    captured = capsys.readouterr()

    assert status == "ERROR"
    assert sleep_calls == [5]
    assert captured.out == (
        "Uploading sample files to GISAID-COV, as a 'TEST' submission. "
        "If this is not intended, interrupt immediately.\n"
        "\n"
        "Submission attempt: 1\n"
        "\n"
        "Submission attempt: 2\n"
        "\n"
        "Submission attempt: 3\n"
        "\n"
        "Submission attempt: 4\n"
    )
    assert captured.err == (
        "Error: 1 sample(s) failed to upload to GISAID\n"
        f"Please check log file at: {submission_dir}/gisaid_upload_log_attempt_{{1,2,3}}.txt\n"
    )
    active_metadata = pd.read_csv(submission_dir / "metadata.csv", dtype=str)
    active_ids = [record.id for record in SeqIO.parse(submission_dir / "sequence.fsa", "fasta")]
    assert active_metadata["covv_virus_name"].tolist() == ["sample2"]
    assert active_ids == ["sample2"]


def test_process_gisaid_log__standard_message_isolate_branch_updates_status(tmp_path: Path, gisaid_handler_module, capsys):
    """Covers the normal first regex branch for EPI_ISL_* accessions."""
    module = gisaid_handler_module
    log_file = tmp_path / "gisaid.log"
    log_file.write_text('"msg": "hCoV-19/USA/CA-10/2024; EPI_ISL_101010"\n')

    returned = module.process_gisaid_log(log_file=str(log_file), submission_dir=str(tmp_path))

    assert returned.empty
    assert capsys.readouterr().out == "GISAID isolates found.\n"
    assert len(module.upload_log.update_calls) == 1
    assert_update_call(
        module.upload_log.update_calls[0],
        tmp_path,
        [
            {
                "gs-sample_name": "hCoV-19/USA/CA-10/2024",
                "gisaid_accession_epi_isl_id": "EPI_ISL_101010",
            }
        ],
    )


def test_process_gisaid_log__standard_message_segment_branch_updates_status_with_epi_digits(tmp_path: Path,
    gisaid_handler_module,
    monkeypatch: pytest.MonkeyPatch,
    capsys,
):
    """Covers the `elif re.match(r"EPI\\d+", accession_string)` branch.

    The source regex accepts `EPI12345` in the broad line matcher, but the later
    `findall` pattern only extracts `EPI_*` tokens. This test patches only the
    extraction call so the intended branch is covered and documented.
    """
    module = gisaid_handler_module
    log_file = tmp_path / "gisaid.log"
    log_file.write_text('epi_id: hCoV-19/USA/CA-11/2024; EPI12345\n')

    def fake_findall(pattern, line):
        return ["hCoV-19/USA/CA-11/2024", "EPI12345"]

    monkeypatch.setattr(module.re, "findall", fake_findall)

    with pytest.raises(KeyError, match="gisaid_accession_epi_isl_id"):
        module.process_gisaid_log(log_file=str(log_file), submission_dir=str(tmp_path))

    assert capsys.readouterr().out == "GISAID segments found.\n"
    assert len(module.upload_log.update_calls) == 1
    assert_update_call(
        module.upload_log.update_calls[0],
        tmp_path,
        [
            {
                "gs-segment_name": "hCoV-19/USA/CA-11/2024",
                "gisaid_accession_epi_id": "EPI12345",
            }
        ],
    )


def test_process_gisaid_log__existing_validation_error_isolate_and_segment_branches(tmp_path: Path, gisaid_handler_module, capsys):
    """Covers both already-exists branches: EPI_ISL_* and EPI_* accessions."""
    module = gisaid_handler_module
    log_file = tmp_path / "gisaid.log"
    log_file.write_text(
        "\n".join(
            [
                '"code": "validation_error" sample already exists; existing_virus_name: hCoV-19/USA/CA-12/2024; [\'EPI_ISL_121212\']',
                '"code": "validation_error" sample already exists; existing_virus_name: hCoV-19/USA/CA-13/2024; [\'EPI_131313\']',
            ]
        )
    )

    returned = module.process_gisaid_log(log_file=str(log_file), submission_dir=str(tmp_path))

    assert returned.empty
    assert capsys.readouterr().out == "GISAID isolates and GISAID segments found.\n"
    assert len(module.upload_log.update_calls) == 2
    assert_update_call(
        module.upload_log.update_calls[0],
        tmp_path,
        [
            {
                "gs-sample_name": "hCoV-19/USA/CA-12/2024",
                "gisaid_accession_epi_isl_id": "EPI_ISL_121212",
            }
        ],
    )

    assert_update_call(
        module.upload_log.update_calls[1],
        tmp_path,
        [
            {
                "gs-segment_name": "hCoV-19/USA/CA-13/2024",
                "gisaid_accession_epi_id": "EPI_131313",
            }
        ],
    )


@pytest.mark.parametrize(
    (
        "organism",
        "metadata_column",
        "fasta_column",
        "metadata_values",
        "fasta_ids",
        "status_frame",
        "expected_remaining_metadata",
        "expected_remaining_fasta",
    ),
    [
        (
            "FLU",
            "Isolate_Name",
            "gs-segment_name",
            ["iso1", "iso2"],
            ["seg1", "seg2"],
            pd.DataFrame(
                {
                    "gs-sample_name": ["iso1", "iso2"],
                    "gs-segment_name": ["seg1", "seg2"],
                    "gisaid_accession_epi_isl_id": ["EPI_ISL_1", ""],
                    "gisaid_accession_epi_id": ["EPI123", ""],
                }
            ),
            ["iso2"],
            ["seg2"],
        ),
        (
            "COV",
            "covv_virus_name",
            "gs-sample_name",
            ["cov1", "cov2"],
            ["cov1", "cov2"],
            pd.DataFrame(
                {
                    "gs-sample_name": ["cov1", "cov2"],
                    "gisaid_accession_epi_isl_id": ["EPI_ISL_1", ""],
                }
            ),
            ["cov2"],
            ["cov2"],
        ),
        (
            "ARBO",
            "arbo_virus_name",
            "gs-sample_name",
            ["arbo1", "arbo2"],
            ["arbo1", "arbo2"],
            pd.DataFrame(
                {
                    "gs-sample_name": ["arbo1", "arbo2"],
                    "gisaid_accession_epi_isl_id": ["EPI_ISL_1", ""],
                }
            ),
            ["arbo2"],
            ["arbo2"],
        ),
    ],
)
def test_submit_gisaid__status_filter_branches_rewrite_remaining_metadata_and_fasta(
    tmp_path: Path,
    gisaid_handler_module,
    monkeypatch: pytest.MonkeyPatch,
    organism: str,
    metadata_column: str,
    fasta_column: str,
    metadata_values: list[str],
    fasta_ids: list[str],
    status_frame: pd.DataFrame,
    expected_remaining_metadata: list[str],
    expected_remaining_fasta: list[str],
):
    """Covers the FLU, COV, and non-COV/non-FLU filtering branches in submit_gisaid()."""
    module = gisaid_handler_module
    submission_dir = tmp_path / organism / "GISAID"
    submission_dir.mkdir(parents=True)

    pd.DataFrame({metadata_column: metadata_values}).to_csv(submission_dir / "orig_metadata.csv", index=False)
    pd.DataFrame({metadata_column: metadata_values}).to_csv(submission_dir / "metadata.csv", index=False)
    SeqIO.write(
        [SeqRecord(Seq("AAAA"), id=fasta_ids[0], description=""), SeqRecord(Seq("CCCC"), id=fasta_ids[1], description="")],
        submission_dir / "orig_sequence.fsa",
        "fasta",
    )
    SeqIO.write(
        [SeqRecord(Seq("AAAA"), id=fasta_ids[0], description=""), SeqRecord(Seq("CCCC"), id=fasta_ids[1], description="")],
        submission_dir / "sequence.fsa",
        "fasta",
    )
    status_frame.to_csv(submission_dir.parent / "submission_status_report.csv", index=False)
    sleep_calls: list[int] = []

    def fake_sleep(seconds):
        sleep_calls.append(seconds)

    def fake_process_gisaid_log(**kwargs):
        return pd.DataFrame()

    monkeypatch.setattr(module.time, "sleep", fake_sleep)
    monkeypatch.setattr(module, "process_gisaid_log", fake_process_gisaid_log)

    run_calls = {"count": 0}

    def fake_run(command, cwd, stdout, stderr):
        run_calls["count"] += 1
        Path(command[command.index("--log") + 1]).write_text("ok\n")
        if run_calls["count"] == 2:
            active_metadata = pd.read_csv(submission_dir / "metadata.csv", dtype=str)
            active_ids = [record.id for record in SeqIO.parse(submission_dir / "sequence.fsa", "fasta")]
            assert active_metadata[metadata_column].tolist() == expected_remaining_metadata
            assert active_ids == expected_remaining_fasta
            raise RuntimeError(f"{organism} retry inputs were rewritten")
        return types.SimpleNamespace(returncode=0, stdout=b"", stderr=b"")

    monkeypatch.setattr(module.subprocess, "run", fake_run)

    with pytest.raises(RuntimeError, match=f"{organism} retry inputs were rewritten"):
        module.submit_gisaid(
            organism=organism,
            submission_dir=str(submission_dir),
            submission_name="sub1",
            config_dict={"Username": "u", "Password": "p", "Client-Id": "cid"},
            submission_type="TEST",
        )
