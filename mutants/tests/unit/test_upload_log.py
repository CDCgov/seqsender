from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path
from typing import Any
import os
from datetime import datetime
import pandas as pd
import pytest

def _source_root() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        mutant_src = parent / "src"
        if parent.name == "mutants" and (mutant_src / "upload_log.py").exists():
            return mutant_src

        normal_src = parent / "src"
        if (normal_src / "upload_log.py").exists() and parent.name != "mutants":
            return normal_src
    raise RuntimeError("Could not find src/upload_log.py")

SOURCE_DIR = _source_root()
MODULE_PATH = SOURCE_DIR / "upload_log.py"

"""Unit tests initially generated with ChatGPT v5.5 "Deep Research" and "Thinking"
then modified for use and clarity with validation and coverage testing via mutmut."""

#*******************************************************************************
#                         create fake schemas
#*******************************************************************************

class FakeSchemaErrors(Exception):
    pass

class FakeSchema:
    def __init__(self, name: str, *, should_raise: bool = False):
        self.name = name
        self.should_raise = should_raise
        self.calls: list[tuple[pd.DataFrame, bool]] = []

    def validate(self, df: pd.DataFrame, lazy: bool = False) -> pd.DataFrame:
        self.calls.append((df.copy(), lazy))
        if lazy is not True:
            raise AssertionError(f"{self.name} schema validate lazy must be True")
        if self.should_raise:
            raise FakeSchemaErrors(self.name)
        return df

#*******************************************************************************
#                    create test upload_log.py connections
#*******************************************************************************

@pytest.fixture()
def upload_log_module(monkeypatch: pytest.MonkeyPatch):
    """Import upload_log.py with all external dependencies replaced by stubs."""
    monkeypatch.syspath_prepend(str(SOURCE_DIR))

    src_pkg: Any = types.ModuleType("src")
    src_pkg.__path__ = [str(SOURCE_DIR)]
    monkeypatch.setitem(sys.modules, "src", src_pkg)


    def alias_src_module(name, module):
        monkeypatch.setitem(sys.modules, f"src.{name}", module)
        monkeypatch.setitem(sys.modules, name, module)
        setattr(src_pkg, name, module)

    # pandera stub: upload_log.py uses both ``from pandera import pandera`` and
    # ``pandera.errors.SchemaErrors``.
    pandera_mod: Any = types.ModuleType("pandera")
    pandera_mod.errors = types.SimpleNamespace(SchemaErrors=FakeSchemaErrors)
    pandera_mod.pandera = pandera_mod
    pandera_mod.DataFrameSchema = object
    pandera_mod.Column = object
    pandera_mod.Check = object
    pandera_mod.Index = object
    pandera_mod.MultiIndex = object
    monkeypatch.setitem(sys.modules, "pandera", pandera_mod)

    # settings constants used by upload_log.py.
    settings_mod: Any = types.ModuleType("settings")
    settings_mod.SAMPLE_NAME_DATABASE_PREFIX = {
        "BIOSAMPLE": "bs-",
        "SRA": "sra-",
        "GENBANK": "gb-",
    }
    settings_mod.BIOSAMPLE_SUBMISSION_STATUS_COLUMNS = [
        "biosample_status",
        "biosample_accession",
        "biosample_message",
    ]
    settings_mod.SRA_SUBMISSION_STATUS_COLUMNS = [
        "sra_status",
        "sra_accession",
        "sra_message",
    ]
    settings_mod.GENBANK_SUBMISSION_STATUS_COLUMNS = [
        "genbank_status",
        "genbank_accession",
        "genbank_message",
    ]
    settings_mod.SUBMISSION_LOG_COLUMNS = [
        "Submission_Name",
        "Organism",
        "Database",
        "Submission_Type",
        "Submission_Date",
        "Submission_ID",
        "Submission_Status",
        "Submission_Directory",
        "Config_File",
        "Update_Date",
    ]
    alias_src_module("settings", settings_mod)

    # config package/schema stubs.
    package_names = [
        "config",
        "config.seqsender",
        "config.seqsender.submission_status_report",
    ]
    for name in package_names:
        pkg: Any = types.ModuleType(name)
        pkg.__path__ = []  # mark as package
        monkeypatch.setitem(sys.modules, name, pkg)

    schema_modules = {
        "config.seqsender.upload_log_schema": "upload",
        "config.seqsender.submission_status_report.biosample_submission_status_report_schema": "biosample",
        "config.seqsender.submission_status_report.sra_submission_status_report_schema": "sra",
        "config.seqsender.submission_status_report.genbank_submission_status_report_schema": "genbank",
    }
    schemas: dict[str, FakeSchema] = {}
    for module_name, schema_name in schema_modules.items():
        schema_mod: Any = types.ModuleType(module_name)
        schema = FakeSchema(schema_name)
        schema_mod.schema = schema
        schemas[schema_name] = schema
        monkeypatch.setitem(sys.modules, module_name, schema_mod)

    # file_handler stub. Individual tests can monkeypatch functions on the
    # imported upload_log_module.file_handler object.
    file_handler: Any = types.ModuleType("file_handler")
    file_handler.saved_csvs = []
    file_handler.loaded_csvs = {}
    file_handler.validated_files = []
    file_handler.validated_directories = []

    def save_csv(df: pd.DataFrame, file_path: str, file_name: str | None = None, sep: str = ",") -> None:
        path = Path(file_path) / file_name if file_name else Path(file_path)
        file_handler.saved_csvs.append(
            {"df": df.copy(), "file_path": str(file_path), "file_name": file_name, "path": str(path), "sep": sep}
        )

    def load_csv(file_path: str, sep: str = ",") -> pd.DataFrame:
        return file_handler.loaded_csvs[str(file_path)].copy()

    def validate_file(file_type: str, file_path: str) -> None:
        file_handler.validated_files.append((file_type, file_path))

    def validate_directory(name: str, path: str) -> None:
        file_handler.validated_directories.append((name, path))

    file_handler.save_csv = save_csv
    file_handler.load_csv = load_csv
    file_handler.validate_file = validate_file
    file_handler.validate_directory = validate_directory
    alias_src_module("file_handler", file_handler)

    # Handler/tool stubs.
    genbank_handler: Any = types.ModuleType("genbank_handler")
    genbank_handler.calls = []

    def create_table2asn(**kwargs):
        genbank_handler.calls.append(("create_table2asn", kwargs))
        return "VALIDATED"

    def create_zip(**kwargs):
        genbank_handler.calls.append(("create_zip", kwargs))
        return None

    def process_genbank_report(**kwargs):
        genbank_handler.calls.append(("process_genbank_report", kwargs))
        return ("PROCESSED", "SUB9")

    def update_genbank_files(**kwargs):
        genbank_handler.calls.append(("update_genbank_files", kwargs))
        return None

    genbank_handler.create_table2asn = create_table2asn
    genbank_handler.create_zip = create_zip
    genbank_handler.process_genbank_report = process_genbank_report
    genbank_handler.update_genbank_files = update_genbank_files
    alias_src_module("genbank_handler", genbank_handler)

    biosample_sra_handler: Any = types.ModuleType("biosample_sra_handler")
    biosample_sra_handler.calls = []

    def process_biosample_sra_report(**kwargs):
        biosample_sra_handler.calls.append(("process_biosample_sra_report", kwargs))
        return ("PROCESSED", "SUB-BS")

    biosample_sra_handler.process_biosample_sra_report = process_biosample_sra_report
    alias_src_module("biosample_sra_handler", biosample_sra_handler)

    ncbi_handler: Any = types.ModuleType("ncbi_handler")
    ncbi_handler.calls = []

    def get_ncbi_report(*args, **kwargs):
        ncbi_handler.calls.append(("get_ncbi_report", args, kwargs))
        return "/tmp/report.xml"

    def email_table2asn(**kwargs):
        ncbi_handler.calls.append(("email_table2asn", kwargs))
        return "EMAILED"

    def submit_ncbi(**kwargs):
        ncbi_handler.calls.append(("submit_ncbi", kwargs))
        return None

    ncbi_handler.get_ncbi_report = get_ncbi_report
    ncbi_handler.email_table2asn = email_table2asn
    ncbi_handler.submit_ncbi = submit_ncbi
    alias_src_module("ncbi_handler", ncbi_handler)

    tools: Any = types.ModuleType("tools")
    tools.pretty_print_calls = []
    tools.config = {"NCBI": {"Link_Sample_Between_NCBI_Databases": False}}

    def pretty_print_pandera_errors(**kwargs):
        tools.pretty_print_calls.append(kwargs)
        return None

    def get_config(config_file, databases, decrypt_key):
        return tools.config

    tools.pretty_print_pandera_errors = pretty_print_pandera_errors
    tools.get_config = get_config
    alias_src_module("tools", tools)

    sys.modules.pop("upload_log", None)
    sys.modules.pop("src.upload_log", None)
    spec = importlib.util.spec_from_file_location("upload_log", MODULE_PATH)
    assert spec and spec.loader
    module: Any = importlib.util.module_from_spec(spec)
    sys.modules["upload_log"] = module
    setattr(src_pkg, "upload_log", module)
    spec.loader.exec_module(module)
    module._schemas = schemas
    return module

#*******************************************************************************
#                      create fake upload_log.csv file
#*******************************************************************************

@pytest.fixture()
def status_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "bs-sample_name": ["bs1"],
            "biosample_status": [""],
            "biosample_accession": [""],
            "biosample_message": [""],
            "sra-sample_name": ["sra1"],
            "sra_status": [""],
            "sra_accession": [""],
            "sra_message": [""],
            "gb-sample_name": ["gb1"],
            "genbank_status": [""],
            "genbank_accession": [""],
            "genbank_message": [""],
        }
    )

def submission_log_df(**overrides: Any) -> pd.DataFrame:
    row = {
        "Submission_Name": "sub1",
        "Organism": "FLU",
        "Database": "BIOSAMPLE",
        "Submission_Type": "TEST",
        "Submission_Date": "2024-01-01",
        "Submission_ID": "PENDING",
        "Submission_Status": "WAITING",
        "Submission_Directory": "/work/sub1/submission_files/BIOSAMPLE",
        "Config_File": "/work/sub1/config.yaml",
        "Update_Date": "2024-01-01",
    }
    row.update(overrides)
    return pd.DataFrame([row])

#*******************************************************************************
#                       create_submission_status_csv
#*******************************************************************************

def test_create_submission_status_csv__all_databases_default(upload_log_module):
    metadata = pd.DataFrame(
        {
            "bs-sample_name": ["bs1"],
            "sra-sample_name": ["sra1"],
            "gb-sample_name": ["gb1"],
        }
    )

    upload_log_module.create_submission_status_csv(database=["BIOSAMPLE", "SRA", "GENBANK"], metadata=metadata, submission_dir="/out")

    saved = upload_log_module.file_handler.saved_csvs[-1]
    assert saved["path"] == "/out/submission_status_report.csv"
    assert saved["df"].columns.tolist() == [
        "bs-sample_name",
        "biosample_status",
        "biosample_accession",
        "biosample_message",
        "sra-sample_name",
        "sra_status",
        "sra_accession",
        "sra_message",
        "gb-sample_name",
        "genbank_status",
        "genbank_accession",
        "genbank_message",
    ]
    assert saved["df"].to_dict("records") == [
        {
            "bs-sample_name": "bs1",
            "biosample_status": "",
            "biosample_accession": "",
            "biosample_message": "",
            "sra-sample_name": "sra1",
            "sra_status": "",
            "sra_accession": "",
            "sra_message": "",
            "gb-sample_name": "gb1",
            "genbank_status": "",
            "genbank_accession": "",
            "genbank_message": "",
        }
    ]
    assert saved["file_path"] == "/out/submission_status_report.csv"
    assert saved["file_name"] is None
    assert saved["sep"] == ","

#*******************************************************************************
#                     validate_submission_status_df
#*******************************************************************************

def test_validate_submission_status_df__calls_requested_schemas_with_lazy_true(upload_log_module):
    df = pd.DataFrame({"x": [1]})
    upload_log_module.validate_submission_status_df(df, ["BIOSAMPLE", "SRA", "GENBANK"])

    assert len(upload_log_module._schemas["biosample"].calls) == 1
    assert len(upload_log_module._schemas["sra"].calls) == 1
    assert len(upload_log_module._schemas["genbank"].calls) == 1

    assert upload_log_module._schemas["biosample"].calls[0][1] is True
    assert upload_log_module._schemas["sra"].calls[0][1] is True
    assert upload_log_module._schemas["genbank"].calls[0][1] is True

def test_validate_submission_status_df__sra_schema_is_called_with_lazy_true(upload_log_module):
    df = pd.DataFrame({"x": [1]})
    upload_log_module.validate_submission_status_df(df, ["SRA"])
    assert len(upload_log_module._schemas["sra"].calls) == 1
    pd.testing.assert_frame_equal(upload_log_module._schemas["sra"].calls[0][0], df)
    assert upload_log_module._schemas["sra"].calls[0][1] is True

def test_validate_submission_status_df__pretty_prints_and_exits_on_schema_error(upload_log_module, capsys):
    upload_log_module.status_report_bs_schema.should_raise = True

    with pytest.raises(SystemExit) as exc:
        upload_log_module.validate_submission_status_df(pd.DataFrame({"x": [1]}), ["BIOSAMPLE"])

    assert exc.value.code == 1
    captured = capsys.readouterr()
    assert captured.out == "   x\n0  1\n"
    assert captured.err == ""
    assert upload_log_module.tools.pretty_print_calls[-1]["file"] == "submission_status_report.csv"

#*******************************************************************************
#                      update_submission_status_csv
#*******************************************************************************

def test_update_submission_status_csv__updates_by_database_sample_name(upload_log_module, status_df):
    upload_log_module.file_handler.loaded_csvs["/logs/submission_status_report.csv"] = status_df
    update = pd.DataFrame(
        {
            "bs-sample_name": ["bs1"],
            "biosample_status": ["PROCESSED"],
            "biosample_accession": ["SAMN1"],
        }
    )

    upload_log_module.update_submission_status_csv("/logs", "BIOSAMPLE", update)

    saved = upload_log_module.file_handler.saved_csvs[-1]
    assert saved["path"] == "/logs/submission_status_report.csv"
    assert saved["df"].loc[0, "biosample_status"] == "PROCESSED"
    assert saved["df"].loc[0, "biosample_accession"] == "SAMN1"

def test_update_submission_status_csv__pops_database_directory(upload_log_module, status_df):
    upload_log_module.file_handler.loaded_csvs["/logs/submission_status_report.csv"] = status_df
    update = pd.DataFrame({"gb-sample_name": ["gb1"], "genbank_status": ["PROCESSED"]})

    upload_log_module.update_submission_status_csv("/logs/GENBANK", "GENBANK", update)

    assert upload_log_module.file_handler.validated_files[-1] == (
        "submission status report",
        "/logs/submission_status_report.csv",
    )

@pytest.mark.parametrize("database_dir", ["BIOSAMPLE", "SRA", "GENBANK"])
def test_update_submission_status_csv__pops_each_database_directory_exactly(upload_log_module, status_df, database_dir):
    upload_log_module.file_handler.loaded_csvs["/logs/submission_status_report.csv"] = status_df
    prefix = {
        "BIOSAMPLE": "bs",
        "SRA": "sra",
        "GENBANK": "gb",
    }[database_dir]

    update_database = database_dir
    update = pd.DataFrame({f"{prefix}-sample_name": [status_df[f"{prefix}-sample_name"].iloc[0]]})

    if database_dir == "BIOSAMPLE":
        update["biosample_status"] = "PROCESSED"
    elif database_dir == "SRA":
        update["sra_status"] = "PROCESSED"
    elif database_dir == "GENBANK":
        update["genbank_status"] = "PROCESSED"

    upload_log_module.update_submission_status_csv(f"/logs/{database_dir}", update_database, update)

    assert upload_log_module.file_handler.validated_files[-1] == ("submission status report", "/logs/submission_status_report.csv")
    assert upload_log_module.file_handler.saved_csvs[-1]["path"] == "/logs/submission_status_report.csv"

def test_update_submission_status_csv__detects_all_databases_and_uses_exact_set_index_args(upload_log_module, status_df, monkeypatch):
    upload_log_module.file_handler.loaded_csvs["/logs/submission_status_report.csv"] = status_df
    observed_validate_calls: list[dict[str, Any]] = []

    def fake_validate_submission_status_df(metadata, database):
        observed_validate_calls.append({"columns": metadata.columns.tolist(), "database": list(database)})

    monkeypatch.setattr(upload_log_module, "validate_submission_status_df", fake_validate_submission_status_df)
    original_set_index = pd.DataFrame.set_index
    observed_set_index: list[dict[str, Any]] = []

    def fake_set_index(self, keys, *args, **kwargs):
        observed_set_index.append({"keys": keys, "args": args, "kwargs": kwargs, "columns": self.columns.tolist()})
        return original_set_index(self, keys, *args, **kwargs)

    monkeypatch.setattr(pd.DataFrame, "set_index", fake_set_index)
    update = pd.DataFrame({"bs-sample_name": ["bs1"], "biosample_status": ["PROCESSED"]})
    upload_log_module.update_submission_status_csv("/logs", "BIOSAMPLE", update)
    assert observed_validate_calls == [
        {
            "columns": status_df.columns.tolist(),
            "database": ["BIOSAMPLE", "SRA", "GENBANK"],
        },
        {
            "columns": status_df.columns.tolist(),
            "database": ["BIOSAMPLE", "SRA", "GENBANK"],
        },
    ]
    assert observed_set_index[0]["keys"] == "bs-sample_name"
    assert observed_set_index[0]["args"] == ()
    assert observed_set_index[0]["kwargs"] == {"drop": False}
    assert observed_set_index[1]["keys"] == "bs-sample_name"
    assert observed_set_index[1]["args"] == ()
    assert observed_set_index[1]["kwargs"] == {}

def test_update_submission_status_csv__empty_update_warns_but_still_processes(upload_log_module, status_df, capsys):
    upload_log_module.file_handler.loaded_csvs["/logs/submission_status_report.csv"] = status_df

    with pytest.raises(UnboundLocalError):
        upload_log_module.update_submission_status_csv("/logs", "BIOSAMPLE", pd.DataFrame())

    assert "Error: Unable to update 'submission_status.csv' for 'BIOSAMPLE' at '/logs'. The log file may be empty.\n" == capsys.readouterr().err

#*******************************************************************************
#                         create_submission_log
#*******************************************************************************

def test_create_submission_log__creates_new_log_when_missing(upload_log_module, monkeypatch):
    def fake_isfile(path):
        return False

    monkeypatch.setattr(upload_log_module.os.path, "isfile", fake_isfile)

    upload_log_module.create_submission_log(
        database="SRA",
        organism="FLU",
        submission_name="sub1",
        submission_dir="/logs",
        database_dir="/logs/submission_files/SRA",
        config_file="/logs/config.yaml",
        submission_status="SUBMITTED",
        submission_id="PENDING",
        submission_type="TEST",
    )

    saved = upload_log_module.file_handler.saved_csvs[-1]
    assert saved["path"] == "/logs/submission_log.csv"
    assert saved["df"].loc[0, "Database"] == "SRA"
    assert saved["df"].loc[0, "Submission_Status"] == "SUBMITTED"

def test_create_submission_log__uses_exact_existing_log_path_dates_and_save_args(upload_log_module, monkeypatch):
    observed_isfile: list[str] = []

    def fake_isfile(path):
        observed_isfile.append(path)
        return False

    monkeypatch.setattr(upload_log_module.os.path, "isfile", fake_isfile)

    class FixedDateTime:
        @staticmethod
        def now():
            return datetime(2026, 6, 23)

    monkeypatch.setattr(upload_log_module, "datetime", FixedDateTime)
    upload_log_module.create_submission_log(
        database="SRA",
        organism="FLU",
        submission_name="sub1",
        submission_dir="/logs",
        database_dir="/logs/submission_files/SRA",
        config_file="/logs/config.yaml",
        submission_status="SUBMITTED",
        submission_id="PENDING",
        submission_type="TEST",
    )

    assert observed_isfile == ["/logs/submission_log.csv"]
    saved = upload_log_module.file_handler.saved_csvs[-1]
    assert saved["file_path"] == "/logs"
    assert saved["file_name"] == "submission_log.csv"
    assert saved["path"] == "/logs/submission_log.csv"
    assert saved["df"].to_dict("records") == [
        {
            "Submission_Name": "sub1",
            "Organism": "FLU",
            "Database": "SRA",
            "Submission_Type": "TEST",
            "Submission_Date": "2026-06-23",
            "Submission_ID": "PENDING",
            "Submission_Status": "SUBMITTED",
            "Submission_Directory": "/logs/submission_files/SRA",
            "Config_File": "/logs/config.yaml",
            "Update_Date": "2026-06-23",
        }
    ]

def test_create_submission_log__loads_existing_and_deduplicates(upload_log_module, monkeypatch):
    existing = submission_log_df(Submission_Status="WAITING")

    def fake_isfile(path):
        return True

    def fake_load_submission_log(submission_dir):
        return existing.copy()

    monkeypatch.setattr(upload_log_module.os.path, "isfile", fake_isfile)
    monkeypatch.setattr(upload_log_module, "load_submission_log", fake_load_submission_log)

    upload_log_module.create_submission_log(
        database="BIOSAMPLE",
        organism="FLU",
        submission_name="sub1",
        submission_dir="/logs",
        database_dir="/work/sub1/submission_files/BIOSAMPLE",
        config_file="/work/sub1/config.yaml",
        submission_status="SUBMITTED",
        submission_id="SUB1",
        submission_type="TEST",
    )

    saved_df = upload_log_module.file_handler.saved_csvs[-1]["df"]
    assert len(saved_df) == 1
    assert saved_df.loc[0, "Submission_Status"] == "SUBMITTED"
    assert saved_df.loc[0, "Submission_ID"] == "SUB1"

#*******************************************************************************
#                         update_submission_log
#*******************************************************************************

def test_update_submission_log__updates_matching_row_with_exact_update_date_and_save_args(upload_log_module, monkeypatch):
    df = submission_log_df(Update_Date="2024-01-01")

    def fake_load_submission_log(submission_dir):
        return df.copy()

    class FixedDateTime:
        @staticmethod
        def now():
            return datetime(2026, 6, 23)

    monkeypatch.setattr(upload_log_module, "load_submission_log", fake_load_submission_log)
    monkeypatch.setattr(upload_log_module, "datetime", FixedDateTime)
    upload_log_module.update_submission_log(
        database="BIOSAMPLE",
        organism="FLU",
        submission_name="sub1",
        submission_log_dir="/logs",
        submission_dir="/work/sub1/submission_files/BIOSAMPLE",
        submission_status="PROCESSED",
        submission_id="SUB99",
        submission_type="TEST",
    )

    saved = upload_log_module.file_handler.saved_csvs[-1]
    assert saved["file_path"] == "/logs"
    assert saved["file_name"] == "submission_log.csv"
    assert saved["path"] == "/logs/submission_log.csv"
    saved_df = saved["df"]
    assert saved_df.to_dict("records") == [
        {
            "Submission_Name": "sub1",
            "Organism": "FLU",
            "Database": "BIOSAMPLE",
            "Submission_Type": "TEST",
            "Submission_Date": "2024-01-01",
            "Submission_ID": "SUB99",
            "Submission_Status": "PROCESSED",
            "Submission_Directory": "/work/sub1/submission_files/BIOSAMPLE",
            "Config_File": "/work/sub1/config.yaml",
            "Update_Date": "2026-06-23",
        }
    ]

def test_update_submission_log__exits_when_row_missing(upload_log_module, monkeypatch, capsys):
    def fake_load_submission_log(submission_dir):
        return submission_log_df(Database="SRA")

    monkeypatch.setattr(upload_log_module, "load_submission_log", fake_load_submission_log)

    with pytest.raises(SystemExit) as exc:
        upload_log_module.update_submission_log(
            database="BIOSAMPLE",
            organism="FLU",
            submission_name="sub1",
            submission_log_dir="/logs",
            submission_dir="/work/sub1/submission_files/BIOSAMPLE",
            submission_status="PROCESSED",
            submission_id="SUB99",
            submission_type="TEST",
        )

    assert exc.value.code == 1
    assert "Error: 'sub1' 'BIOSAMPLE' is not present in the submission log at '/logs'.\n" == capsys.readouterr().err

#*******************************************************************************
#                         load_submission_log
#*******************************************************************************

def test_load_submission_log__drops_old_columns_uppercases_and_validates(upload_log_module):
    df = submission_log_df(
        Organism="flu",
        Database="biosample",
        Submission_Type="test",
        Submission_Status="waiting",
        Submission_ID="pending",
        Table2asn="False",
        GFF_File="x.gff",
    )
    upload_log_module.file_handler.loaded_csvs["/logs/submission_log.csv"] = df

    loaded = upload_log_module.load_submission_log("/logs")

    assert loaded.loc[0, "Organism"] == "FLU"
    assert loaded.loc[0, "Database"] == "BIOSAMPLE"
    assert upload_log_module._schemas["upload"].calls[-1][1] is True

def test_load_submission_log__uses_exact_validate_drop_deduplicate_and_uppercase_columns(upload_log_module, monkeypatch):
    df = pd.DataFrame(
        [
            {
                "Submission_Name": "sub1",
                "Organism": "flu",
                "Database": "biosample",
                "Submission_Type": "test",
                "Submission_Date": "2024-01-01",
                "Submission_ID": "pending",
                "Submission_Status": "waiting",
                "Submission_Directory": "/dir/old",
                "Config_File": "cfg",
                "Update_Date": "2024-01-01",
                "Table2asn": "old",
                "GFF_File": "old",
            },
            {
                "Submission_Name": "sub1",
                "Organism": "flu",
                "Database": "biosample",
                "Submission_Type": "test",
                "Submission_Date": "2024-01-02",
                "Submission_ID": "submitted",
                "Submission_Status": "submitted",
                "Submission_Directory": "/dir/new",
                "Config_File": "cfg",
                "Update_Date": "2024-01-02",
                "Table2asn": "old",
                "GFF_File": "old",
            },
        ]
    )
    upload_log_module.file_handler.loaded_csvs["/logs/submission_log.csv"] = df
    original_drop = pd.DataFrame.drop
    observed_drop: list[dict[str, Any]] = []

    def fake_drop(self, *args, **kwargs):
        observed_drop.append({"args": args, "kwargs": kwargs})
        return original_drop(self, *args, **kwargs)

    monkeypatch.setattr(pd.DataFrame, "drop", fake_drop)
    original_drop_duplicates = pd.DataFrame.drop_duplicates
    observed_drop_duplicates: list[dict[str, Any]] = []

    def fake_drop_duplicates(self, *args, **kwargs):
        observed_drop_duplicates.append({"args": args, "kwargs": kwargs})
        return original_drop_duplicates(self, *args, **kwargs)

    monkeypatch.setattr(pd.DataFrame, "drop_duplicates", fake_drop_duplicates)
    loaded = upload_log_module.load_submission_log("/logs")
    assert upload_log_module.file_handler.validated_files[-1] == ("submission_log", "/logs/submission_log.csv")
    assert observed_drop[0] == {
        "args": (),
        "kwargs": {
            "columns": ["Table2asn", "GFF_File"],
            "errors": "ignore",
        },
    }
    assert observed_drop_duplicates[0] == {
        "args": (),
        "kwargs": {
            "subset": [
                "Submission_Name",
                "Organism",
                "Database",
                "Submission_Type",
                "Config_File",
            ],
            "keep": "last",
            "ignore_index": True,
        },
    }
    assert loaded.to_dict("records") == [
        {
            "Submission_Name": "sub1",
            "Organism": "FLU",
            "Database": "BIOSAMPLE",
            "Submission_Type": "TEST",
            "Submission_Date": "2024-01-02",
            "Submission_ID": "SUBMITTED",
            "Submission_Status": "SUBMITTED",
            "Submission_Directory": "/dir/new",
            "Config_File": "cfg",
            "Update_Date": "2024-01-02",
        }
    ]
    assert "Table2asn" not in loaded.columns
    assert "GFF_File" not in loaded.columns
    assert "submission_type" not in loaded.columns
    assert "SUBMISSION_TYPE" not in loaded.columns
    assert "submission_status" not in loaded.columns
    assert "SUBMISSION_STATUS" not in loaded.columns
    assert "submission_id" not in loaded.columns
    assert "SUBMISSION_ID" not in loaded.columns

def test_load_submission_log__pretty_prints_and_exits_on_schema_error(upload_log_module, capsys):
    upload_log_module.upload_schema.should_raise = True
    upload_log_module.file_handler.loaded_csvs["/logs/submission_log.csv"] = submission_log_df()

    with pytest.raises(SystemExit) as exc:
        upload_log_module.load_submission_log("/logs")

    assert exc.value.code == 1
    assert "Error: Upload log columns are incorrect. Cannot process submissions.\n" == capsys.readouterr().err
    assert upload_log_module.tools.pretty_print_calls[-1]["file"] == "/logs/submission_log.csv"

#*******************************************************************************
#                         validate_fields_exist
#*******************************************************************************

def test_validate_fields_exist__checks_directory_status_report_and_config(upload_log_module):
    df = submission_log_df(
        Submission_Directory="/logs/submission_files/BIOSAMPLE",
        Config_File="/logs/config.yaml",
    )

    upload_log_module.validate_fields_exist(df)

    assert upload_log_module.file_handler.validated_directories == [
        ("directory", "/logs/submission_files/BIOSAMPLE")
    ]
    assert ("submission_status_report", "/logs/submission_files/submission_status_report.csv") in upload_log_module.file_handler.validated_files
    assert ("config file", "/logs/config.yaml") in upload_log_module.file_handler.validated_files

#*******************************************************************************
#                         process_biosample_sra
#*******************************************************************************

def test_process_biosample_sra__processed(upload_log_module):
    done, status = upload_log_module.process_biosample_sra(
        submission_name="sub1",
        database="BIOSAMPLE",
        organism="FLU",
        submission_log_dir="/logs",
        submission_dir="/dbdir",
        curr_status="PROCESSED",
        config_dict={},
        submission_type="TEST",
    )

    assert (done, status) == (True, "PROCESSED")
    assert upload_log_module.ncbi_handler.calls == []

def test_process_biosample_sra__report_updates_and_returns_done(upload_log_module, monkeypatch):
    updates = []

    def fake_update_submission_log(**kwargs):
        return updates.append(kwargs)

    monkeypatch.setattr(upload_log_module, "update_submission_log", fake_update_submission_log)

    done, status = upload_log_module.process_biosample_sra(
        submission_name="sub1",
        database="BIOSAMPLE",
        organism="FLU",
        submission_log_dir="/logs",
        submission_dir="/dbdir",
        curr_status="SUBMITTED",
        config_dict={"Username": "u"},
        submission_type="TEST",
    )

    assert (done, status) == (True, "PROCESSED")
    assert updates[-1]["submission_id"] == "SUB-BS"

def test_process_biosample_sra__no_report_keeps_current_status(upload_log_module, monkeypatch):
    def fake_get_ncbi_report(**kwargs):
        return None

    updates = []

    def fake_update_submission_log(**kwargs):
        updates.append(kwargs)

    monkeypatch.setattr(upload_log_module.ncbi_handler, "get_ncbi_report", fake_get_ncbi_report)
    monkeypatch.setattr(upload_log_module, "update_submission_log", fake_update_submission_log)

    done, status = upload_log_module.process_biosample_sra(
        submission_name="sub1",
        database="SRA",
        organism="FLU",
        submission_log_dir="/logs",
        submission_dir="/dbdir",
        curr_status="SUBMITTED",
        config_dict={},
        submission_type="TEST",
    )

    assert (done, status) == (False, "SUBMITTED")
    assert updates[-1]["submission_id"] == "PENDING"

#*******************************************************************************
#                         upload_log_submit_genbank
#*******************************************************************************

def test_upload_log_submit_genbank__table2asn_validated_emails(upload_log_module, monkeypatch):
    updates = []

    def fake_update_submission_log(**kwargs):
        updates.append(kwargs)

    monkeypatch.setattr(upload_log_module, "update_submission_log", fake_update_submission_log)

    status = upload_log_module.upload_log_submit_genbank(
        genbank_type="GENBANK-TBL2ASN",
        submission_name="sub1",
        organism="FLU",
        submission_log_dir="/logs",
        submission_dir="/gb",
        config_dict={"Description": {}},
        submission_type="TEST",
    )

    assert status == "EMAILED"
    assert ("create_table2asn", {"submission_name": "sub1", "submission_dir": "/gb"}) in upload_log_module.genbank_handler.calls
    assert updates[-1]["submission_status"] == "EMAILED"
    assert updates[-1]["submission_id"] == "VALIDATED"

def test_upload_log_submit_genbank__table2asn_invalid_stays_pending(upload_log_module, monkeypatch):
    def fake_create_table2asn(**kwargs):
        return "ERROR"

    updates = []

    def fake_update_submission_log(**kwargs):
        updates.append(kwargs)

    monkeypatch.setattr(upload_log_module.genbank_handler, "create_table2asn", fake_create_table2asn)
    monkeypatch.setattr(upload_log_module, "update_submission_log", fake_update_submission_log)

    status = upload_log_module.upload_log_submit_genbank(
        genbank_type="GENBANK-TBL2ASN",
        submission_name="sub1",
        organism="FLU",
        submission_log_dir="/logs",
        submission_dir="/gb",
        config_dict={},
        submission_type="TEST",
    )

    assert status == "PENDING"
    assert updates[-1]["submission_status"] == "PENDING"
    assert updates[-1]["submission_id"] == "ERROR"

def test_upload_log_submit_genbank__ftp_zips_and_submits(upload_log_module, monkeypatch):
    updates = []

    def fake_update_submission_log(**kwargs):
        updates.append(kwargs)

    monkeypatch.setattr(upload_log_module, "update_submission_log", fake_update_submission_log)

    status = upload_log_module.upload_log_submit_genbank(
        genbank_type="GENBANK-FTP",
        submission_name="sub1",
        organism="COV",
        submission_log_dir="/logs",
        submission_dir="/gb",
        config_dict={},
        submission_type="PRODUCTION",
    )

    assert status == "SUBMITTED"
    assert upload_log_module.genbank_handler.calls == [
        (
            "create_zip",
            {
                "submission_name": "sub1",
                "submission_dir": "/gb",
            },
        )
    ]
    assert upload_log_module.ncbi_handler.calls == [
        (
            "submit_ncbi",
            {
                "database": "GENBANK",
                "submission_name": "sub1",
                "submission_dir": "/gb",
                "config_dict": {},
                "submission_type": "PRODUCTION",
            },
        )
    ]
    assert updates[-1]["submission_status"] == "SUBMITTED"
    assert updates[-1]["submission_id"] == "PENDING"

def test_upload_log_submit_genbank__invalid_type_exits(upload_log_module, capsys):
    with pytest.raises(SystemExit) as exc:
        upload_log_module.upload_log_submit_genbank(
            genbank_type="GENBANK-BAD",
            submission_name="sub1",
            organism="FLU",
            submission_log_dir="/logs",
            submission_dir="/gb",
            config_dict={},
            submission_type="TEST",
        )

    assert exc.value.code == 1
    assert "Error: GENBANK-BAD is not a valid GenBank submission option.\n" == capsys.readouterr().err

#*******************************************************************************
#                         process_genbank
#*******************************************************************************

def test_process_genbank__processed_or_emailed(upload_log_module):
    for curr_status in ["PROCESSED", "EMAILED"]:
        assert upload_log_module.process_genbank(
            genbank_type="GENBANK-FTP",
            submission_name="sub1",
            submission_log_dir="/logs",
            submission_dir="/gb",
            curr_status=curr_status,
            organism="FLU",
            config_dict={"Link_Sample_Between_NCBI_Databases": True},
            submission_type="TEST",
            linking_databases={"BIOSAMPLE": True, "SRA": True},
        ) == (True, curr_status)

def test_process_genbank__uppercase_emailed_is_terminal_exactly(upload_log_module):
    result = upload_log_module.process_genbank(
        genbank_type="GENBANK-TBL2ASN",
        submission_name="sub1",
        submission_log_dir="/logs",
        submission_dir="/gb",
        curr_status="EMAILED",
        organism="FLU",
        config_dict={"Link_Sample_Between_NCBI_Databases": False},
        submission_type="TEST",
        linking_databases={"BIOSAMPLE": False, "SRA": False},
    )
    assert result == (True, "EMAILED")
    assert upload_log_module.genbank_handler.calls == []
    assert upload_log_module.ncbi_handler.calls == []

def test_process_genbank__waiting_calls_submission_ready_with_exact_genbank_database(upload_log_module, monkeypatch):
    ready_calls: list[dict[str, Any]] = []

    def fake_submission_ready(**kwargs):
        ready_calls.append(kwargs)
        return False

    monkeypatch.setattr(upload_log_module, "submission_ready", fake_submission_ready)
    result = upload_log_module.process_genbank(
        genbank_type="GENBANK-FTP",
        submission_name="sub1",
        submission_log_dir="/logs",
        submission_dir="/gb",
        curr_status="WAITING",
        organism="FLU",
        config_dict={"Link_Sample_Between_NCBI_Databases": False},
        submission_type="TEST",
        linking_databases={"BIOSAMPLE": True, "SRA": True},
    )
    assert result == (False, "WAITING")
    assert ready_calls == [
        {
            "submission_requirements": {"BIOSAMPLE": True, "SRA": True},
            "config_dict": {"Link_Sample_Between_NCBI_Databases": False},
            "database": "GENBANK",
        }
    ]

def test_process_genbank__waiting_ready_updates_files_then_submits(upload_log_module, monkeypatch):
    def fake_submission_ready(**kwargs):
        return True

    def fake_upload_log_submit_genbank(**kwargs):
        return "SUBMITTED"

    monkeypatch.setattr(upload_log_module, "submission_ready", fake_submission_ready)
    monkeypatch.setattr(upload_log_module, "upload_log_submit_genbank", fake_upload_log_submit_genbank)

    done, status = upload_log_module.process_genbank(
        genbank_type="GENBANK-FTP",
        submission_name="sub1",
        submission_log_dir="/logs",
        submission_dir="/gb",
        curr_status="WAITING",
        organism="FLU",
        config_dict={"Link_Sample_Between_NCBI_Databases": True},
        submission_type="TEST",
        linking_databases={"BIOSAMPLE": True, "SRA": True},
    )

    assert (done, status) == (False, "SUBMITTED")
    assert upload_log_module.genbank_handler.calls[-1][0] == "update_genbank_files"

def test_process_genbank__not_ready_stays_waiting(upload_log_module, monkeypatch):
    def fake_submission_ready(**kwargs):
        return False

    monkeypatch.setattr(upload_log_module, "submission_ready", fake_submission_ready)

    assert upload_log_module.process_genbank(
        genbank_type="GENBANK-FTP",
        submission_name="sub1",
        submission_log_dir="/logs",
        submission_dir="/gb",
        curr_status="WAITING",
        organism="FLU",
        config_dict={"Link_Sample_Between_NCBI_Databases": False},
        submission_type="TEST",
        linking_databases={"BIOSAMPLE": False, "SRA": True},
    ) == (False, "WAITING")

def test_process_genbank__pending_non_tbl2asn_checks_report_not_table2asn(upload_log_module, monkeypatch):
    updates: list[dict[str, Any]] = []

    def fake_update_submission_log(**kwargs):
        updates.append(kwargs)

    monkeypatch.setattr(upload_log_module, "update_submission_log", fake_update_submission_log)
    result = upload_log_module.process_genbank(
        genbank_type="GENBANK-FTP",
        submission_name="sub1",
        submission_log_dir="/logs",
        submission_dir="/gb",
        curr_status="PENDING",
        organism="FLU",
        config_dict={"Link_Sample_Between_NCBI_Databases": False},
        submission_type="TEST",
        linking_databases={"BIOSAMPLE": True, "SRA": True},
    )
    assert result == (True, "PROCESSED")
    assert upload_log_module.genbank_handler.calls == [
        (
            "process_genbank_report",
            {
                "report_file": "/tmp/report.xml",
                "submission_dir": "/gb",
            },
        )
    ]
    assert upload_log_module.ncbi_handler.calls == [
        (
            "get_ncbi_report",
            ("GENBANK", "sub1", "/gb", {"Link_Sample_Between_NCBI_Databases": False}, "TEST"),
            {},
        )
    ]
    assert updates[-1]["submission_status"] == "PROCESSED"
    assert updates[-1]["submission_id"] == "SUB9"

def test_process_genbank__pending_table2asn_validated_emails(upload_log_module, monkeypatch):
    updates = []

    def fake_update_submission_log(**kwargs):
        return updates.append(kwargs)

    monkeypatch.setattr(upload_log_module, "update_submission_log", fake_update_submission_log)
    done, status = upload_log_module.process_genbank(
        genbank_type="GENBANK-TBL2ASN",
        submission_name="sub1",
        submission_log_dir="/logs",
        submission_dir="/gb",
        curr_status="PENDING",
        organism="FLU",
        config_dict={"Link_Sample_Between_NCBI_Databases": False},
        submission_type="TEST",
        linking_databases={"BIOSAMPLE": True, "SRA": True},
    )
    assert (done, status) == (False, "EMAILED")
    assert updates[-1]["submission_id"] == "VALIDATED"

def test_process_genbank__pending_table2asn_invalid_returns_exact_pending_and_updates_log(upload_log_module, monkeypatch):
    updates: list[dict[str, Any]] = []

    def fake_create_table2asn(**kwargs):
        upload_log_module.genbank_handler.calls.append(("create_table2asn_invalid", kwargs))
        return "ERROR"

    def fake_update_submission_log(**kwargs):
        updates.append(kwargs)

    monkeypatch.setattr(upload_log_module.genbank_handler, "create_table2asn", fake_create_table2asn)
    monkeypatch.setattr(upload_log_module, "update_submission_log", fake_update_submission_log)
    result = upload_log_module.process_genbank(
        genbank_type="GENBANK-TBL2ASN",
        submission_name="sub1",
        submission_log_dir="/logs",
        submission_dir="/gb",
        curr_status="PENDING",
        organism="FLU",
        config_dict={"Link_Sample_Between_NCBI_Databases": False},
        submission_type="TEST",
        linking_databases={"BIOSAMPLE": True, "SRA": True},
    )
    assert result == (False, "PENDING")
    assert upload_log_module.ncbi_handler.calls == []
    assert upload_log_module.genbank_handler.calls == [
        (
            "create_table2asn_invalid",
            {
                "submission_name": "sub1",
                "submission_dir": "/gb",
            },
        )
    ]
    assert updates == [
        {
            "database": "GENBANK-TBL2ASN",
            "organism": "FLU",
            "submission_name": "sub1",
            "submission_log_dir": "/logs",
            "submission_dir": "/gb",
            "submission_status": "PENDING",
            "submission_id": "ERROR",
            "submission_type": "TEST",
        }
    ]

def test_process_genbank__ftp_report_updates_and_returns_processed(upload_log_module, monkeypatch):
    updates = []

    def fake_update_submission_log(**kwargs):
        return updates.append(kwargs)

    monkeypatch.setattr(upload_log_module, "update_submission_log", fake_update_submission_log)

    done, status = upload_log_module.process_genbank(
        genbank_type="GENBANK-FTP",
        submission_name="sub1",
        submission_log_dir="/logs",
        submission_dir="/gb",
        curr_status="SUBMITTED",
        organism="FLU",
        config_dict={"Link_Sample_Between_NCBI_Databases": False},
        submission_type="TEST",
        linking_databases={"BIOSAMPLE": True, "SRA": True},
    )

    assert (done, status) == (True, "PROCESSED")
    assert updates[-1]["submission_id"] == "SUB9"

def test_process_genbank__no_report_keeps_status(upload_log_module, monkeypatch):
    def fake_get_ncbi_report(*args, **kwargs):
        return None

    monkeypatch.setattr(upload_log_module.ncbi_handler, "get_ncbi_report", fake_get_ncbi_report)
    updates = []

    def fake_update_submission_log(**kwargs):
        return updates.append(kwargs)

    monkeypatch.setattr(upload_log_module, "update_submission_log", fake_update_submission_log)

    done, status = upload_log_module.process_genbank(
        genbank_type="GENBANK-FTP",
        submission_name="sub1",
        submission_log_dir="/logs",
        submission_dir="/gb",
        curr_status="SUBMITTED",
        organism="FLU",
        config_dict={"Link_Sample_Between_NCBI_Databases": False},
        submission_type="TEST",
        linking_databases={"BIOSAMPLE": True, "SRA": True},
    )

    assert (done, status) == (False, "SUBMITTED")
    assert updates[-1]["submission_id"] == "PENDING"

#*******************************************************************************
#                    create_submission_requirements_dict
#*******************************************************************************

def test_create_submission_requirements_dict__handles_present_and_absent_databases(upload_log_module):
    group = pd.DataFrame(
        {
            "Database": ["BIOSAMPLE", "SRA", "GENBANK-TBL2ASN"],
            "Submission_Status": ["PROCESSED", "WAITING", "PROCESSED"],
        }
    )

    assert upload_log_module.create_submission_requirements_dict(group) == {
        "BIOSAMPLE": True,
        "SRA": False,
        "GENBANK": True,
    }

def test_create_submission_requirements_dict__genbank_ftp_not_processed(upload_log_module):
    group = pd.DataFrame({"Database": ["GENBANK-FTP"], "Submission_Status": ["SUBMITTED"]})
    assert upload_log_module.create_submission_requirements_dict(group)["GENBANK"] is False

def test_create_submission_requirements_dict__genbank_ftp_processed_exact_key_and_value(upload_log_module):
    group = pd.DataFrame({"Database": ["GENBANK-FTP"], "Submission_Status": ["PROCESSED"]})
    assert upload_log_module.create_submission_requirements_dict(group) == {
        "BIOSAMPLE": True,
        "SRA": True,
        "GENBANK": True,
    }


def test_create_submission_requirements_dict__genbank_tbl2asn_processed_and_failed_exact_key_and_value(upload_log_module):
    processed = pd.DataFrame({"Database": ["GENBANK-TBL2ASN"], "Submission_Status": ["PROCESSED"]})
    failed = pd.DataFrame({"Database": ["GENBANK-TBL2ASN"], "Submission_Status": ["PENDING"]})
    assert upload_log_module.create_submission_requirements_dict(processed) == {
        "BIOSAMPLE": True,
        "SRA": True,
        "GENBANK": True,
    }
    assert upload_log_module.create_submission_requirements_dict(failed) == {
        "BIOSAMPLE": True,
        "SRA": True,
        "GENBANK": False,
    }


def test_create_submission_requirements_dict__no_genbank_exact_default_key_and_value(upload_log_module):
    group = pd.DataFrame({"Database": ["BIOSAMPLE"], "Submission_Status": ["PROCESSED"]})
    assert upload_log_module.create_submission_requirements_dict(group) == {
        "BIOSAMPLE": True,
        "SRA": True,
        "GENBANK": True,
    }

#*******************************************************************************
#                         update_grouped_submission
#*******************************************************************************

def test_update_grouped_submission__processes_all_databases(upload_log_module, monkeypatch, capsys):
    group = pd.DataFrame(
        {
            "Submission_Name": ["sub1", "sub1", "sub1"],
            "Organism": ["FLU"] * 3,
            "Submission_Type": ["TEST"] * 3,
            "Config_File": ["/logs/config.yaml"] * 3,
            "Database": ["BIOSAMPLE", "SRA", "GENBANK-FTP"],
            "Submission_Status": ["BS-STATUS", "SRA-STATUS", "GB-STATUS"],
            "Submission_Directory": ["/bs", "/sra", "/gb"],
        }
    )
    calls: list[tuple[str, dict[str, Any]]] = []
    monkeypatch.setattr(upload_log_module, "validate_fields_exist", lambda df: calls.append(("validate", {"df": df.copy()})))

    def fake_process_biosample_sra(**kwargs):
        calls.append(("biosample_sra", kwargs))
        return (True, f"{kwargs['database']}-DONE")

    def fake_process_genbank(**kwargs):
        calls.append(("genbank", kwargs))
        return (True, "GENBANK-DONE")

    monkeypatch.setattr(upload_log_module, "process_biosample_sra", fake_process_biosample_sra)
    monkeypatch.setattr(upload_log_module, "process_genbank", fake_process_genbank)
    upload_log_module.tools.config = {"NCBI": {"Link_Sample_Between_NCBI_Databases": False}}
    upload_log_module.update_grouped_submission(group, "/logs", "test-key")
    assert [name for name, _ in calls] == [
        "validate",
        "biosample_sra",
        "biosample_sra",
        "genbank",
    ]
    assert calls[1][1] == {
        "submission_name": "sub1",
        "organism": "FLU",
        "database": "BIOSAMPLE",
        "curr_status": "BS-STATUS",
        "submission_log_dir": "/logs",
        "submission_dir": "/bs",
        "config_dict": {"Link_Sample_Between_NCBI_Databases": False},
        "submission_type": "TEST",
    }
    assert calls[2][1] == {
        "submission_name": "sub1",
        "organism": "FLU",
        "database": "SRA",
        "curr_status": "SRA-STATUS",
        "submission_log_dir": "/logs",
        "submission_dir": "/sra",
        "config_dict": {"Link_Sample_Between_NCBI_Databases": False},
        "submission_type": "TEST",
    }
    assert calls[3][1]["genbank_type"] == "GENBANK-FTP"
    assert calls[3][1]["submission_name"] == "sub1"
    assert calls[3][1]["submission_log_dir"] == "/logs"
    assert calls[3][1]["submission_dir"] == "/gb"
    assert calls[3][1]["curr_status"] == "GB-STATUS"
    assert calls[3][1]["organism"] == "FLU"
    assert calls[3][1]["config_dict"] == {"Link_Sample_Between_NCBI_Databases": False}
    assert calls[3][1]["submission_type"] == "TEST"
    assert calls[3][1]["linking_databases"] == {
        "BIOSAMPLE": True,
        "SRA": True,
        "GENBANK": True,
    }
    assert capsys.readouterr().out == (
        "\tBioSample: BIOSAMPLE-DONE\n"
        "\tSRA: SRA-DONE\n"
        "\tGenBank: GENBANK-DONE\n"
    )

def test_update_grouped_submission__invalid_genbank_option_exits(upload_log_module, monkeypatch, capsys):
    group = pd.DataFrame(
        {
            "Submission_Name": ["sub1"],
            "Organism": ["FLU"],
            "Submission_Type": ["TEST"],
            "Config_File": ["/logs/config.yaml"],
            "Database": ["GENBANK"],
            "Submission_Status": ["WAITING"],
            "Submission_Directory": ["/gb"],
        }
    )

    def fake_validate_fields_exist(df):
        return None

    monkeypatch.setattr(upload_log_module, "validate_fields_exist", fake_validate_fields_exist)

    with pytest.raises(SystemExit) as exc:
        upload_log_module.update_grouped_submission(group, "/logs", "test-key")

    assert exc.value.code == 1
    assert "Error: Incorrect database option for GenBank in 'submission_log.csv' databases '['GENBANK']' for 'sub1'.\n" == capsys.readouterr().out

def test_update_submission_status__updates_incomplete_groups_and_named_completed_group(upload_log_module, monkeypatch, capsys):
    df = pd.DataFrame(
        {
            "Submission_Name": ["sub1", "sub2"],
            "Organism": ["FLU", "FLU"],
            "Submission_Type": ["TEST", "TEST"],
            "Config_File": ["cfg", "cfg"],
            "Database": ["BIOSAMPLE", "BIOSAMPLE"],
            "Submission_Status": ["WAITING", "PROCESSED"],
        }
    )

    def fake_load_submission_log(submission_dir):
        return df.copy()

    calls = []

    def fake_update_grouped_submission(group_df, submission_log_dir, decrypt_key):
        calls.append(group_df["Submission_Name"].iloc[0])

    monkeypatch.setattr(upload_log_module, "load_submission_log", fake_load_submission_log)
    monkeypatch.setattr(upload_log_module, "update_grouped_submission", fake_update_grouped_submission)

    upload_log_module.update_submission_status("/logs", submission_name=None, decrypt_key={"test-key"})
    assert calls == ["sub1"]

    upload_log_module.update_submission_status("/logs", submission_name="sub2", decrypt_key={"test-key"})
    assert calls == ["sub1", "sub2"]
    assert capsys.readouterr().out == (
        "Checking Submissions:\n"
        "Submission: sub1\n"
        "\n"
        "Updating submissions complete.\n"
        "Checking Submissions:\n"
        "Submission: sub2\n"
        "\n"
        "Updating submissions complete.\n"
    )

def test_update_submission_status__skips_groups_that_are_all_processed_or_emailed(upload_log_module, monkeypatch, capsys):
    df = pd.DataFrame(
        {
            "Submission_Name": ["done1", "done1", "todo1"],
            "Organism": ["FLU", "FLU", "FLU"],
            "Submission_Type": ["TEST", "TEST", "TEST"],
            "Config_File": ["cfg", "cfg", "cfg"],
            "Database": ["GENBANK-TBL2ASN", "SRA", "BIOSAMPLE"],
            "Submission_Status": ["EMAILED", "PROCESSED", "WAITING"],
        }
    )
    monkeypatch.setattr(upload_log_module, "load_submission_log", lambda submission_dir: df.copy())
    calls: list[str] = []

    def fake_update_grouped_submission(group_df, submission_log_dir, decrypt_key):
        calls.append(group_df["Submission_Name"].iloc[0])

    monkeypatch.setattr(upload_log_module, "update_grouped_submission", fake_update_grouped_submission)
    upload_log_module.update_submission_status("/logs", submission_name=None, decrypt_key={"test-key"})
    assert calls == ["todo1"]
    assert capsys.readouterr().out == (
        "Checking Submissions:\n"
        "Submission: todo1\n"
        "\n"
        "Updating submissions complete.\n"
    )
