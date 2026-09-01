# tests/unit/test_tools.py
from __future__ import annotations
import pandas as pd
import pytest
import builtins
import importlib
import sys
import os
import types
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace
from typing import Any
from unittest.mock import Mock, call

def _source_root() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        mutant_src = parent / "src"
        if parent.name == "mutants" and (mutant_src / "tools.py").exists():
            return mutant_src

        normal_src = parent / "src"
        if (normal_src / "tools.py").exists() and parent.name != "mutants":
            return normal_src

    raise RuntimeError("Could not find src/tools.py")


SOURCE_DIR = _source_root()
MODULE_PATH = SOURCE_DIR / "tools.py"
if str(SOURCE_DIR) not in sys.path:
    sys.path.insert(0, str(SOURCE_DIR))

"""Unit tests initially generated with ChatGPT v5.5 "Deep Research" and "Thinking"
then modified for use and clarity with validation and coverage testing via mutmut."""

#*******************************************************************************
#                     create fake schemas / configs
#*******************************************************************************

class _ImportTimeSeqSenderSchema:
    def __init__(self) -> None:
        self.update_columns_calls: list[Any] = []
        self.validate_calls: list[Any] = []

    def update_columns(self, payload: Any) -> None:
        self.update_columns_calls.append(payload)

    def validate(self, df: pd.DataFrame, lazy: bool = True) -> pd.DataFrame:
        self.validate_calls.append((df.copy(), lazy))
        return df

class _ImportTimeSchemaErrors(Exception):
    pass


class _ImportTimeCheck:
    @staticmethod
    def str_matches(*args: Any, **kwargs: Any) -> tuple[str, tuple[Any, ...], dict[str, Any]]:
        return ("str_matches", args, kwargs)

    @staticmethod
    def isin(*args: Any, **kwargs: Any) -> tuple[str, tuple[Any, ...], dict[str, Any]]:
        return ("isin", args, kwargs)

    @staticmethod
    def str_length(*args: Any, **kwargs: Any) -> tuple[str, tuple[Any, ...], dict[str, Any]]:
        return ("str_length", args, kwargs)

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        self.args = args
        self.kwargs = kwargs


class _ImportTimeValidator:
    errors: dict[str, Any] = {}

    def __init__(self, schema: dict[str, Any]) -> None:
        self.schema = schema

    def validate(self, config_dict: dict[str, Any], schema: dict[str, Any]) -> bool:
        return True

def _install_import_stubs() -> None:
    settings_stub: Any = types.ModuleType("settings")
    settings_stub.PROG_DIR = str(SOURCE_DIR)
    settings_stub.SCHEMA_EXCLUSIONS = [
        "config.seqsender.upload_log_schema",
        "config_file.ncbi_schema",
        "config_file.ncbi_gisaid_schema",
        "config_file.gisaid_schema",
    ]
    settings_stub.BIOSAMPLE_REGEX = r"^bs-|^bioproject$|^organism$|^collection_date$"
    settings_stub.SRA_REGEX = r"^sra-|^bioproject$|bs-sample_name|^organism$|^collection_date$"
    settings_stub.GISAID_REGEX = r"^gs-|^collection_date$|^authors$"
    settings_stub.GENBANK_REGEX = r"^gb-sample_name$"
    settings_stub.GENBANK_REGEX_CMT = r"^gb-sample_name$|^cmt-"
    settings_stub.GENBANK_REGEX_SRC = r"^gb-sample_name$|^src-|^bioproject$|^organism$|^collection_date$"
    settings_stub.GENBANK_DEPRECATED_COLUMNS = [
        "src-Authority",
        "src-Biotype",
        "src-Biovar",
        "src-Chemovar",
        "src-Forma",
        "src-Forma_specialis",
        "src-Identified_by",
        "src-Pathovar",
        "src-Pop_variant",
        "src-Serogroup",
        "src-Subclone",
        "src-Subtype",
        "src-Substrain",
        "src-Type",
    ]
    sys.modules["settings"] = settings_stub
    sys.modules["src.settings"] = settings_stub

    ncbi_handler_stub: Any = types.ModuleType("src.ncbi_handler")
    sys.modules["src.ncbi_handler"] = ncbi_handler_stub
    sys.modules["ncbi_handler"] = ncbi_handler_stub

    file_handler_stub: Any = types.ModuleType("file_handler")
    file_handler_stub.load_yaml_calls = []

    def load_yaml(yaml_type, yaml_path):
        file_handler_stub.load_yaml_calls.append(
            {
                "yaml_type": yaml_type,
                "yaml_path": yaml_path,
            }
        )
        return {}

    def load_csv(file_path, sep=","):
        return pd.DataFrame()

    def load_fasta_file(fasta_file):
        return pd.DataFrame()

    def save_yaml(config_dict, yaml_path):
        return None

    file_handler_stub.load_yaml = load_yaml
    file_handler_stub.load_csv = load_csv
    file_handler_stub.load_fasta_file = load_fasta_file
    file_handler_stub.save_yaml = save_yaml
    sys.modules["file_handler"] = file_handler_stub
    sys.modules["src.file_handler"] = file_handler_stub

    pandera_stub: Any = types.ModuleType("pandera")
    pandera_stub.Check = _ImportTimeCheck
    pandera_stub.errors = SimpleNamespace(SchemaErrors=_ImportTimeSchemaErrors)
    sys.modules["pandera"] = pandera_stub

    cerberus_stub: Any = types.ModuleType("cerberus")
    cerberus_stub.Validator = _ImportTimeValidator
    sys.modules["cerberus"] = cerberus_stub

    config_stub: Any = types.ModuleType("config")
    config_stub.__path__ = []  # mark as package for nested imports
    seqsender_pkg_stub: Any = types.ModuleType("config.seqsender")
    seqsender_pkg_stub.__path__ = []
    seqsender_schema_module: Any = types.ModuleType("config.seqsender.seqsender_schema")
    seqsender_schema_module.schema = _ImportTimeSeqSenderSchema()
    sys.modules["config"] = config_stub
    sys.modules["config.seqsender"] = seqsender_pkg_stub
    sys.modules["config.seqsender.seqsender_schema"] = seqsender_schema_module

def _load_tools_module():
    _install_import_stubs()
    sys.modules.pop("src.tools", None)
    sys.modules.pop("tools", None)
    spec = importlib.util.spec_from_file_location("tools", MODULE_PATH)
    assert spec and spec.loader
    module: Any = importlib.util.module_from_spec(spec)
    sys.modules["tools"] = module
    spec.loader.exec_module(module)
    return module

tools = _load_tools_module()

@pytest.fixture
def base_config() -> dict[str, Any]:
    return {
        "Submission": {
            "NCBI": {
                "Submission_Position": 1,
                "Specified_Release_Date": None,
                "BioSample_Package": "Pathogen.cl.1.0",
            },
            "GISAID": {
                "Submission_Position": 2,
            },
        }
    }


@pytest.fixture
def configurable_schema() -> dict[str, Any]:
    return {
        "Submission": {
            "schema": {
                "NCBI": {
                    "schema": {
                        "BioSample_Package": {"required": False, "nullable": True},
                        "Publication_Title": {"required": False, "nullable": True},
                        "Publication_Status": {"required": False, "nullable": True},
                    }
                }
            }
        }
    }


class FakeValidator:
    """Small Cerberus Validator double used by get_config tests."""

    errors = {"field": ["bad"]}
    result = True
    init_calls: list[dict[str, Any]] = []
    calls: list[tuple[Any, Any]] = []

    def __init__(self, schema: dict[str, Any]) -> None:
        type(self).init_calls.append({"schema": schema})
        self.schema = schema

    def validate(self, config_dict: dict[str, Any], schema: dict[str, Any]) -> bool:
        type(self).calls.append((config_dict, schema))
        return type(self).result


class DummySchema:
    def __init__(self, name: str = "schema", raise_on_validate: Exception | None = None) -> None:
        self.name = name
        self.raise_on_validate = raise_on_validate
        self.validate_calls: list[tuple[pd.DataFrame, bool]] = []
        self.update_columns_calls: list[Any] = []
        self.columns: dict[str, Any] = {}

    def update_columns(self, payload: Any) -> None:
        self.update_columns_calls.append(payload)

    def validate(self, df: pd.DataFrame, lazy: bool = True) -> pd.DataFrame:
        self.validate_calls.append((df.copy(), lazy))
        if lazy is not True:
            raise AssertionError(f"{self.name} validate lazy must be True")
        if self.raise_on_validate is not None:
            raise self.raise_on_validate
        return df


class FakeSchemaErrors(Exception):
    pass


class FixedTimestamp:
    @staticmethod
    def now() -> pd.Timestamp:
        return pd._libs.tslibs.timestamps.Timestamp("2026-05-29")


@dataclass
class FakeSchemaColumn:
    required: bool
    description: str | None

# sys.modules.setdefault("config", types.ModuleType("config"))
# sys.modules.setdefault("config.seqsender", types.ModuleType("config.seqsender"))
# _seqsender_schema_module = types.ModuleType("config.seqsender.seqsender_schema")
# _seqsender_schema_module.schema = _ImportTimeSeqSenderSchema()
# sys.modules.setdefault("config.seqsender.seqsender_schema", _seqsender_schema_module)


#*******************************************************************************
#                        determine_parent_database
#*******************************************************************************


@pytest.mark.parametrize(
    ("databases", "expected"),
    [
        (["BIOSAMPLE"], {"ncbi"}),
        (["SRA"], {"ncbi"}),
        (["GENBANK"], {"ncbi"}),
        (["GISAID"], {"gisaid"}),
        (["BIOSAMPLE", "SRA", "GENBANK"], {"ncbi"}),
        (["GENBANK", "GISAID"], {"ncbi", "gisaid"}),
    ],
)
def test_determine_parent_database__maps_databases_to_parent_portals(databases: list[str], expected: set[str]) -> None:
    assert tools.determine_parent_database(databases) == expected

def test_determine_parent_database__empty_or_unknown_database_exits(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc:
        tools.determine_parent_database([])
    assert exc.value.code == 1
    assert capsys.readouterr().err == "Error: Submission portals list cannot be empty.\n"

def test_determine_parent_database__unknown_database_exits(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc:
        tools.determine_parent_database(["BADDB"])
    assert exc.value.code == 1
    assert capsys.readouterr().err == "Error: Submission portals list cannot be empty.\n"

#*******************************************************************************
#                          decrypt_passwords
#*******************************************************************************

def test_decrypt_passwords__decrypts_ncbi_password() -> None:
    key = tools.Fernet.generate_key()
    encrypted_password = tools.Fernet(key).encrypt(b"secret")
    config = {"Submission": {"NCBI": {"Password": encrypted_password}}}
    result = tools.decrypt_passwords(config_dict=config, submission_portals={"ncbi"}, key=key.decode())
    assert result["Submission"]["NCBI"]["Password"] == b"secret"

def test_decrypt_passwords__decrypts_gisaid_password_and_client_id() -> None:
    key = tools.Fernet.generate_key()
    encrypter = tools.Fernet(key)
    config = {
        "Submission": { "GISAID": {
            "Password": encrypter.encrypt(b"secret"),
            "Client-Id": encrypter.encrypt(b"client-id"),
    }}}
    result = tools.decrypt_passwords(config_dict=config, submission_portals={"gisaid"}, key=key.decode())
    assert result["Submission"]["GISAID"]["Password"] == b"secret"
    assert result["Submission"]["GISAID"]["Client-Id"] == b"client-id"

def test_decrypt_passwords__invalid_unencrypted_password_prints_helpful_error(capsys: pytest.CaptureFixture[str]) -> None:
    key = tools.Fernet.generate_key()
    config = {"Submission": {"NCBI": {"Password": "plain-text-password"}}}
    with pytest.raises(tools.InvalidToken):
        tools.decrypt_passwords(config_dict=config, submission_portals={"NCBI"}, key=key.decode())
    assert capsys.readouterr().err == (
        "Passwords field does not appear to be encrypted. "
        "Use SeqSender command 'load_credentials' to encrypt your credentials before submission.\n"
    )

def test_decrypt_passwords__invalid_encrypted_looking_password_does_not_print_plaintext_warning(capsys: pytest.CaptureFixture[str]) -> None:
    key = tools.Fernet.generate_key()
    config = {"Submission": {"NCBI": {"Password": "not-a-valid-token="}}}
    with pytest.raises(tools.InvalidToken):
        tools.decrypt_passwords(config_dict=config, submission_portals={"NCBI"}, key=key.decode())
    assert capsys.readouterr().err == ""

#*******************************************************************************
#                             encrypt_passwords
#*******************************************************************************

def test_encrypt_passwords__uses_supplied_key_and_saves_encrypted_gisaid_credentials(monkeypatch: pytest.MonkeyPatch) -> None:
    key = tools.Fernet.generate_key()
    config = {"GISAID": {"Username": "user", "Password": None, "Client-Id": None}}
    get_config_mock = Mock(side_effect=[config, config])
    save_yaml_mock = Mock()

    monkeypatch.setattr(tools, "get_config", get_config_mock)
    monkeypatch.setattr(tools, "determine_parent_database", Mock(return_value={"gisaid"}))
    monkeypatch.setattr(tools, "getpass", Mock(side_effect=["password123", "client123"]))
    monkeypatch.setattr(tools.file_handler, "save_yaml", save_yaml_mock)
    tools.encrypt_passwords(config_file="config.yaml", databases=["GISAID"], encryption_key=key.decode())
    encrypted_password = config["GISAID"]["Password"]
    encrypted_client_id = config["GISAID"]["Client-Id"]

    assert tools.Fernet(key).decrypt(encrypted_password) == b"password123"
    assert tools.Fernet(key).decrypt(encrypted_client_id) == b"client123"
    save_yaml_mock.assert_called_once_with(config_dict=config, yaml_path="config.yaml")
    assert get_config_mock.call_args_list == [
        call(config_file="config.yaml", databases=["GISAID"], passwords_validation=False),
        call(config_file="config.yaml", databases=["GISAID"], decrypt_key=key.decode())
    ]

def test_encrypt_passwords__generated_key_is_printed(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    config = {"GISAID": {"Password": None, "Client-Id": None}}
    monkeypatch.setattr(tools, "determine_parent_database", Mock(return_value={"gisaid"}))
    monkeypatch.setattr(tools, "get_config", Mock(side_effect=[config, config]))
    monkeypatch.setattr(tools, "getpass", Mock(side_effect=["password", "client-id"]))
    monkeypatch.setattr(tools.file_handler, "save_yaml", Mock())
    tools.encrypt_passwords(config_file="config.yaml", databases=["GISAID"], encryption_key=None)
    captured = capsys.readouterr()
    assert ("Save this key somewhere secure. It will be required for performing submission." in captured.out)
    assert "key: " in captured.out
    assert ("Credentials successfully loaded/encrypted into config file." in captured.out)

def test_encrypt_passwords__supplied_key_is_not_printed(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    key = tools.Fernet.generate_key()
    config = {"GISAID": {"Password": None, "Client-Id": None}}
    monkeypatch.setattr(tools, "determine_parent_database", Mock(return_value={"gisaid"}))
    monkeypatch.setattr(tools, "get_config", Mock(side_effect=[config, config]))
    monkeypatch.setattr(tools, "getpass", Mock(side_effect=["password", "client-id"]))
    monkeypatch.setattr(tools.file_handler, "save_yaml", Mock())
    tools.encrypt_passwords(config_file="config.yaml", databases=["GISAID"], encryption_key=key.decode())
    captured = capsys.readouterr()
    assert "Save this key somewhere secure." not in captured.out
    assert f"key: {key.decode()}" not in captured.out
    assert ("Credentials successfully loaded/encrypted into config file." in captured.out)

#*******************************************************************************
#                      get_submission_schema_config_name
#*******************************************************************************

@pytest.mark.parametrize(
    ("submission_portals", "expected"),
    [
        ({"ncbi"}, "ncbi_schema.py"),
        ({"gisaid"}, "gisaid_schema.py"),
        ({"ncbi", "gisaid"}, "ncbi_gisaid_schema.py"),
    ],
)
def test_get_submission_schema_config_name__maps_portals_to_schema(submission_portals, expected):
    assert tools.get_submission_schema_config_name(submission_portals) == expected

#*******************************************************************************
#                         get_submission_position
#*******************************************************************************

@pytest.mark.parametrize(
    ("database", "expected"),
    [
        ("BIOSAMPLE", 1),
        ("SRA", 1),
        ("GENBANK", 1),
        ("GISAID", 2),
    ],
)
def test_get_submission_position__reads_nested_submission_config(base_config: dict[str, Any], database: str, expected: int) -> None:
    assert tools.get_submission_position(base_config, database) == expected

def test_get_submission_position__accepts_already_nested_parent_config() -> None:
    assert tools.get_submission_position({"NCBI": {"Submission_Position": 2}}, "GENBANK") == 2
    assert tools.get_submission_position({"GISAID": {"Submission_Position": 1}}, "GISAID") == 1


@pytest.mark.parametrize(
    "config_dict",
    [
        {"Submission": {"NCBI": {}}},
        {"Submission": {"NCBI": {"Submission_Position": "1"}}},
        {"Submission": {"GISAID": {"Submission_Position": None}}},
    ],
)
def test_get_submission_position__returns_none_when_position_missing_or_non_int(config_dict: dict[str, Any]) -> None:
    assert tools.get_submission_position(config_dict, "GENBANK") is None


def test_get_submission_position__invalid_database_exits(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc:
        tools.get_submission_position({"Submission": {}}, "BADDB")
    assert exc.value.code == 1
    assert "Error: database BADDB is not a valid selection.\n" == capsys.readouterr().err

@pytest.mark.parametrize(
    ("config_dict", "expected_message"),
    [
        (
            {"Submission": {"NCBI": {"Submission_Position": 1}, "GISAID": {"Submission_Position": 1}}},
            "Error: Config file is incorrect. Submission position for GISAID '1' and GenBank '1' must both be either left empty, or set to '1' and '2' based on submission preference.\n",
        ),
        (
            {"Submission": {"NCBI": {"Submission_Position": 1}, "GISAID": {}}},
            "Error: Config file is incorrect. Submission position for GISAID 'None' and GenBank '1' must both be either left empty, or set to '1' and '2' based on submission preference.\n",
        ),
        (
            {"Submission": {"NCBI": {}, "GISAID": {"Submission_Position": 2}}},
            "Error: Config file is incorrect. Submission position for GISAID '2' and GenBank 'None' must both be either left empty, or set to '1' and '2' based on submission preference.\n",
        ),
    ],
)
def test_validate_submission_position__exits_for_inconsistent_or_duplicate_positions(config_dict: dict[str, Any], expected_message: str, capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc:
        tools.validate_submission_position(config_dict)

    assert exc.value.code == 1
    assert capsys.readouterr().err == expected_message


@pytest.mark.parametrize(
    "config_dict",
    [
        {"Submission": {"NCBI": {}, "GISAID": {}}},
        {"Submission": {"NCBI": {"Submission_Position": 1}, "GISAID": {"Submission_Position": 2}}},
        {"Submission": {"NCBI": {"Submission_Position": 2}, "GISAID": {"Submission_Position": 1}}},
    ],
)
def test_validate_submission_position__accepts_empty_or_distinct_positions(config_dict: dict[str, Any]) -> None:
    tools.validate_submission_position(config_dict)

#*******************************************************************************
#                         get_submission_schema
#*******************************************************************************

@pytest.mark.parametrize(
    ("portals", "expected"),
    [
        ({"ncbi"}, "ncbi_schema.py"),
        ({"gisaid"}, "gisaid_schema.py"),
        ({"ncbi", "gisaid"}, "ncbi_gisaid_schema.py"),
        (set(), "schema.py"),
    ],
)
def test_get_submission_schema_config_name__builds_expected_name(portals: set[str], expected: str) -> None:
    assert tools.get_submission_schema_config_name(portals) == expected

#*******************************************************************************
#                         get_submission_type
#*******************************************************************************

@pytest.mark.parametrize(
    ("test_flag", "expected"),
    [(True, "TEST"), (False, "PRODUCTION")],
)
def test_get_submission_type__maps_boolean_to_submission_type(test_flag: bool, expected: str) -> None:
    assert tools.get_submission_type(test_flag) == expected

#*******************************************************************************
#                           check_credentials
#*******************************************************************************

@pytest.mark.parametrize(
    ("config_dict", "expected_message"),
    [
        (
            {"Password": "p"},
            "Error: there is no Submission > NCBI > Username information in config file.\n",
        ),
        (
            {"Username": "", "Password": "p"},
            "Error: Submission > NCBI > Username in the config file cannot be empty.\n",
        ),
        (
            {"Username": None, "Password": "p"},
            "Error: Submission > NCBI > Username in the config file cannot be empty.\n",
        ),
        (
            {"Username": "u"},
            "Error: there is no Submission > NCBI > Password information in config file.\n",
        ),
        (
            {"Username": "u", "Password": ""},
            "Error: Submission > NCBI > Password in the config file cannot be empty.\n",
        ),
        (
            {"Username": "u", "Password": None},
            "Error: Submission > NCBI > Password in the config file cannot be empty.\n",
        ),
    ],
)
def test_check_credentials__ncbi_missing_or_empty_fields_exit(config_dict: dict[str, Any], expected_message: str, capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc:
        tools.check_credentials(config_dict, "NCBI")

    assert exc.value.code == 1
    assert expected_message == capsys.readouterr().err

@pytest.mark.parametrize(
    ("config_dict", "expected_message"),
    [
        (
            {"Username": "u", "Password": "p"},
            "Error: there is no Submission > GISAID > Client-Id information in config file.\n",
        ),
        (
            {"Username": "u", "Password": "p", "Client-Id": ""},
            "Error: Submission > GISAID > Client-Id in the config file cannot be empty.\n",
        ),
        (
            {"Username": "u", "Password": "p", "Client-Id": None},
            "Error: Submission > GISAID > Client-Id in the config file cannot be empty.\n",
        ),
    ],
)
def test_check_credentials__gisaid_requires_client_id(
    config_dict: dict[str, Any], expected_message: str, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(SystemExit) as exc:
        tools.check_credentials(config_dict, "GISAID")

    assert exc.value.code == 1
    assert expected_message == capsys.readouterr().err

@pytest.mark.parametrize(
    ("database", "config_dict"),
    [
        ("NCBI", {"Username": "u", "Password": "p"}),
        ("GISAID", {"Username": "u", "Password": "p", "Client-Id": "cid"}),
    ],
)
def test_check_credentials__valid_configs_do_not_exit(database: str, config_dict: dict[str, Any]) -> None:
    tools.check_credentials(config_dict, database)

@pytest.mark.parametrize(
    ("config_dict", "expected_message"),
    [
        (
            {"Username": "u"},
            "Error: there is no Submission > NCBI > Password information in config file.\n",
        ),
        (
            {"Username": "u", "Password": ""},
            "Error: Submission > NCBI > Password in the config file cannot be empty.\n",
        ),
        (
            {"Username": "u", "Password": None},
            "Error: Submission > NCBI > Password in the config file cannot be empty.\n",
        ),
    ],
)
def test_check_credentials__password_key_and_value_are_required_exactly(config_dict: dict[str, Any], expected_message: str, capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc:
        tools.check_credentials(config_dict, "NCBI")
    assert exc.value.code == 1
    assert capsys.readouterr().err == expected_message

#*******************************************************************************
#                           process_fasta_samples
#*******************************************************************************

def test_process_fasta_samples__returns_successful_one_to_one_merge(monkeypatch: pytest.MonkeyPatch) -> None:
    metadata = pd.DataFrame([{"sequence_name": "seq1", "other": "value"}])
    fasta_df = pd.DataFrame(
        [{"fasta_name_orig": "seq1", "fasta_sequence_orig": "ACTG", "fasta_description_orig": "seq1 desc"}]
    )
    monkeypatch.setattr(tools.file_handler, "load_fasta_file", Mock(return_value=fasta_df))

    merged = tools.process_fasta_samples(metadata, "sequence.fasta")

    assert merged.loc[0, "sequence_name"] == "seq1"
    assert merged.loc[0, "fasta_sequence_orig"] == "ACTG"


def test_process_fasta_samples__duplicate_fasta_names_exit(monkeypatch: pytest.MonkeyPatch, capsys) -> None:
    metadata = pd.DataFrame([{"sequence_name": "seq1"}])
    fasta_df = pd.DataFrame(
        [
            {"fasta_name_orig": "seq1", "fasta_sequence_orig": "ACTG", "fasta_description_orig": "seq1"},
            {"fasta_name_orig": "seq1", "fasta_sequence_orig": "TGCA", "fasta_description_orig": "seq1 duplicate"},
        ]
    )
    monkeypatch.setattr(tools.file_handler, "load_fasta_file", Mock(return_value=fasta_df))

    with pytest.raises(SystemExit) as exc:
        tools.process_fasta_samples(metadata, "sequence.fasta")

    assert exc.value.code == 1
    assert capsys.readouterr().err == (
        "Error: Sequences in fasta file must be unique at: sequence.fasta\n"
        "Duplicate Sequences\n"
        "ACTG\n"
        "TGCA\n"
    )

def test_process_fasta_samples__duplicate_metadata_merge_failure_exits(monkeypatch: pytest.MonkeyPatch, capsys) -> None:
    metadata = pd.DataFrame([{"sequence_name": "seq1"}, {"sequence_name": "seq1"}])
    fasta_df = pd.DataFrame(
        [{"fasta_name_orig": "seq1", "fasta_sequence_orig": "ACTG", "fasta_description_orig": "seq1"}]
    )
    monkeypatch.setattr(tools.file_handler, "load_fasta_file", Mock(return_value=fasta_df))

    with pytest.raises(SystemExit) as exc:
        tools.process_fasta_samples(metadata, "sequence.fasta")

    assert exc.value.code == 1
    assert "Error: Unable to merge fasta file to metadata file. Please validate there are not duplicate sequences in both files.\n" == capsys.readouterr().err


def test_process_fasta_samples__fasta_sequence_missing_from_metadata_exits(monkeypatch: pytest.MonkeyPatch, capsys) -> None:
    metadata = pd.DataFrame([{"sequence_name": "seq1"}])
    fasta_df = pd.DataFrame(
        [
            {"fasta_name_orig": "seq1", "fasta_sequence_orig": "ACTG", "fasta_description_orig": "seq1"},
            {"fasta_name_orig": "seq2", "fasta_sequence_orig": "TGCA", "fasta_description_orig": "seq2"},
        ]
    )
    monkeypatch.setattr(tools.file_handler, "load_fasta_file", Mock(return_value=fasta_df))

    with pytest.raises(SystemExit) as exc:
        tools.process_fasta_samples(metadata, "sequence.fasta")

    assert exc.value.code == 1
    assert capsys.readouterr().err == (
        "Error: Sequences in fasta file do not have an associated sequence in metadata file. Please update sequences below:\n"
        "1    seq2\n"
    )


def test_process_fasta_samples__metadata_sequence_missing_from_fasta_exits(monkeypatch: pytest.MonkeyPatch, capsys) -> None:
    metadata = pd.DataFrame([{"sequence_name": "seq1"}, {"sequence_name": "seq2"}])
    fasta_df = pd.DataFrame(
        [{"fasta_name_orig": "seq1", "fasta_sequence_orig": "ACTG", "fasta_description_orig": "seq1"}]
    )
    monkeypatch.setattr(tools.file_handler, "load_fasta_file", Mock(return_value=fasta_df))

    with pytest.raises(SystemExit) as exc:
        tools.process_fasta_samples(metadata, "sequence.fasta")

    assert exc.value.code == 1
    assert capsys.readouterr().err == (
        "Error: Sequences in metadata file do not have an associated sequence in fasta file. Please update sequences below:\n"
        "1    seq2\n"
    )

#*******************************************************************************
#                           parse_hold_date
#*******************************************************************************

@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (None, None),
        ("", ""),
        ("   ", "   "),
        ("2099-01-01", "2099-01-01"),
    ],
)
def test_parse_hold_date__preserves_empty_or_future_absolute_dates(value: str | None, expected: str | None) -> None:
    config = {"Submission": {"NCBI": {"Specified_Release_Date": value}}}
    parsed = tools.parse_hold_date(config)
    print(parsed["Submission"]["NCBI"]["Specified_Release_Date"])
    assert parsed["Submission"]["NCBI"]["Specified_Release_Date"] == expected


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("1 days", "2026-05-30"),
        ("2 weeks", "2026-06-12"),
        ("3 months", "2026-08-29"),
    ],
)
def test_parse_hold_date__parses_relative_offsets(
    monkeypatch: pytest.MonkeyPatch, value: str, expected: str
) -> None:
    monkeypatch.setattr(tools.pd, "Timestamp", FixedTimestamp)
    config = {"Submission": {"NCBI": {"Specified_Release_Date": value}}}
    parsed = tools.parse_hold_date(config)
    assert parsed["Submission"]["NCBI"]["Specified_Release_Date"] == expected

def test_parse_hold_date__future_absolute_date_updates_exact_original_key() -> None:
    config = {"Submission": {"NCBI": {"Specified_Release_Date": "2099-01-01"}}}
    parsed = tools.parse_hold_date(config)
    assert parsed is config
    assert parsed["Submission"]["NCBI"] == {"Specified_Release_Date": "2099-01-01"}

def test_parse_hold_date__today_is_not_a_valid_hold_date(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    class FixedDateTime:
        @staticmethod
        def today():
            return datetime(2026, 5, 29)

        @staticmethod
        def now():
            return datetime(2026, 5, 29)

        @staticmethod
        def strptime(value: str, fmt: str):
            return datetime.strptime(value, fmt)

    monkeypatch.setattr(tools, "datetime", FixedDateTime)

    config = {"Submission": {"NCBI": {"Specified_Release_Date": "2026-05-29"}}}

    with pytest.raises(SystemExit) as exc:
        tools.parse_hold_date(config)

    assert exc.value.code == 1
    assert capsys.readouterr().err == "Error: Config file field 'Specified_Release_Date', with value '2026-05-29' is not a valid date. To be a valid date format, it must be formatted as 'YYYY-MM-DD', with zero padding and it must be later than today 2026-05-29.\n"

@pytest.mark.parametrize(
    ("value", "expected_kind"),
    [
        ("2020-01-01", "past_date"),
        ("2026/01/01", "unable_parse"),
        ("tomorrow", "unable_parse"),
        ("1 years", "unable_parse"),
        ("abc days", "unable_parse"),
        ("2026-99-99", "unable_parse"),
    ],
)
def test_parse_hold_date__invalid_values_exit(value: str, expected_kind: str, capsys: pytest.CaptureFixture[str]) -> None:
    config = {"Submission": {"NCBI": {"Specified_Release_Date": value}}}
    with pytest.raises(SystemExit) as exc:
        tools.parse_hold_date(config)

    today = tools.datetime.now().strftime("%Y-%m-%d")

    if expected_kind == "past_date":
        expected = (
            f"Error: Config file field 'Specified_Release_Date', with value '{value}' is not a valid date. "
            "To be a valid date format, it must be formatted as 'YYYY-MM-DD', with zero padding "
            f"and it must be later than today {today}.\n"
        )
    else:
        expected = (
            f"Error: Unable to parse config file field 'Specified_Release_Date', with value '{value}'. "
            "For field to be valid it must be 'None', formatted as 'YYYY-MM-DD', with zero padding "
            f"and it must be later than today {today}.\n"
        )

    assert exc.value.code == 1
    assert capsys.readouterr().err == expected

def test_parse_hold_date__preserves_config_when_ncbi_or_field_missing() -> None:
    no_ncbi: dict[str, Any] = {"Submission": {"GISAID": {}}}
    no_field: dict[str, Any] = {"Submission": {"NCBI": {}}}
    assert tools.parse_hold_date(no_ncbi) is no_ncbi
    assert tools.parse_hold_date(no_field) is no_field

#*******************************************************************************
#                  database_specific_config_schema_updates
#*******************************************************************************

def test_database_specific_config_schema_updates__biosample_requires_package(configurable_schema: dict[str, Any]) -> None:
    updated = tools.database_specific_config_schema_updates(configurable_schema, ["BIOSAMPLE"])
    package = updated["Submission"]["schema"]["NCBI"]["schema"]["BioSample_Package"]
    assert package == {"required": True, "nullable": False}

def test_database_specific_config_schema_updates__genbank_requires_publication_fields(configurable_schema: dict[str, Any]) -> None:
    updated = tools.database_specific_config_schema_updates(configurable_schema, ["GENBANK"])
    ncbi = updated["Submission"]["schema"]["NCBI"]["schema"]
    assert ncbi["Publication_Title"] == {"required": True, "nullable": False}
    assert ncbi["Publication_Status"] == {"required": True, "nullable": False}

def test_database_specific_config_schema_updates__biosample_and_genbank_update_only_exact_ncbi_schema_keys(configurable_schema: dict[str, Any]) -> None:
    updated = tools.database_specific_config_schema_updates(configurable_schema, ["BIOSAMPLE", "GENBANK"])
    assert updated == {
        "Submission": {
            "schema": {
                "NCBI": {
                    "schema": {
                        "BioSample_Package": {
                            "required": True,
                            "nullable": False,
                        },
                        "Publication_Title": {
                            "required": True,
                            "nullable": False,
                        },
                        "Publication_Status": {
                            "required": True,
                            "nullable": False,
                        },
                    }
                }
            }
        }
    }
    assert set(updated["Submission"]["schema"]) == {"NCBI"}
    assert set(updated["Submission"]["schema"]["NCBI"]) == {"schema"}
    assert set(updated["Submission"]["schema"]["NCBI"]["schema"]) == {"BioSample_Package", "Publication_Title", "Publication_Status"}

#*******************************************************************************
#                        warn_deprecated_columns
#*******************************************************************************

def test_warn_deprecated_columns__does_not_exit_for_non_genbank() -> None:
    metadata = pd.DataFrame([{"src-Authority": "deprecated"}])
    tools.warn_deprecated_columns(["SRA"], metadata)

def test_warn_deprecated_columns__exits_for_genbank_deprecated_column(capsys: pytest.CaptureFixture[str]) -> None:
    metadata = pd.DataFrame([{"src-Authority": "deprecated"}])
    with pytest.raises(SystemExit) as exc:
        tools.warn_deprecated_columns(["GENBANK"], metadata)

    assert exc.value.code == 1
    assert "Error: GenBank columns '['src-Authority']' are deprecated and are no longer supported by GenBank. Please remove them before submission.\n" == capsys.readouterr().err

#*******************************************************************************
#                        pretty_print_pandera_errors
#*******************************************************************************

def _fake_schema_error(rows: list[dict[str, Any]]) -> SimpleNamespace:
    return SimpleNamespace(failure_cases=pd.DataFrame(rows))

def test_pretty_print_pandera_errors__prints_common_error_shapes(capsys: pytest.CaptureFixture[str]) -> None:
    error = _fake_schema_error(
        [
            {
                "schema_context": "DataFrame",
                "column": None,
                "check": "column_in_dataframe",
                "failure_case": "required_col",
                "index": 0,
            },
            {
                "schema_context": "Column",
                "column": "choice",
                "check": "isin(['A', 'B'])",
                "failure_case": "C",
                "index": 1,
            },
            {
                "schema_context": "Column",
                "column": "required",
                "check": "str_matches('^(?!\\s*$).+')",
                "failure_case": "",
                "index": 2,
            },
            {
                "schema_context": "DataFrame",
                "column": None,
                "check": '(lambda df: ~(df["a"].isnull() & df["b"].isnull()), ignore_na = False)',
                "failure_case": False,
                "index": 0,
            },
            {
                "schema_context": "Column",
                "column": "length",
                "check": "str_length(2, 10)",
                "failure_case": "x",
                "index": 4,
            },
            {
                "schema_context": "Column",
                "column": "status",
                "check": "str_matches('^(PENDING|SUBMITTED|\\WSUB\\d*\\W)$')",
                "failure_case": "bad",
                "index": 5,
            },
            {
                "schema_context": "Column",
                "column": "date",
                "check": "invalid_date_format",
                "failure_case": "01/01/2026",
                "index": 6,
            },
            {
                "schema_context": "Column",
                "column": "bs-title",
                "check": "same_value",
                "failure_case": "two values",
                "index": 7,
            },
            {
                "schema_context": "DataFrame",
                "column": None,
                "check": "column_ordered",
                "failure_case": "col_x",
                "index": 0,
            },
            {
                "schema_context": "DataFrame",
                "column": None,
                "check": "no_regex_column_match('sra-file_[2-9]\\d*')",
                "failure_case": None,
                "index": 0,
            },
        ]
    )

    tools.pretty_print_pandera_errors("metadata.csv", [error])

    captured = capsys.readouterr()
    assert captured.out == (
        "Error: file metadata.csv has the following error('s):\n"
        "(Note: Index position is calculated excluding column headers and the first row index value starting at '1'.)\n"
        "\n"
    )
    assert captured.err == (
        "Error: Missing required column 'required_col', ensure the file has not been modified and retry.\n"
        "\n"
        "Error: Column 'choice' has an incorrect value at index '2'. This field can only contain the values '['A', 'B']', you provided 'C'.\n"
        "\n"
        "Error: Column 'required' has an empty field at index '3' that is required. This field cannot be left blank.\n"
        "\n"
        "Error: In column group, every sample must have at least one non-null value in at least one of the following columns: '['a', 'b']'.\n"
        "\n"
        "Error: Column 'length' has a character limit of minimum '2', maximum '10'. The value 'x' at index '5' does not meet these requirements.\n"
        "\n"
        "Error: Column 'status' at index '6' has a value 'bad'. This field must be one of the accepted values: '['PENDING', 'SUBMITTED', 'SUB<numeric_values>']'.\n"
        "\n"
        "Error: Column 'date' at index '7' has a value '01/01/2026'. This field must be a valid date format based on ISO 8601: '[\"YYYY-MM-DD\", \"YYYY-MM\", or \"YYYY\"]'.\n"
        "\n"
        "Error: Column 'bs-title' must have the same value for every row as it is only used once and applies to the entire submission. This field is an internal NCBI field for the NCBI submission portal website (https://submit.ncbi.nlm.nih.gov/subs/) to aid you in identifying your submissions.\n"
        "\n"
        "Error: Column 'col_x' is incorrectly ordered for file 'metadata.csv'.\n"
        "\n"
        "Error: Column 'sra-file_#' is required, where # is the numeric value of the file for the SRA sample. (i.e. sra-file_1)\n"
        "\n"
    )

def test_pretty_print_pandera_errors__case_insensitive_accepted_values_are_split_exactly(capsys: pytest.CaptureFixture[str]) -> None:
    error = _fake_schema_error(
        [
            {
                "schema_context": "Column",
                "column": "database",
                "check": "str_matches('(?i)(\\\\W|^)(BIOSAMPLE|SRA|GENBANK|GISAID)(\\\\W|$)')",
                "failure_case": "bad",
                "index": 2,
            }
        ]
    )
    tools.pretty_print_pandera_errors("metadata.csv", [error])
    captured = capsys.readouterr()
    assert captured.out == "Error: file metadata.csv has the following error('s):\n(Note: Index position is calculated excluding column headers and the first row index value starting at '1'.)\n\n"
    assert captured.err == "Error: Column 'database' at index '3' has the value 'bad'. This field must be one of the accepted values: ['BIOSAMPLE', 'SRA', 'GENBANK', 'GISAID'].\n\n"

@pytest.mark.parametrize("column", ["bs-title", "bs-comment", "sra-title", "sra-comment", "gb-title", "gb-comment"])
def test_pretty_print_pandera_errors__same_value_columns_are_all_recognized(column: str, capsys: pytest.CaptureFixture[str]) -> None:
    error = _fake_schema_error(
        [
            {
                "schema_context": "Column",
                "column": column,
                "check": "same_value",
                "failure_case": "two values",
                "index": 0,
            }
        ]
    )

    tools.pretty_print_pandera_errors("metadata.csv", [error])
    captured = capsys.readouterr()
    assert captured.out == (
        "Error: file metadata.csv has the following error('s):\n"
        "(Note: Index position is calculated excluding column headers and the first row index value starting at '1'.)\n"
        "\n"
    )
    assert captured.err == (
        f"Error: Column '{column}' must have the same value for every row as it is only used once and applies to the entire submission. "
        "This field is an internal NCBI field for the NCBI submission portal website "
        "(https://submit.ncbi.nlm.nih.gov/subs/) to aid you in identifying your submissions.\n"
        "\n"
    )

def test_pretty_print_pandera_errors__groups_duplicate_errors(capsys: pytest.CaptureFixture[str]) -> None:
    error = _fake_schema_error(
        [
            {"schema_context": "Column", "column": "sample", "check": "field_uniqueness", "failure_case": "S1", "index": 0},
            {"schema_context": "Column", "column": "sample", "check": "field_uniqueness", "failure_case": "S1", "index": 3},
        ]
    )

    tools.pretty_print_pandera_errors("metadata.csv", [error])
    captured = capsys.readouterr()
    assert captured.out == (
        "Error: file metadata.csv has the following error('s):\n"
        "(Note: Index position is calculated excluding column headers and the first row index value starting at '1'.)\n"
        "\n"
    )
    assert captured.err == "Error: Column 'sample' with value 'S1' is duplicated at indices: '[1, 4]'.\n\n"


def test_pretty_print_pandera_errors__unknown_error_falls_back(capsys: pytest.CaptureFixture[str]) -> None:
    error = _fake_schema_error(
        [
            {
                "schema_context": "Column",
                "column": "x",
                "check": "custom_check",
                "failure_case": "bad",
                "index": 1,
            },
        ]
    )

    tools.pretty_print_pandera_errors("metadata.csv", [error])
    captured = capsys.readouterr()
    assert captured.out == (
        "Error: file metadata.csv has the following error('s):\n"
        "(Note: Index position is calculated excluding column headers and the first row index value starting at '1'.)\n"
        "\n"
        "Pandas(Index=0, schema_context='Column', column='x', check='custom_check', failure_case='bad', index=1)\n"
    )
    assert captured.err == (
        "Error: Unable to pretty print pandera error message for data. Error message is:\n"
        "Column: x\n"
        "Index: 2\n"
        "Value: bad\n"
        "Validator: custom_check\n"
        "If you would like to contribute to SeqSender. Make a issue on github reporting this error case to have a descriptive version of this error added.\n"
        "\n"
    )

#*******************************************************************************
#                             get_config
#*******************************************************************************

def _schema_file_text() -> str:
    return repr(
        {
            "Submission": {
                "schema": {
                    "NCBI": {
                        "schema": {
                            "BioSample_Package": {"required": False, "nullable": True},
                            "Publication_Title": {"required": False, "nullable": True},
                            "Publication_Status": {"required": False, "nullable": True},
                        }
                    }
                }
            }
        }
    )


@pytest.fixture
def patched_get_config_boundaries(monkeypatch: pytest.MonkeyPatch, base_config: dict[str, Any]) -> dict[str, Any]:
    FakeValidator.init_calls = []
    FakeValidator.calls = []
    FakeValidator.result = True
    FakeValidator.errors = {"field": ["bad"]}
    monkeypatch.setattr(tools, "Validator", FakeValidator)
    monkeypatch.setattr(tools.file_handler, "load_yaml", Mock(return_value=base_config))

    def fake_parse_hold_date(config_dict):
        return config_dict

    monkeypatch.setattr(tools, "parse_hold_date", Mock(side_effect=fake_parse_hold_date))
    monkeypatch.setattr(tools, "validate_submission_position", Mock())
    monkeypatch.setattr(tools, "password_encryption_config_schema_updates", Mock(side_effect=lambda schema, submission_portals: schema))
    monkeypatch.setattr(tools, "decrypt_passwords", Mock(side_effect=lambda config_dict, submission_portals, key: config_dict))

    open_calls: list[dict[str, Any]] = []

    class FakeSchemaFile:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def read(self):
            return _schema_file_text()

    def fake_open(file, mode="r", *args: Any, **kwargs: Any):
        open_calls.append(
            {
                "file": file,
                "mode": mode,
                "args": args,
                "kwargs": kwargs,
            }
        )
        return FakeSchemaFile()

    monkeypatch.setattr(builtins, "open", fake_open)
    return {
        "load_yaml": tools.file_handler.load_yaml,
        "parse_hold_date": tools.parse_hold_date,
        "validate_submission_position": tools.validate_submission_position,
        "password_encryption_config_schema_updates": tools.password_encryption_config_schema_updates,
        "decrypt_passwords": tools.decrypt_passwords,
        "open_calls": open_calls,
    }

@pytest.mark.parametrize("databases", [["GENBANK"], ["GISAID"], ["SRA"], ["BIOSAMPLE"]])
def test_get_config__does_not_validate_submission_position_unless_genbank_and_gisaid_both_selected(databases: list[str], patched_get_config_boundaries: dict[str, Any]) -> None:
    tools.get_config("config.yaml", databases)
    patched_get_config_boundaries["validate_submission_position"].assert_not_called()

def test_get_config__empty_database_list_exits(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc:
        tools.get_config("config.yaml", [])
    assert exc.value.code == 1
    assert "Error: Submission portals list cannot be empty.\n" == capsys.readouterr().err

@pytest.mark.parametrize(
    ("databases", "expected_schema_file"),
    [
        (["BIOSAMPLE"], "ncbi_schema.py"),
        (["SRA"], "ncbi_schema.py"),
        (["GENBANK"], "ncbi_schema.py"),
        (["GISAID"], "gisaid_schema.py"),
        (["GENBANK", "GISAID"], "ncbi_gisaid_schema.py"),
    ],
)
def test_get_config__selects_exact_schema_file_and_loads_yaml_with_exact_label(databases: list[str], expected_schema_file: str, patched_get_config_boundaries: dict[str, Any], base_config: dict[str, Any]) -> None:
    result = tools.get_config("config.yaml", databases)
    assert result is base_config["Submission"]
    patched_get_config_boundaries["load_yaml"].assert_called_once_with(yaml_type="Config file", yaml_path="config.yaml")
    assert patched_get_config_boundaries["open_calls"] == [
        {
            "file": os.path.join(
                tools.PROG_DIR,
                "config",
                "seqsender",
                "config_file",
                expected_schema_file,
            ),
            "mode": "r",
            "args": (),
            "kwargs": {},
        }
    ]
    expected_schema = eval(_schema_file_text())
    expected_schema = tools.database_specific_config_schema_updates(expected_schema, databases)
    assert FakeValidator.init_calls[0]["schema"] == expected_schema
    assert FakeValidator.calls == [(base_config, expected_schema)]

def test_get_config__opens_schema_file_with_explicit_read_mode(patched_get_config_boundaries: dict[str, Any]):
    tools.get_config("config.yaml", ["GENBANK", "GISAID"])
    assert patched_get_config_boundaries["open_calls"] == [
        {
            "file": os.path.join(
                tools.PROG_DIR,
                "config",
                "seqsender",
                "config_file",
                "ncbi_gisaid_schema.py",
            ),
            "mode": "r",
            "args": (),
            "kwargs": {},
        }
    ]

def test_get_config__password_validation_false_updates_schema_for_encryption(patched_get_config_boundaries: dict[str, Any]) -> None:
    tools.get_config("config.yaml", ["GISAID"], passwords_validation=False)
    patched_get_config_boundaries["password_encryption_config_schema_updates"].assert_called_once()
    call_args = patched_get_config_boundaries["password_encryption_config_schema_updates"].call_args
    assert call_args.args[1] == {"gisaid"}

def test_get_config__password_validation_true_does_not_update_encryption_schema(patched_get_config_boundaries: dict[str, Any]) -> None:
    tools.get_config("config.yaml", ["GISAID"], passwords_validation=True)
    patched_get_config_boundaries["password_encryption_config_schema_updates"].assert_not_called()

def test_get_config__decrypt_key_decrypts_credentials(patched_get_config_boundaries: dict[str, Any], base_config: dict[str, Any]) -> None:
    result = tools.get_config("config.yaml", ["GISAID"], decrypt_key="secret-key")
    patched_get_config_boundaries["decrypt_passwords"].assert_called_once_with(config_dict=base_config, submission_portals={"gisaid"}, key="secret-key")
    assert result is base_config["Submission"]

def test_get_config__without_decrypt_key_does_not_decrypt_credentials(patched_get_config_boundaries: dict[str, Any]) -> None:
    tools.get_config("config.yaml", ["GISAID"])
    patched_get_config_boundaries["decrypt_passwords"].assert_not_called()

def test_get_config__non_dict_yaml_exits(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    monkeypatch.setattr(tools.file_handler, "load_yaml", Mock(return_value=["not", "dict"]))
    with pytest.raises(SystemExit) as exc:
        tools.get_config("config.yaml", ["SRA"])

    assert exc.value.code == 1
    assert "Error: Config file is incorrect. File must be a valid yaml format.\n" == capsys.readouterr().err

def test_get_config__validator_failure_exits(patched_get_config_boundaries: dict[str, Any], capsys: pytest.CaptureFixture[str]) -> None:
    FakeValidator.result = False
    with pytest.raises(SystemExit) as exc:
        tools.get_config("config.yaml", ["SRA"])

    assert exc.value.code == 1
    assert capsys.readouterr().err == (
        "Error: Config file is not properly setup. Please correct config file based on issue below:\n"
        "{\n"
        '    "field": [\n'
        '        "bad"\n'
        "    ]\n"
        "}\n"
    )

def test_get_config__returns_submission_and_parses_hold_date(patched_get_config_boundaries: dict[str, Any], base_config: dict[str, Any]) -> None:
    result = tools.get_config("config.yaml", ["SRA"])
    assert result is base_config["Submission"]
    patched_get_config_boundaries["parse_hold_date"].assert_called_once_with(config_dict=base_config)
    patched_get_config_boundaries["validate_submission_position"].assert_not_called()

def test_get_config__validates_submission_position_when_genbank_and_gisaid_selected(patched_get_config_boundaries: dict[str, Any]) -> None:
    tools.get_config("config.yaml", ["GENBANK", "GISAID"])
    patched_get_config_boundaries["validate_submission_position"].assert_called_once()

#*******************************************************************************
#                            get_metadata
#*******************************************************************************

def test_get_metadata__validates_all_database_specific_schemas_with_expected_names_and_columns(monkeypatch: pytest.MonkeyPatch) -> None:
    metadata = pd.DataFrame(
        [
            {
                "bioproject": "PRJNA1",
                "organism": "Virus",
                "collection_date": "2026-01-01",
                "bs-sample_name": "BS1",
                "sra-sample_name": "SRA1",
                "sra-file_1": "reads.fastq.gz",
                "gb-sample_name": "GB1",
                "gs-sample_name": "GS1",
                "sequence_name": "seq1",
                "src-country": "USA",
            }
        ]
    )
    seqsender = DummySchema("seqsender")
    schemas: dict[str, DummySchema] = {}

    monkeypatch.setattr(tools.file_handler, "load_csv", Mock(return_value=metadata))
    monkeypatch.setattr(tools, "warn_deprecated_columns", Mock())
    monkeypatch.setattr(tools, "seqsender_schema", seqsender)

    def fake_import_module(name: str):
        schema = DummySchema(name)
        schemas[name] = schema
        return SimpleNamespace(schema=schema)

    monkeypatch.setattr(tools.importlib, "import_module", fake_import_module)

    result = tools.get_metadata(database=["BIOSAMPLE", "SRA", "GENBANK", "GISAID"],organism="FLU",metadata_file="metadata.csv",config_dict={"NCBI": {"BioSample_Package": "Pathogen.cl.1.0"}},skip_validation=False)

    assert result is metadata
    assert list(schemas) == [
        "config.biosample.Pathogen_cl_1_0",
        "config.sra.sra_schema",
        "config.genbank.genbank_schema",
        "config.genbank.genbank_flu_src_schema",
        "config.gisaid.gisaid_FLU_schema",
    ]

    assert len(seqsender.validate_calls) == 1
    assert seqsender.validate_calls[0][1] is True

    validated_by_schema = {
        name: schema.validate_calls[0][0].columns.tolist()
        for name, schema in schemas.items()
    }

    assert validated_by_schema["config.biosample.Pathogen_cl_1_0"] == [
        "bioproject",
        "organism",
        "collection_date",
        "bs-sample_name",
    ]
    assert validated_by_schema["config.sra.sra_schema"] == [
        "bioproject",
        "organism",
        "collection_date",
        "bs-sample_name",
        "sra-sample_name",
        "sra-file_1",
    ]
    assert validated_by_schema["config.genbank.genbank_schema"] == [
        "gb-sample_name",
        "sequence_name",
    ]
    assert validated_by_schema["config.genbank.genbank_flu_src_schema"] == [
        "bioproject",
        "organism",
        "collection_date",
        "gb-sample_name",
        "src-country",
    ]
    assert validated_by_schema["config.gisaid.gisaid_FLU_schema"] == [
        "collection_date",
        "gs-sample_name",
        "sequence_name",
    ]
    assert all(schema.validate_calls[0][1] is True for schema in schemas.values())

def test_get_metadata__schema_error_preserves_exact_database_schema_display_names(monkeypatch: pytest.MonkeyPatch) -> None:
    metadata = pd.DataFrame(
        [
            {
                "bioproject": "PRJNA1",
                "organism": "Virus",
                "collection_date": "2026-01-01",
                "bs-sample_name": "BS1",
                "sra-sample_name": "SRA1",
                "sra-file_1": "reads.fastq.gz",
                "gb-sample_name": "GB1",
                "sequence_name": "seq1",
                "src-country": "USA",
                "gs-sample_name": "GS1",
            }
        ]
    )
    schema_errors: dict[str, FakeSchemaErrors] = {}
    pretty_calls: list[dict[str, Any]] = []

    class NamedFailingSchema(DummySchema):
        def validate(self, df: pd.DataFrame, lazy: bool = True) -> pd.DataFrame:
            self.validate_calls.append((df.copy(), lazy))
            raise schema_errors[self.name]

    seqsender = DummySchema("seqsender")
    monkeypatch.setattr(tools.pandera.errors, "SchemaErrors", FakeSchemaErrors)
    monkeypatch.setattr(tools.file_handler, "load_csv", Mock(return_value=metadata))
    monkeypatch.setattr(tools, "warn_deprecated_columns", Mock())
    monkeypatch.setattr(tools, "seqsender_schema", seqsender)
    imported_schemas: dict[str, NamedFailingSchema] = {}

    def fake_import_module(name: str):
        schema = NamedFailingSchema(name)
        imported_schemas[name] = schema
        schema_errors[name] = FakeSchemaErrors(name)
        return SimpleNamespace(schema=schema)

    def fake_pretty_print_pandera_errors(**kwargs):
        pretty_calls.append(kwargs)

    monkeypatch.setattr(tools.importlib, "import_module", fake_import_module)
    monkeypatch.setattr(tools, "pretty_print_pandera_errors", fake_pretty_print_pandera_errors)

    with pytest.raises(SystemExit) as exc:
        tools.get_metadata(
            database=["BIOSAMPLE", "SRA", "GENBANK", "GISAID"],
            organism="FLU",
            metadata_file="metadata.csv",
            config_dict={"NCBI": {"BioSample_Package": "Pathogen.cl.1.0"}},
            skip_validation=False,
        )
    assert exc.value.code == 1
    assert list(imported_schemas) == [
        "config.biosample.Pathogen_cl_1_0",
        "config.sra.sra_schema",
        "config.genbank.genbank_schema",
        "config.genbank.genbank_flu_src_schema",
        "config.gisaid.gisaid_FLU_schema",
    ]
    assert pretty_calls == [
        {
            "file": "metadata.csv",
            "error_msgs": [
                schema_errors["config.biosample.Pathogen_cl_1_0"],
                schema_errors["config.sra.sra_schema"],
                schema_errors["config.genbank.genbank_schema"],
                schema_errors["config.genbank.genbank_flu_src_schema"],
                schema_errors["config.gisaid.gisaid_FLU_schema"],
            ],
        }
    ]

def test_get_metadata__default_skip_validation_is_false_and_validates_seqsender(monkeypatch: pytest.MonkeyPatch) -> None:
    metadata = pd.DataFrame(
        [
            {
                "sra-sample_name": "SRA1",
                "sra-file_1": "reads.fastq.gz",
                "bioproject": "PRJNA1",
                "organism": "Virus",
                "collection_date": "2026-01-01",
                "bs-sample_name": "BS1",
            }
        ]
    )
    seqsender = DummySchema("seqsender")
    sra_schema = DummySchema("sra")
    monkeypatch.setattr(tools.file_handler, "load_csv", Mock(return_value=metadata))
    monkeypatch.setattr(tools, "warn_deprecated_columns", Mock())
    monkeypatch.setattr(tools, "seqsender_schema", seqsender)

    def fake_import_module(name: str):
        return SimpleNamespace(schema=sra_schema)

    monkeypatch.setattr(tools.importlib, "import_module", fake_import_module)
    tools.get_metadata(database=["SRA"], organism="COV", metadata_file="metadata.csv", config_dict={"NCBI": {"BioSample_Package": "Pathogen.cl.1.0"}})
    assert len(seqsender.validate_calls) == 1
    assert seqsender.validate_calls[0][1] is True
    assert len(sra_schema.validate_calls) == 1
    assert sra_schema.validate_calls[0][1] is True

def test_get_metadata__non_flu_genbank_source_uses_general_source_schema_name(monkeypatch: pytest.MonkeyPatch) -> None:
    metadata = pd.DataFrame(
        [
            {
                "bioproject": "PRJNA1",
                "organism": "Virus",
                "collection_date": "2026-01-01",
                "gb-sample_name": "GB1",
                "sequence_name": "seq1",
                "src-country": "USA",
            }
        ]
    )
    imported: list[str] = []
    monkeypatch.setattr(tools.file_handler, "load_csv", Mock(return_value=metadata))
    monkeypatch.setattr(tools, "warn_deprecated_columns", Mock())
    monkeypatch.setattr(tools, "seqsender_schema", DummySchema("seqsender"))

    def fake_import_module(name: str):
        imported.append(name)
        return SimpleNamespace(schema=DummySchema(name))

    monkeypatch.setattr(tools.importlib, "import_module", fake_import_module)
    result = tools.get_metadata(database=["GENBANK"], organism="COV", metadata_file="metadata.csv", config_dict={"NCBI": {"BioSample_Package": "Pathogen.cl.1.0"}}, skip_validation=True)
    assert result is metadata
    assert imported == [
        "config.genbank.genbank_schema",
        "config.genbank.genbank_src_schema",
    ]

def test_get_metadata__skip_validation_loads_metadata_and_imports_needed_schemas(monkeypatch: pytest.MonkeyPatch) -> None:
    metadata = pd.DataFrame(
        [
            {
                "bioproject": "PRJNA1",
                "organism": "Virus",
                "collection_date": "2026-01-01",
                "bs-sample_name": "BS1",
                "sra-sample_name": "SRA1",
                "gb-sample_name": "GB1",
                "gs-sample_name": "GS1",
                "sequence_name": "seq1",
                "src-country": "USA",
            }
        ]
    )
    imported: list[str] = []

    monkeypatch.setattr(tools.file_handler, "load_csv", Mock(return_value=metadata))
    monkeypatch.setattr(tools, "warn_deprecated_columns", Mock())
    monkeypatch.setattr(tools, "seqsender_schema", DummySchema("seqsender"))

    def fake_import_module(name: str) -> types.SimpleNamespace:
        imported.append(name)
        return SimpleNamespace(schema=DummySchema(name))

    monkeypatch.setattr(tools.importlib, "import_module", fake_import_module)

    result = tools.get_metadata(
        database=["BIOSAMPLE", "SRA", "GENBANK", "GISAID"],
        organism="FLU",
        metadata_file="metadata.csv",
        config_dict={"NCBI": {"BioSample_Package": "Pathogen.cl.1.0"}},
        skip_validation=True,
    )

    assert result is metadata
    assert "config.biosample.Pathogen_cl_1_0" in imported
    assert "config.sra.sra_schema" in imported
    assert "config.genbank.genbank_schema" in imported
    assert "config.genbank.genbank_flu_src_schema" in imported
    assert "config.gisaid.gisaid_FLU_schema" in imported
    assert tools.seqsender_schema.update_columns_calls == [
        {
            "bioproject": {
                "checks": ("str_matches", (r"^(?!\s*$).+",), {}),
                "nullable": False,
                "required": True,
            }
        }
    ]

@pytest.mark.parametrize("database", [["BIOSAMPLE"], ["SRA"], ["BIOSAMPLE", "SRA"]])
def test_get_metadata__biosample_or_sra_requires_exact_bioproject_schema_update(monkeypatch: pytest.MonkeyPatch, database: list[str]) -> None:
    metadata = pd.DataFrame(
        [
            {
                "bioproject": "PRJNA1",
                "organism": "Virus",
                "collection_date": "2026-01-01",
                "bs-sample_name": "BS1",
                "sra-sample_name": "SRA1",
                "sra-file_1": "reads.fastq.gz",
            }
        ]
    )
    seqsender = DummySchema("seqsender")

    monkeypatch.setattr(tools.file_handler, "load_csv", Mock(return_value=metadata))
    monkeypatch.setattr(tools, "warn_deprecated_columns", Mock())
    monkeypatch.setattr(tools, "seqsender_schema", seqsender)

    def fake_import_module(name: str):
        return SimpleNamespace(schema=DummySchema(name))

    monkeypatch.setattr(tools.importlib, "import_module", fake_import_module)
    tools.get_metadata(database=database,organism="COV",metadata_file="metadata.csv",config_dict={"NCBI": {"BioSample_Package": "Pathogen.cl.1.0"}},skip_validation=True)

    assert len(seqsender.update_columns_calls) == 1
    payload = seqsender.update_columns_calls[0]
    assert set(payload) == {"bioproject"}
    assert set(payload["bioproject"]) == {"checks", "nullable", "required"}
    assert payload["bioproject"]["nullable"] is False
    assert payload["bioproject"]["required"] is True
    assert payload["bioproject"]["checks"][0] == "str_matches"
    assert payload["bioproject"]["checks"][1] == (r"^(?!\s*$).+",)
    assert payload["bioproject"]["checks"][2] == {}

@pytest.mark.parametrize("database", [["GENBANK"], ["GISAID"], ["GENBANK", "GISAID"]])
def test_get_metadata__does_not_require_bioproject_update_without_biosample_or_sra(monkeypatch: pytest.MonkeyPatch, database: list[str]) -> None:
    metadata = pd.DataFrame(
        [
            {
                "organism": "Virus",
                "collection_date": "2026-01-01",
                "gb-sample_name": "GB1",
                "gs-sample_name": "GS1",
                "sequence_name": "seq1",
            }
        ]
    )
    seqsender = DummySchema("seqsender")

    monkeypatch.setattr(tools.file_handler, "load_csv", Mock(return_value=metadata))
    monkeypatch.setattr(tools, "warn_deprecated_columns", Mock())
    monkeypatch.setattr(tools, "seqsender_schema", seqsender)

    def fake_import_module(name: str):
        return SimpleNamespace(schema=DummySchema(name))

    monkeypatch.setattr(tools.importlib, "import_module", fake_import_module)

    tools.get_metadata(database=database,organism="COV",metadata_file="metadata.csv",config_dict={"NCBI": {"BioSample_Package": "Pathogen.cl.1.0"}},skip_validation=True)

    assert seqsender.update_columns_calls == []

def test_get_metadata__validates_seqsender_and_database_specific_data(monkeypatch: pytest.MonkeyPatch) -> None:
    metadata = pd.DataFrame(
        [
            {
                "bioproject": "PRJNA1",
                "organism": "Virus",
                "collection_date": "2026-01-01",
                "sra-sample_name": "SRA1",
                "sra-file_1": "reads.fastq.gz",
                "bs-sample_name": "BS1",
            }
        ]
    )
    seqsender = DummySchema("seqsender")
    sra = DummySchema("sra")

    monkeypatch.setattr(tools.file_handler, "load_csv", Mock(return_value=metadata))
    monkeypatch.setattr(tools, "warn_deprecated_columns", Mock())
    monkeypatch.setattr(tools, "seqsender_schema", seqsender)
    def fake_import_module(name):
        return SimpleNamespace(schema=sra)

    monkeypatch.setattr(tools.importlib, "import_module", fake_import_module)

    result = tools.get_metadata(
        database=["SRA"],
        organism="COV",
        metadata_file="metadata.csv",
        config_dict={"NCBI": {"BioSample_Package": "Pathogen.cl.1.0"}},
        skip_validation=False,
    )

    assert result is metadata
    assert len(seqsender.validate_calls) == 1
    assert len(sra.validate_calls) == 1
    validated_sra_df = sra.validate_calls[0][0]
    assert set(validated_sra_df.columns) == {
        "bioproject",
        "organism",
        "collection_date",
        "sra-sample_name",
        "sra-file_1",
        "bs-sample_name",
    }

def test_get_metadata__default_validation_passes_lazy_true_to_seqsender_and_database_schema(monkeypatch: pytest.MonkeyPatch) -> None:
    metadata = pd.DataFrame(
        [
            {
                "bioproject": "PRJNA1",
                "organism": "Virus",
                "collection_date": "2026-01-01",
                "bs-sample_name": "BS1",
                "sra-sample_name": "SRA1",
                "sra-file_1": "reads.fastq.gz",
            }
        ]
    )
    seqsender = DummySchema("seqsender")
    sra_schema = DummySchema("sra")
    monkeypatch.setattr(tools.file_handler, "load_csv", Mock(return_value=metadata))
    monkeypatch.setattr(tools, "warn_deprecated_columns", Mock())
    monkeypatch.setattr(tools, "seqsender_schema", seqsender)

    def fake_import_module(name: str):
        assert name == "config.sra.sra_schema"
        return SimpleNamespace(schema=sra_schema)

    monkeypatch.setattr(tools.importlib, "import_module", fake_import_module)
    result = tools.get_metadata(database=["SRA"], organism="COV", metadata_file="metadata.csv", config_dict={"NCBI": {"BioSample_Package": "Pathogen.cl.1.0"}})
    assert result is metadata
    assert len(seqsender.validate_calls) == 1
    assert seqsender.validate_calls[0][1] is True
    assert len(sra_schema.validate_calls) == 1
    assert sra_schema.validate_calls[0][1] is True

def test_get_metadata__schema_error_pretty_prints_and_exits(monkeypatch: pytest.MonkeyPatch) -> None:
    metadata = pd.DataFrame([{"sra-sample_name": "SRA1"}])
    fake_error = FakeSchemaErrors("bad schema")
    seqsender = DummySchema("seqsender", raise_on_validate=fake_error)

    monkeypatch.setattr(tools.pandera.errors, "SchemaErrors", FakeSchemaErrors)
    monkeypatch.setattr(tools.file_handler, "load_csv", Mock(return_value=metadata))
    monkeypatch.setattr(tools, "warn_deprecated_columns", Mock())
    monkeypatch.setattr(tools, "seqsender_schema", seqsender)
    monkeypatch.setattr(tools, "pretty_print_pandera_errors", Mock())

    with pytest.raises(SystemExit) as exc:
        tools.get_metadata(
            database=[],
            organism="COV",
            metadata_file="metadata.csv",
            config_dict={},
            skip_validation=False,
        )

    assert exc.value.code == 1
    tools.pretty_print_pandera_errors.assert_called_once_with(file="metadata.csv", error_msgs=[fake_error])

#*******************************************************************************
#                            get_all_schema_files
#*******************************************************************************

def test_get_all_schema_files__walks_prog_config_and_returns_importable_schema_names(monkeypatch: pytest.MonkeyPatch) -> None:
    observed_walk_paths: list[str] = []

    def fake_walk(path: str):
        observed_walk_paths.append(path)
        return iter(
            [
                (
                    "/tmp/seqsender/config",
                    ["seqsender", "sra"],
                    ["root_schema.py", "not_python.txt"],
                ),
                (
                    "/tmp/seqsender/config/seqsender",
                    [],
                    ["seqsender_schema.py"],
                ),
                (
                    "/tmp/seqsender/config/sra",
                    [],
                    ["sra_schema.py"],
                ),
            ]
        )
    monkeypatch.setattr(tools, "PROG_DIR", "/tmp/seqsender")
    monkeypatch.setattr(tools.os, "walk", fake_walk)
    assert tools.get_all_schema_files() == ["config.root_schema", "config.seqsender.seqsender_schema", "config.sra.sra_schema"]
    assert observed_walk_paths == ["/tmp/seqsender/config"]

#*******************************************************************************
#                            load_schema
#*******************************************************************************

def test_load_schema__imports_schema_name_and_returns_schema_attribute(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    imported: list[str] = []
    schema = object()

    def fake_import_module(name: str):
        imported.append(name)
        return SimpleNamespace(schema=schema)

    monkeypatch.setattr(tools.importlib, "import_module", fake_import_module)

    assert tools.load_schema("config.sra.sra_schema") is schema
    assert imported == ["config.sra.sra_schema"]

#*******************************************************************************
#                            process_schema
#*******************************************************************************

def test_process_schema__returns_exact_template_dataframe_for_required_optional_group_and_sra_wildcard() -> None:
    schema = SimpleNamespace(
        columns={
            "required_col": FakeSchemaColumn(required=True, description="Required field"),
            "optional_col": FakeSchemaColumn(required=False, description=None),
            "group_col": FakeSchemaColumn(
                required=False,
                description='At least one required: Group: "host-data".',
            ),
            r"sra-file_[2-9]\d*": FakeSchemaColumn(
                required=False,
                description="Additional SRA files",
            ),
        }
    )

    result = tools.process_schema(schema)

    assert result.to_dict("records") == [
        {
            "column_name": "required_col",
            "required_column": "Required",
            "description": "Required field",
        },
        {
            "column_name": "optional_col",
            "required_column": "Optional",
            "description": None,
        },
        {
            "column_name": "group_col",
            "required_column": "At least one field required. Group: host-data",
            "description": 'At least one required: Group: "host-data".',
        },
        {
            "column_name": "sra-file_#",
            "required_column": "Optional",
            "description": "Additional SRA files",
        },
    ]

def test_process_schema__group_description_uses_last_group_marker() -> None:
    schema = SimpleNamespace(columns={"group_col": FakeSchemaColumn(required=False, description='Prefix Group: "wrong". At least one required: Group: "right".')})
    result = tools.process_schema(schema)
    assert result.to_dict("records") == [
        {
            "column_name": "group_col",
            "required_column": "At least one field required. Group: right",
            "description": 'Prefix Group: "wrong". At least one required: Group: "right".',
        }
    ]

#*******************************************************************************
#                          update_all_schema_templates
#*******************************************************************************
def test_update_all_schema_templates__continues_after_process_failure_and_writes_later_templates(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    calls: list[tuple[str, Any]] = []
    first_schema = SimpleNamespace(name="first")
    bad_schema = SimpleNamespace(name="bad")
    later_schema = SimpleNamespace(name="later")
    first_template = pd.DataFrame([{"column_name": "first", "required_column": "Required", "description": "desc"}])
    later_template = pd.DataFrame([{"column_name": "later", "required_column": "Optional", "description": "desc2"}])

    def fake_get_all_schema_files():
        return ["config.seqsender.upload_log_schema", "config.good_schema", "config.bad_load_schema", "config.bad_process_schema", "config.later_schema"]

    def fake_load_schema(schema_name: str):
        calls.append(("load_schema", schema_name))
        if schema_name == "config.bad_load_schema":
            raise RuntimeError("load failed")
        if schema_name == "config.good_schema":
            return first_schema
        if schema_name == "config.bad_process_schema":
            return bad_schema
        return later_schema

    def fake_process_schema(schema):
        calls.append(("process_schema", schema))
        if schema is bad_schema:
            raise RuntimeError("process failed")
        if schema is first_schema:
            return first_template
        return later_template

    to_csv_calls: list[dict[str, Any]] = []

    def fake_to_csv(self, path, *args, **kwargs):
        to_csv_calls.append({"records": self.to_dict("records"), "path": path, "args": args, "kwargs": kwargs})

    monkeypatch.setattr(tools, "get_all_schema_files", fake_get_all_schema_files)
    monkeypatch.setattr(tools, "load_schema", fake_load_schema)
    monkeypatch.setattr(tools, "process_schema", fake_process_schema)
    monkeypatch.setattr(pd.DataFrame, "to_csv", fake_to_csv)
    tools.update_all_schema_templates()
    assert calls == [
        ("load_schema", "config.good_schema"),
        ("process_schema", first_schema),
        ("load_schema", "config.bad_load_schema"),
        ("load_schema", "config.bad_process_schema"),
        ("process_schema", bad_schema),
        ("load_schema", "config.later_schema"),
        ("process_schema", later_schema),
    ]
    assert to_csv_calls == [
        {
            "records": [
                {
                    "column_name": "first",
                    "required_column": "Required",
                    "description": "desc",
                }
            ],
            "path": os.path.join(
                tools.PROG_DIR,
                "shiny",
                "templates",
                "config.good.schema_template.csv",
            ),
            "args": (),
            "kwargs": {"header": True, "index": False},
        },
        {
            "records": [
                {
                    "column_name": "later",
                    "required_column": "Optional",
                    "description": "desc2",
                }
            ],
            "path": os.path.join(
                tools.PROG_DIR,
                "shiny",
                "templates",
                "config.later.schema_template.csv",
            ),
            "args": (),
            "kwargs": {"header": True, "index": False},
        },
    ]
    assert capsys.readouterr().err == (
        'Warning: Unable to load schema "config.bad_load_schema".\n'
        "load failed\n"
        "Warning: Unable to process schema into metadata template.\n"
        "process failed\n"
    )
