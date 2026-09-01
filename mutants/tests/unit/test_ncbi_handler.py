from __future__ import annotations

import ftplib
import importlib.util
import os
import sys
import types
from pathlib import Path
from typing import Any

import pytest

def _source_root() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        # During mutmut, prefer the mutated source tree.
        mutant_src = parent / "src"
        if parent.name == "mutants" and (mutant_src / "ncbi_handler.py").exists():
            return mutant_src

        # Normal pytest run.
        normal_src = parent / "src"
        if (normal_src / "ncbi_handler.py").exists() and parent.name != "mutants":
            return normal_src

    raise RuntimeError("Could not find src/ncbi_handler.py")


SOURCE_DIR = _source_root()
MODULE_PATH = SOURCE_DIR / "ncbi_handler.py"

"""Unit tests initially generated with ChatGPT v5.5 "Deep Research" and "Thinking"
then modified for use and clarity with validation and coverage testing via mutmut."""

#*******************************************************************************
#                         create fake FTP server
#*******************************************************************************

class FakeFTP:
    def __init__(self, listings: list[list[str]] | None = None):
        self.listings = listings or [[]]
        self.nlst_calls = 0
        self.cwd_calls: list[str] = []
        self.mkd_calls: list[str] = []
        self.storbinary_calls: list[tuple[str, bytes]] = []
        self.storlines_calls: list[dict[str, Any]] = []
        self.retrbinary_calls: list[dict[str, Any]] = []
        self.login_calls: list[tuple[str, str]] = []
        self.storbinary_response = "226 Transfer complete"
        self.storlines_response = "226 Transfer complete"

    def nlst(self):
        idx = min(self.nlst_calls, len(self.listings) - 1)
        self.nlst_calls += 1
        return self.listings[idx]

    def cwd(self, folder: str):
        self.cwd_calls.append(folder)

    def mkd(self, folder: str):
        self.mkd_calls.append(folder)

    def storbinary(self, command: str, file_obj):
        self.storbinary_calls.append((command, file_obj.read()))
        return self.storbinary_response

    def storlines(self, command: str, file_obj):
        self.storlines_calls.append(
            {
                "command": command,
                "file_name": getattr(file_obj, "name", None),
                "file_mode": getattr(file_obj, "mode", None),
                "file_contents": file_obj.read(),
            }
        )
        return self.storlines_response

    def retrbinary(self, command: str, callback, blocksize: int):
        self.retrbinary_calls.append({"command": command, "callback": callback, "blocksize": blocksize})
        callback(b"<SubmissionStatus status='processed-ok' submission_id='SUB1'/>")

    def login(self, user: str, passwd: str):
        self.login_calls.append((user, passwd))

#*******************************************************************************
#                 create test ncbi_handler.py connections
#*******************************************************************************

@pytest.fixture()
def ncbi_handler_module(monkeypatch: pytest.MonkeyPatch):
    src_pkg: Any = types.ModuleType("src")
    src_pkg.__path__ = [str(SOURCE_DIR)]
    monkeypatch.setitem(sys.modules, "src", src_pkg)

    def alias_src_module(name, module):
        monkeypatch.setitem(sys.modules, f"src.{name}", module)
        monkeypatch.setitem(sys.modules, name, module)
        setattr(src_pkg, name, module)

    settings_stub: Any = types.ModuleType("settings")
    settings_stub.NCBI_FTP_HOST = "ftp.example.test"
    settings_stub.TABLE2ASN_EMAIL = "gb-admin@example.test"

    tools_stub: Any = types.ModuleType("tools")
    tools_stub.check_credentials_calls = []

    def check_credentials(**kwargs):
        tools_stub.check_credentials_calls.append(kwargs)

    tools_stub.check_credentials = check_credentials

    setup_stub: Any = types.ModuleType("setup")
    setup_stub.test_internet_connection_calls = []

    def test_internet_connection(databases):
        setup_stub.test_internet_connection_calls.append(databases)

    setup_stub.test_internet_connection = test_internet_connection

    xmltodict_stub: Any = types.ModuleType("xmltodict")

    def parse_xml(xml_bytes):
        import xml.etree.ElementTree as _ET
        root = _ET.fromstring(xml_bytes)

        def convert(elem):
            data = {f"@{key}": value for key, value in elem.attrib.items()}
            children_by_tag = {}
            for child in elem:
                children_by_tag.setdefault(child.tag, []).append(convert(child))
            for tag, values in children_by_tag.items():
                data[tag] = values if len(values) > 1 else values[0]
            text = (elem.text or "").strip()
            if text and not data:
                return text
            if text:
                data["#text"] = text
            return data

        return {root.tag: convert(root)}

    xmltodict_stub.parse = parse_xml
    for name in [
        "ncbi_handler",
        "src.ncbi_handler",
        "settings",
        "src.settings",
        "tools",
        "src.tools",
        "setup",
        "src.setup",
        "xmltodict",
    ]:
        sys.modules.pop(name, None)

    monkeypatch.setitem(sys.modules, "src", src_pkg)
    alias_src_module("settings", settings_stub)
    alias_src_module("tools", tools_stub)
    alias_src_module("setup", setup_stub)
    monkeypatch.setitem(sys.modules, "xmltodict", xmltodict_stub)

    sys.modules.pop("ncbi_handler", None)
    sys.modules.pop("src.ncbi_handler", None)
    spec = importlib.util.spec_from_file_location("ncbi_handler", MODULE_PATH)
    assert spec and spec.loader
    module: Any = importlib.util.module_from_spec(spec)
    sys.modules["ncbi_handler"] = module
    setattr(src_pkg, "ncbi_handler", module)
    spec.loader.exec_module(module)

    return module

#*******************************************************************************
#                     create fake ncbi config file
#*******************************************************************************

@pytest.fixture()
def config_dict() -> dict[str, Any]:
    return {
        "Username": "user",
        "Password": "pass",
        "Description": {
            "Organization": {
                "Submitter": {
                    "Email": "submitter@example.test",
                    "Alt_Email": "alt@example.test",
                }
            }
        },
    }

#*******************************************************************************
#                      standardize_submission_status
#*******************************************************************************

@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("SUBMITTED", "SUBMITTED"),
        (" Submitted ", "SUBMITTED"),
        ("CREATED", "CREATED"),
        ("Queued", "QUEUED"),
        ("Processing", "PROCESSING"),
        ("Failed", "FAILED"),
        ("Processed-OK", "PROCESSED"),
        ("Processed-Error", "ERROR"),
        ("Deleted", "DELETED"),
        ("Waiting", "WAITING"),
        ("Retried", "RETRIED"),
    ],
)
def test_standardize_submission_status__check_all_statuses(ncbi_handler_module, raw, expected):
    assert ncbi_handler_module.standardize_submission_status(raw) == expected

@pytest.mark.parametrize("raw", ["", "unknown", "not-ready", "not-made", "pending-review", "unprocessed-error"])
def test_standardize_submission_status__unknown_values_return_error(ncbi_handler_module, raw):
    assert ncbi_handler_module.standardize_submission_status(raw) == "ERROR"

#*******************************************************************************
#                           process_report_header
#*******************************************************************************

def test_process_report_header__reads_status_and_submission_id(monkeypatch, tmp_path, ncbi_handler_module):
    report = tmp_path / "report.xml"
    report.write_text('<SubmissionStatus status="processed-ok" submission_id="SUB123"/>')

    original_parse = ncbi_handler_module.xmltodict.parse
    observed_xml: list[bytes] = []

    def fake_parse(xml_bytes):
        observed_xml.append(xml_bytes)
        return original_parse(xml_bytes)

    monkeypatch.setattr(ncbi_handler_module.xmltodict, "parse", fake_parse)

    report_dict, status, submission_id = ncbi_handler_module.process_report_header(str(report))

    assert len(observed_xml) == 1
    assert isinstance(observed_xml[0], bytes)
    assert observed_xml[0].startswith(b"<?xml version='1.0' encoding='utf8'?>\n")
    assert b'<SubmissionStatus status="processed-ok" submission_id="SUB123" />' in observed_xml[0]
    assert report_dict["SubmissionStatus"]["@status"] == "processed-ok"
    assert status == "PROCESSED"
    assert submission_id == "SUB123"

def test_process_report_header__serializes_xml_with_exact_tostring_options(tmp_path, ncbi_handler_module, monkeypatch):
    report = tmp_path / "report.xml"
    report.write_text('<SubmissionStatus status="processed-ok" submission_id="SUB123" />', encoding="utf-8")

    original_tostring = ncbi_handler_module.ET.tostring
    observed_calls: list[dict[str, Any]] = []

    def fake_tostring(element, *args, **kwargs):
        observed_calls.append(
            {
                "element_tag": element.tag,
                "args": args,
                "kwargs": kwargs,
            }
        )
        return original_tostring(element, *args, **kwargs)

    monkeypatch.setattr(ncbi_handler_module.ET, "tostring", fake_tostring)

    _, status, submission_id = ncbi_handler_module.process_report_header(str(report))

    assert status == "PROCESSED"
    assert submission_id == "SUB123"
    assert observed_calls == [
        {
            "element_tag": "SubmissionStatus",
            "args": (),
            "kwargs": {
                "encoding": "utf8",
                "method": "xml",
            },
        }
    ]

# Validate default status and ID are assigned if missing file info
def test_process_report_header__defaults_missing_status_and_id(monkeypatch, tmp_path, ncbi_handler_module):
    report = tmp_path / "report.xml"
    report.write_text("<SubmissionStatus />")
    observed_statuses: List[str] = []

    def fake_standardize_submission_status(submission_status: str):
        observed_statuses.append(submission_status)
        return submission_status

    monkeypatch.setattr(ncbi_handler_module, "standardize_submission_status", fake_standardize_submission_status)
    _, status, submission_id = ncbi_handler_module.process_report_header(str(report))

    assert observed_statuses == ["SUBMITTED"]
    assert status == "SUBMITTED"
    assert submission_id == "PENDING"

#*******************************************************************************
#                                ncbi_login
#*******************************************************************************

def test_ncbi_login__success(monkeypatch, ncbi_handler_module, config_dict):
    ftp = FakeFTP()

    def fake_ftp_entry(host):
        return ftp

    monkeypatch.setattr(ncbi_handler_module.ftplib, "FTP", fake_ftp_entry)

    result = ncbi_handler_module.ncbi_login(config_dict)

    assert result is ftp
    assert ftp.login_calls == [("user", "pass")]


def test_ncbi_login__permission_error(monkeypatch, ncbi_handler_module, config_dict, capsys):
    def raise_perm(host):
        raise ftplib.error_perm("530 login incorrect")

    monkeypatch.setattr(ncbi_handler_module.ftplib, "FTP", raise_perm)

    with pytest.raises(UnboundLocalError):
        ncbi_handler_module.ncbi_login(config_dict)

    assert "Error: login error. Possible incorrect credentials for NCBI FTP site in config file. \nException 530 login incorrect\n" == capsys.readouterr().err

def test_ncbi_login__network_error_runs_network_test_and_exits(monkeypatch, ncbi_handler_module, config_dict, capsys):
    def raise_network(host):
        raise OSError("network down")

    monkeypatch.setattr(ncbi_handler_module.ftplib, "FTP", raise_network)

    with pytest.raises(SystemExit) as exc:
        ncbi_handler_module.ncbi_login(config_dict)

    assert exc.value.code == 1
    assert sys.modules["setup"].test_internet_connection_calls == [["NCBI"]]
    assert "Error unable to connect to FTP site. Running network test...\nException: network down\n" == capsys.readouterr().err

#*******************************************************************************
#                              ftp_upload_file
#*******************************************************************************

def test_ftp_upload_file__uses_basename_by_default(tmp_path, ncbi_handler_module):
    upload_file = tmp_path / "payload.txt"
    upload_file.write_text("hello")
    ftp = FakeFTP()

    assert ncbi_handler_module.ftp_upload_file(ftp, str(upload_file)) is ftp
    assert ftp.storbinary_calls == [("STOR payload.txt", b"hello")]

def test_ftp_upload_file__allows_custom_upload_name(tmp_path, ncbi_handler_module):
    upload_file = tmp_path / "payload.txt"
    upload_file.write_text("hello")
    ftp = FakeFTP()

    ncbi_handler_module.ftp_upload_file(ftp, str(upload_file), upload_name="remote.dat")

    assert ftp.storbinary_calls[0][0] == "STOR remote.dat"

def test_ftp_upload_file__exits_on_failed_transfer(tmp_path, ncbi_handler_module, capsys):
    upload_file = tmp_path / "payload.txt"
    upload_file.write_text("hello")
    ftp = FakeFTP()
    ftp.storbinary_response = "500 failed"

    with pytest.raises(SystemExit) as exc:
        ncbi_handler_module.ftp_upload_file(ftp, str(upload_file))

    assert exc.value.code == 1
    assert capsys.readouterr().err == f"Error: Uploading {upload_file} failed.\n"

#*******************************************************************************
#                           ftp_navigate_to_folder
#*******************************************************************************

def test_ftp_navigate_to_folder__direct_production_folder(ncbi_handler_module):
    ftp = FakeFTP(listings=[["Production"], ["submission_A"]])

    result = ncbi_handler_module.ftp_navigate_to_folder(
        ftp, "submission_A", "production", make_folder=False
    )

    assert result is ftp
    assert ftp.cwd_calls == ["Production", "submission_A"]

def test_ftp_navigate_to_folder__nested_submit_folder(ncbi_handler_module):
    ftp = FakeFTP(listings=[["submit"], ["submit"], ["submit"], ["Test"], ["submission_A"]])

    ncbi_handler_module.ftp_navigate_to_folder(ftp, "submission_A", "test")

    assert ftp.cwd_calls == ["submit", "Test", "submission_A"]

def test_ftp_navigate_to_folder__makes_submission_folder_when_requested(ncbi_handler_module):
    ftp = FakeFTP(listings=[["Production"], []])

    ncbi_handler_module.ftp_navigate_to_folder(
        ftp, "new_submission", "PRODUCTION", make_folder=True
    )

    assert ftp.mkd_calls == ["new_submission"]
    assert ftp.cwd_calls == ["Production", "new_submission"]

def test_ftp_navigate_to_folder__default_make_folder_is_false(ncbi_handler_module, capsys):
    ftp = FakeFTP(listings=[["Production"], []])

    with pytest.raises(SystemExit) as exc:
        ncbi_handler_module.ftp_navigate_to_folder(ftp, "missing_submission", "PRODUCTION")

    assert exc.value.code == 1
    assert ftp.mkd_calls == []
    assert ftp.cwd_calls == ["Production"]
    assert capsys.readouterr().err == "Error: Cannot find submission folder on NCBI FTP site.\n"

def test_ftp_navigate_to_folder__default_does_not_create_missing_folder(ncbi_handler_module, capsys):
    ftp = FakeFTP(listings=[["Production"], []])

    with pytest.raises(SystemExit) as exc:
        ncbi_handler_module.ftp_navigate_to_folder(ftp, "missing_submission", "PRODUCTION")

    assert exc.value.code == 1
    assert ftp.mkd_calls == []
    assert ftp.cwd_calls == ["Production"]
    assert capsys.readouterr().err == "Error: Cannot find submission folder on NCBI FTP site.\n"

def test_ftp_navigate_to_folder__exits_when_top_level_directory_missing(ncbi_handler_module, capsys):
    ftp = FakeFTP(listings=[["other"]])

    with pytest.raises(SystemExit) as exc:
        ncbi_handler_module.ftp_navigate_to_folder(ftp, "submission", "PRODUCTION")

    assert exc.value.code == 1
    assert "Error: Cannot find submission folder on NCBI FTP site.\n" == capsys.readouterr().err

def test_ftp_navigate_to_folder__exits_when_submit_exists_but_type_folder_missing(ncbi_handler_module, capsys):
    ftp = FakeFTP(listings=[["submit"], ["submit"], ["submit"], ["other"]])
    with pytest.raises(SystemExit) as exc:
        ncbi_handler_module.ftp_navigate_to_folder(ftp, "submission", "PRODUCTION")

    assert exc.value.code == 1
    assert ftp.cwd_calls == ["submit"]
    assert capsys.readouterr().err == "Error: Cannot find submission folder on NCBI FTP site.\n"

def test_ftp_navigate_to_folder__exits_when_submission_folder_missing(ncbi_handler_module, capsys):
    ftp = FakeFTP(listings=[["Production"], []])

    with pytest.raises(SystemExit) as exc:
        ncbi_handler_module.ftp_navigate_to_folder(ftp, "missing", "PRODUCTION")

    assert exc.value.code == 1
    assert "Error: Cannot find submission folder on NCBI FTP site.\n" == capsys.readouterr().err

def test_ftp_navigate_to_folder__does_not_create_existing_folder_even_when_make_folder_true(ncbi_handler_module):
    ftp = FakeFTP(listings=[["Production"], ["existing_submission"]])

    result = ncbi_handler_module.ftp_navigate_to_folder(ftp, "existing_submission", "PRODUCTION", make_folder=True)

    assert result is ftp
    assert ftp.mkd_calls == []
    assert ftp.cwd_calls == ["Production", "existing_submission"]

#*******************************************************************************
#                      create_submit_ready_file
#*******************************************************************************

def test_create_submit_ready_file__success(tmp_path, ncbi_handler_module):
    ftp = FakeFTP()

    assert ncbi_handler_module.create_submit_ready_file(ftp, str(tmp_path)) is ftp
    submit_ready_file = tmp_path / "submit.ready"
    assert submit_ready_file.exists()
    assert submit_ready_file.read_text() == ""
    assert ftp.storlines_calls == [
        {
            "command": "STOR submit.ready",
            "file_name": str(submit_ready_file),
            "file_mode": "rb",
            "file_contents": b"",
        }
    ]

def test_create_submit_ready_file__permission_denied_is_nonfatal(tmp_path, ncbi_handler_module, capsys):
    class PermissionDeniedFTP(FakeFTP):
        def storlines(self, command, file_obj):
            raise Exception("Error:550 submit.ready: Permission denied")

    ftp = PermissionDeniedFTP()

    result = ncbi_handler_module.create_submit_ready_file(ftp, str(tmp_path))

    assert result is ftp
    captured = capsys.readouterr()
    assert "" == captured.err
    assert "The submission has already been made and is currently processing.\n" == captured.out

def test_create_submit_ready_file__exits_on_other_error(tmp_path, ncbi_handler_module, capsys):
    class BrokenFTP(FakeFTP):
        def storlines(self, command, file_obj):
            raise RuntimeError("boom")

    with pytest.raises(SystemExit) as exc:
        ncbi_handler_module.create_submit_ready_file(BrokenFTP(), str(tmp_path))

    assert exc.value.code == 1
    assert capsys.readouterr().err == "Error: Unable to upload submit.ready file. boom\n"

def test_create_submit_ready_file__exits_on_bad_transfer_response(tmp_path, ncbi_handler_module, capsys):
    ftp = FakeFTP()
    ftp.storlines_response = "500 failed"

    with pytest.raises(SystemExit) as exc:
        ncbi_handler_module.create_submit_ready_file(ftp, str(tmp_path))

    assert exc.value.code == 1
    assert capsys.readouterr().err == "Error: submit.ready upload failed.\n"

#*******************************************************************************
#                             upload_raw_reads
#*******************************************************************************

def test_upload_raw_reads__uploads_existing_files_and_skips_blanks(tmp_path, monkeypatch, ncbi_handler_module):
    read1 = tmp_path / "r1.fastq.gz"
    read2 = tmp_path / "r2.fastq.gz"
    read1.write_text("r1")
    read2.write_text("r2")
    (tmp_path / "raw_reads_location.txt").write_text(f"{read1}\n\n{read2}\n")
    ftp = FakeFTP()
    calls: list[str] = []

    def fake_upload(ftp, upload_file, upload_name=None):
        calls.append(upload_file)
        return ftp

    monkeypatch.setattr(ncbi_handler_module, "ftp_upload_file", fake_upload)

    assert ncbi_handler_module.upload_raw_reads(ftp, str(tmp_path), "sub") is ftp
    assert calls == [str(read1), str(read2)]

def test_upload_raw_reads__opens_location_file_in_read_mode(tmp_path, ncbi_handler_module, monkeypatch):
    read1 = tmp_path / "r1.fastq.gz"
    read1.write_text("r1", encoding="utf-8")
    raw_reads_location = tmp_path / "raw_reads_location.txt"
    raw_reads_location.write_text(f"{read1}\n", encoding="utf-8")
    ftp = FakeFTP()
    original_open = open
    observed_open_calls = []

    def fake_open(file, mode="r", *args, **kwargs):
        observed_open_calls.append((str(file), mode, args, kwargs))
        return original_open(file, mode, *args, **kwargs)

    def fake_upload(**kwargs):
        return kwargs["ftp"]

    monkeypatch.setattr("builtins.open", fake_open)
    monkeypatch.setattr(ncbi_handler_module, "ftp_upload_file", fake_upload)
    assert ncbi_handler_module.upload_raw_reads(ftp, str(tmp_path), "sub") is ftp
    assert observed_open_calls == [(str(raw_reads_location), "r", (), {})]

def test_upload_raw_reads__exits_when_location_file_missing(tmp_path, ncbi_handler_module, capsys):
    with pytest.raises(SystemExit) as exc:
        ncbi_handler_module.upload_raw_reads(FakeFTP(), str(tmp_path), "sub")

    raw_reads_file = tmp_path / "raw_reads_location.txt"
    assert exc.value.code == 1
    assert  f"Error: Submission sub is missing raw reads file at {raw_reads_file}\n" == capsys.readouterr().err

def test_upload_raw_reads__exits_when_referenced_file_missing(tmp_path, ncbi_handler_module, capsys):
    missing = tmp_path / "missing.fastq.gz"
    (tmp_path / "raw_reads_location.txt").write_text(str(tmp_path / "missing.fastq.gz"))

    with pytest.raises(SystemExit) as exc:
        ncbi_handler_module.upload_raw_reads(FakeFTP(), str(tmp_path), "sub")

    assert exc.value.code == 1
    assert f"Error: Uploading files to SRA database failed. Possibly files have been moved or this is not a valid file: {missing}\n" == capsys.readouterr().err

#*******************************************************************************
#                            get_ncbi_report
#*******************************************************************************

def test_get_ncbi_report__downloads_report_when_present(tmp_path, monkeypatch, ncbi_handler_module, config_dict, capsys):
    ftp = FakeFTP(listings=[["Production"], ["sub_BIOSAMPLE"], ["report.xml"]])

    def fake_ncbi_entry(config_dict):
        return ftp

    monkeypatch.setattr(ncbi_handler_module, "ncbi_login", fake_ncbi_entry)

    report_path = ncbi_handler_module.get_ncbi_report("BIOSAMPLE", "sub", str(tmp_path), config_dict, "PRODUCTION")

    assert report_path == str(tmp_path / "report.xml")
    assert (tmp_path / "report.xml").read_bytes().startswith(b"<SubmissionStatus")
    assert len(ftp.retrbinary_calls) == 1
    assert ftp.retrbinary_calls[0]["command"] == "RETR report.xml"
    assert callable(ftp.retrbinary_calls[0]["callback"])
    assert ftp.retrbinary_calls[0]["blocksize"] == 262144
    assert "Downloading report.xml\n" == capsys.readouterr().out
    assert sys.modules["tools"].check_credentials_calls == [{"config_dict": config_dict, "database": "NCBI"}]

def test_get_ncbi_report__returns_none_when_report_absent(tmp_path, monkeypatch, ncbi_handler_module, config_dict, capsys):
    ftp = FakeFTP(listings=[["Production"], ["sub_BIOSAMPLE"], []])

    def fake_ncbi_entry(config_dict):
        return ftp

    monkeypatch.setattr(ncbi_handler_module, "ncbi_login", fake_ncbi_entry)

    assert ncbi_handler_module.get_ncbi_report("BIOSAMPLE", "sub", str(tmp_path), config_dict, "PRODUCTION") is None
    assert "The report.xml has not yet been generated.\n" == capsys.readouterr().out
    assert sys.modules["tools"].check_credentials_calls == [{"config_dict": config_dict, "database": "NCBI"}]

def test_get_ncbi_report__exits_on_ftp_error(tmp_path, monkeypatch, ncbi_handler_module, config_dict, capsys):
    def raise_ftp(config_dict):
        raise ftplib.error_temp("temporarily down")

    monkeypatch.setattr(ncbi_handler_module, "ncbi_login", raise_ftp)

    with pytest.raises(SystemExit) as exc:
        ncbi_handler_module.get_ncbi_report(
            "BIOSAMPLE", "sub", str(tmp_path), config_dict, "PRODUCTION"
        )

    assert exc.value.code == 1
    assert "\nError: temporarily down\n" == capsys.readouterr().err

#*******************************************************************************
#                              submit_ncbi
#*******************************************************************************

def test_submit_ncbi__sra_submission(tmp_path, monkeypatch, ncbi_handler_module, config_dict, capsys):
    (tmp_path / "submission.xml").write_text("xml")
    ftp = FakeFTP()
    events: list[tuple[str, Any]] = []

    def fake_sleep(seconds):
        events.append({"function": "sleep", "seconds": seconds})

    def fake_ncbi_login(config_dict):
        events.append({"function": "login", "config_dict": config_dict})
        return ftp

    def fake_ftp_navigate_to_folder(**kwargs):
        events.append({"function": "navigate", **kwargs})
        assert kwargs["ftp"] is ftp
        return kwargs["ftp"]

    def fake_ftp_upload_file(**kwargs):
        events.append({"function": "upload", **kwargs})
        assert kwargs["ftp"] is ftp
        return kwargs["ftp"]

    def fake_ftp_upload_raw_reads(**kwargs):
        events.append({"function": "raw", **kwargs})
        assert kwargs["ftp"] is ftp
        return kwargs["ftp"]

    def fake_create_submit_ready_file(**kwargs):
        events.append({"function": "ready", **kwargs})
        assert kwargs["ftp"] is ftp
        return kwargs["ftp"]

    monkeypatch.setattr(ncbi_handler_module.time, "sleep", fake_sleep)
    monkeypatch.setattr(ncbi_handler_module, "ncbi_login", fake_ncbi_login)
    monkeypatch.setattr(ncbi_handler_module, "ftp_navigate_to_folder", fake_ftp_navigate_to_folder)
    monkeypatch.setattr(ncbi_handler_module, "ftp_upload_file", fake_ftp_upload_file)
    monkeypatch.setattr(ncbi_handler_module, "upload_raw_reads", fake_ftp_upload_raw_reads)
    monkeypatch.setattr(ncbi_handler_module, "create_submit_ready_file", fake_create_submit_ready_file)

    ncbi_handler_module.submit_ncbi("SRA", "sub", str(tmp_path), config_dict, "TEST")
    assert sys.modules["tools"].check_credentials_calls == [{"config_dict": config_dict, "database": "NCBI"}]
    assert events == [
        {"function": "sleep", "seconds": 5},
        {"function": "login", "config_dict": config_dict},
        {
            "function": "navigate",
            "ftp": ftp,
            "folder_name": "sub_SRA",
            "submission_type": "TEST",
            "make_folder": True,
        },
        {
            "function": "upload",
            "ftp": ftp,
            "upload_file": str(tmp_path / "submission.xml"),
        },
        {
            "function": "raw",
            "ftp": ftp,
            "submission_dir": str(tmp_path),
            "submission_name": "sub",
        },
        {
            "function": "ready",
            "ftp": ftp,
            "submission_dir": str(tmp_path),
        },
    ]
    assert capsys.readouterr().out == (
        "Uploading sample files to NCBI-SRA, as a 'TEST' submission. If this is not intended, interrupt immediately.\n"
        "Connecting to NCBI FTP Server\n"
        "Submission name: sub_SRA\n"
        "Submitting 'sub'\n"
    )

def test_submit_ncbi__genbank_uploads_zip_not_raw_reads(tmp_path, monkeypatch, ncbi_handler_module, config_dict, capsys):
    (tmp_path / "submission.xml").write_text("xml")
    (tmp_path / "sub.zip").write_text("zip")
    ftp = FakeFTP()
    events: list[tuple[str, Any]] = []
    raw_called = False

    def fake_sleep(seconds):
        events.append({"function": "sleep", "seconds": seconds})

    def fake_ncbi_login(config_dict):
        events.append({"function": "login", "config_dict": config_dict})
        return ftp

    def fake_ftp_navigate_to_folder(**kwargs):
        events.append({"function": "navigate", **kwargs})
        assert kwargs["ftp"] is ftp
        return kwargs["ftp"]

    def fake_ftp_upload_file(**kwargs):
        events.append({"function": "upload", **kwargs})
        assert kwargs["ftp"] is ftp
        return kwargs["ftp"]

    def fake_create_submit_ready_file(**kwargs):
        events.append({"function": "ready", **kwargs})
        assert kwargs["ftp"] is ftp
        return kwargs["ftp"]

    def fake_raw(*args, **kwargs):
        nonlocal raw_called
        raw_called = True
        return ftp

    monkeypatch.setattr(ncbi_handler_module.time, "sleep", fake_sleep)
    monkeypatch.setattr(ncbi_handler_module, "ncbi_login", fake_ncbi_login)
    monkeypatch.setattr(ncbi_handler_module, "ftp_navigate_to_folder", fake_ftp_navigate_to_folder)
    monkeypatch.setattr(ncbi_handler_module, "ftp_upload_file", fake_ftp_upload_file)
    monkeypatch.setattr(ncbi_handler_module, "upload_raw_reads", fake_raw)
    monkeypatch.setattr(ncbi_handler_module, "create_submit_ready_file", fake_create_submit_ready_file)

    ncbi_handler_module.submit_ncbi("GENBANK", "sub", str(tmp_path), config_dict, "PRODUCTION")

    assert raw_called is False
    assert sys.modules["tools"].check_credentials_calls == [{"config_dict": config_dict, "database": "NCBI"}]
    assert events == [
        {"function": "sleep", "seconds": 5},
        {"function": "login", "config_dict": config_dict},
        {
            "function": "navigate",
            "ftp": ftp,
            "folder_name": "sub_GENBANK",
            "submission_type": "PRODUCTION",
            "make_folder": True,
        },
        {
            "function": "upload",
            "ftp": ftp,
            "upload_file": str(tmp_path / "submission.xml"),
        },
        {
            "function": "upload",
            "ftp": ftp,
            "upload_file": str(tmp_path / "sub.zip"),
        },
        {
            "function": "ready",
            "ftp": ftp,
            "submission_dir": str(tmp_path),
        },
    ]
    assert capsys.readouterr().out == (
        "Uploading sample files to NCBI-GENBANK, as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.\n"
        "Connecting to NCBI FTP Server\n"
        "Submission name: sub_GENBANK\n"
        "Submitting 'sub'\n"
    )

def test_submit_ncbi__exits_on_ftp_error(tmp_path, monkeypatch, ncbi_handler_module, config_dict, capsys):
    def fake_sleep(seconds):
        return None

    def fake_ncbi_login(config_dict):
        raise ftplib.error_temp("ftp down")

    monkeypatch.setattr(ncbi_handler_module.time, "sleep", fake_sleep)
    monkeypatch.setattr(ncbi_handler_module, "ncbi_login", fake_ncbi_login)

    with pytest.raises(SystemExit) as exc:
        ncbi_handler_module.submit_ncbi("SRA", "sub", str(tmp_path), config_dict, "TEST")

    assert exc.value.code == 1
    captured = capsys.readouterr()
    assert "Uploading sample files to NCBI-SRA, as a 'TEST' submission. If this is not intended, interrupt immediately.\n" == captured.out
    assert captured.err == "\nError: ftp down\n"

#*******************************************************************************
#                       create fake email server
#*******************************************************************************

class FakeSMTP:
    instances: list["FakeSMTP"] = []

    def __init__(self, host: str):
        self.host = host
        self.sent: list[tuple[str, list[str], str]] = []
        FakeSMTP.instances.append(self)

    def sendmail(self, from_email: str, to_email: list[str], message: str):
        self.sent.append((from_email, to_email, message))

#*******************************************************************************
#                             email_table2asn
#*******************************************************************************

def test_email_table2asn__test_submission_sends_to_submitter_only(tmp_path, monkeypatch, ncbi_handler_module, config_dict, capsys):
    (tmp_path / "sub.sqn").write_bytes(b"sqn")
    FakeSMTP.instances.clear()
    sleep_calls: list[int] = []

    def fake_sleep(seconds):
        return sleep_calls.append(seconds)

    monkeypatch.setattr(ncbi_handler_module.time, "sleep", fake_sleep)
    monkeypatch.setattr(ncbi_handler_module.smtplib, "SMTP", FakeSMTP)

    status = ncbi_handler_module.email_table2asn("sub", str(tmp_path), config_dict, "TEST")

    captured = capsys.readouterr()
    assert status == "PROCESSED"
    assert sleep_calls == [5]
    assert captured.out == "Emailing table2asn sqn file to submitter 'submitter@example.test' as a 'TEST' submission. If this is not intended, interrupt immediately.\n"
    assert captured.err == ""

    smtp = FakeSMTP.instances[0]
    assert smtp.host == "localhost"
    from_email, to_email, message = smtp.sent[0]
    assert from_email == "submitter@example.test"
    assert to_email == ["submitter@example.test"]
    assert "Content-Type: multipart/multipart;" in message
    assert "Subject: sub table2asn submission" in message
    assert "From: submitter@example.test" in message
    assert "To: submitter@example.test" in message
    assert "Cc: alt@example.test" in message
    assert "Content-Disposition: attachment; filename=sub.sqn" in message
    assert "Name=\"sub.sqn\"" in message

def test_email_table2asn__production_submission_sends_to_table2asn_and_ccs_submitter(tmp_path, monkeypatch, ncbi_handler_module, config_dict, capsys):
    (tmp_path / "sub.sqn").write_bytes(b"sqn")
    FakeSMTP.instances.clear()
    sleep_calls: list[int] = []

    def fake_sleep(seconds):
        return sleep_calls.append(seconds)

    monkeypatch.setattr(ncbi_handler_module.time, "sleep", fake_sleep)
    monkeypatch.setattr(ncbi_handler_module.smtplib, "SMTP", FakeSMTP)

    status = ncbi_handler_module.email_table2asn("sub", str(tmp_path), config_dict, "PRODUCTION")

    assert status == "PROCESSED"
    assert sleep_calls == [5]
    captured = capsys.readouterr()
    assert "Emailing table2asn sqn file to NCBI-GENBANK 'gb-admin@example.test', as a 'PRODUCTION' submission. If this is not intended, interrupt immediately.\n" == captured.out
    assert captured.err == ""

    smtp = FakeSMTP.instances[0]
    assert smtp.host == "localhost"
    from_email, to_email, message = FakeSMTP.instances[0].sent[0]
    assert "Content-Type: multipart/multipart;" in message
    assert "Subject: sub table2asn submission" in message
    assert "From: submitter@example.test" in message
    assert "To: gb-admin@example.test" in message
    assert "Cc: submitter@example.test, alt@example.test" in message
    assert "Content-Disposition: attachment; filename=sub.sqn" in message
    assert "Name=\"sub.sqn\"" in message

def test_email_table2asn__opens_sqn_file_in_binary_mode(tmp_path, monkeypatch, ncbi_handler_module, config_dict):
    sqn_file = tmp_path / "sub.sqn"
    sqn_file.write_bytes(b"sqn")

    original_open = open
    observed_open_calls: list[tuple[str, str]] = []

    def fake_open(file, mode="r", *args, **kwargs):
        observed_open_calls.append((str(file), mode))
        return original_open(file, mode, *args, **kwargs)

    def fake_sleep(seconds):
        return None

    monkeypatch.setattr("builtins.open", fake_open)
    monkeypatch.setattr(ncbi_handler_module.time, "sleep", fake_sleep)
    monkeypatch.setattr(ncbi_handler_module.smtplib, "SMTP", FakeSMTP)

    status = ncbi_handler_module.email_table2asn("sub", str(tmp_path), config_dict, "TEST")

    assert status == "PROCESSED"
    assert (str(sqn_file), "rb") in observed_open_calls

def test_email_table2asn__exits_for_invalid_submission_type(tmp_path, monkeypatch, ncbi_handler_module, config_dict, capsys):
    def fake_sleep(seconds):
        return None

    monkeypatch.setattr(ncbi_handler_module.time, "sleep", fake_sleep)

    with pytest.raises(SystemExit) as exc:
        ncbi_handler_module.email_table2asn("sub", str(tmp_path), config_dict, "BAD")

    assert exc.value.code == 1
    assert "Error: Submission type 'BAD' is not a valid option.\n" == capsys.readouterr().err

def test_email_table2asn__returns_error_on_send_failure(tmp_path, monkeypatch, ncbi_handler_module, config_dict, capsys):
    class BrokenSMTP:
        def __init__(self, host):
            pass
        def sendmail(self, from_email, to_email, message):
            raise RuntimeError("send failed")

    (tmp_path / "sub.sqn").write_bytes(b"sqn")

    def fake_sleep(seconds):
        return None

    monkeypatch.setattr(ncbi_handler_module.time, "sleep", fake_sleep)
    monkeypatch.setattr(ncbi_handler_module.smtplib, "SMTP", BrokenSMTP)

    assert ncbi_handler_module.email_table2asn("sub", str(tmp_path), config_dict, "TEST") == "ERROR"
    captured = capsys.readouterr()
    assert "Emailing table2asn sqn file to submitter 'submitter@example.test' as a 'TEST' submission. If this is not intended, interrupt immediately.\n" == captured.out
    assert captured.err == (
        "Error: Unable to send mail automatically. If unable to email, submission can be made manually using the sqn file.\n"
        f"sqn_file:{tmp_path / 'sub.sqn'}\n"
        "send failed\n"
    )
