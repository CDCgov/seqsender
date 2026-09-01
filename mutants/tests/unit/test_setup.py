from __future__ import annotations

import gzip
import importlib.util
import io
import os
import socket
import subprocess
import sys
import types
import zipfile
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from typing import Any

def _source_root() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        mutant_src = parent / "src"
        if parent.name == "mutants" and (mutant_src / "setup.py").exists():
            return mutant_src

        normal_src = parent / "src"
        if (normal_src / "setup.py").exists() and parent.name != "mutants":
            return normal_src

    raise RuntimeError("Could not find src/setup.py")


SOURCE_DIR = _source_root()
MODULE_PATH = SOURCE_DIR / "setup.py"
if str(SOURCE_DIR) not in sys.path:
    sys.path.insert(0, str(SOURCE_DIR))


"""Unit tests initially generated with ChatGPT v5.5 "Deep Research" and "Thinking"
then modified for use and clarity with validation and coverage testing via mutmut."""

#*******************************************************************************
#                      create test setup.py connections
#*******************************************************************************

@pytest.fixture()
def setup_module(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.syspath_prepend(str(SOURCE_DIR))
    src_pkg: Any = types.ModuleType("src")
    src_pkg.__path__ = [str(SOURCE_DIR)]
    monkeypatch.setitem(sys.modules, "src", src_pkg)

    def alias_src_module(name, module):
        monkeypatch.setitem(sys.modules, name, module)
        monkeypatch.setitem(sys.modules, f"src.{name}", module)
        setattr(src_pkg, name, module)

    tools_stub: Any = types.ModuleType("tools")

    def update_all_schema_templates():
        return None

    tools_stub.update_all_schema_templates = update_all_schema_templates
    alias_src_module("tools", tools_stub)
    settings_stub: Any = types.ModuleType("settings")
    settings_stub.NCBI_FTP_HOST = "ftp-private.ncbi.nlm.nih.gov"
    alias_src_module("settings", settings_stub)

    xmltodict_stub: Any = types.ModuleType("xmltodict")

    def _parse_xmltodict(xml_bytes):
        import xml.etree.ElementTree as ET

        root = ET.fromstring(xml_bytes)
        attributes = []
        for attr in root.findall(".//Attribute"):
            item = {f"@{key}": value for key, value in attr.attrib.items()}
            for child in attr:
                item[child.tag] = child.text or ""
                if child.attrib:
                    item[child.tag] = {f"@{k}": v for k, v in child.attrib.items()} | {
                        grand.tag: grand.text or "" for grand in child
                    }
            attributes.append(item)
        return {"BioSamplePackages": {"Package": {"Attribute": attributes}}}

    xmltodict_stub.parse = _parse_xmltodict
    monkeypatch.setitem(sys.modules, "xmltodict", xmltodict_stub)

    sys.modules.pop("setup", None)
    sys.modules.pop("src.setup", None)
    spec = importlib.util.spec_from_file_location("setup", MODULE_PATH)
    assert spec and spec.loader
    module: Any = importlib.util.module_from_spec(spec)
    sys.modules["setup"] = module
    setattr(src_pkg, "setup", module)
    spec.loader.exec_module(module)
    return module

def _write_text(path: Path, text: str = "x") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")

def _make_test_data_tree(root: Path, organism: str = "FLU") -> None:
    base = root / "test_data" / organism
    lower = organism.lower()
    base.mkdir(parents=True)
    _write_text(base / f"{lower}_config.yaml", "Submission: {}\n")
    _write_text(base / f"{lower}_sequence.fasta", ">seq1\nACGT\n")
    for name in ["fastq_1_R1", "fastq_1_R2", "fastq_2_R1", "fastq_2_R2"]:
        _write_text(base / f"{lower}_{name}.fastq.gz", f"{name}\n")
    pd.DataFrame(
        {
            "bs-sample_name": ["BS1"],
            "sample_name": ["S1"],
            "organism": ["Influenza A virus"],
            "collection_date": ["2024-01-01"],
            "authors": ["Doe, Jane"],
            "bioproject": ["PRJNA1"],
        }
    ).to_csv(base / f"{lower}_biosample_metadata.csv", index=False)
    pd.DataFrame(
        {
            "sra-sample_name": ["SRA1"],
            "sample_name": ["S1"],
            "organism": ["Influenza A virus"],
            "collection_date": ["2024-01-01"],
            "authors": ["Doe, Jane"],
            "bioproject": ["PRJNA1"],
            "bs-sample_name": ["BS1"],
        }
    ).to_csv(base / f"{lower}_sra_metadata.csv", index=False)
    pd.DataFrame(
        {
            "gb-sample_name": ["GB1"],
            "sequence_name": ["seq1"],
            "sample_name": ["S1"],
            "organism": ["Influenza A virus"],
            "collection_date": ["2024-01-01"],
            "authors": ["Doe, Jane"],
            "bioproject": ["PRJNA1"],
            "bs-sample_name": ["BS1"],
        }
    ).to_csv(base / f"{lower}_genbank_metadata.csv", index=False)

#*******************************************************************************
#                          create_test_data
#*******************************************************************************

def test_create_test_data__cov_is_valid_organism(setup_module, tmp_path, monkeypatch):
    prog_dir = tmp_path / "seqsender"
    _make_test_data_tree(prog_dir, "COV")
    monkeypatch.setattr(setup_module, "PROG_DIR", str(prog_dir))
    setup_module.create_test_data("COV", ["BIOSAMPLE"], str(tmp_path))
    assert (tmp_path / "COV_TEST_DATA" / "metadata.csv").is_file()
    assert (tmp_path / "COV_TEST_DATA" / "config.yaml").is_file()

def test_create_test_data__invalid_organism_exits(setup_module, capsys):
    with pytest.raises(SystemExit) as exc:
        setup_module.create_test_data("OTHER", ["BIOSAMPLE"], "/tmp/out")
    assert exc.value.code == 0
    assert capsys.readouterr().out == (
        'SeqSender currently only has test data available for the organisms "FLU" and "COV" '
        "currently, more test sets will be added with later versions. \n"
    )

def test_create_test_data__drops_repeat_columns_from_later_databases_but_keeps_unique_columns(setup_module, tmp_path, monkeypatch):
    prog_dir = tmp_path / "seqsender"
    _make_test_data_tree(prog_dir, "FLU")
    base = prog_dir / "test_data" / "FLU"
    pd.DataFrame(
        {
            "gb-sample_name": ["GB1"],
            "sequence_name": ["seq1"],
            "gb-extra": ["keep-me"],
            "sample_name": ["S1"],
            "organism": ["Influenza A virus"],
            "collection_date": ["2024-01-01"],
            "authors": ["Doe, Jane"],
            "bioproject": ["PRJNA1"],
            "bs-sample_name": ["BS1"],
        }
    ).to_csv(base / "flu_genbank_metadata.csv", index=False)
    monkeypatch.setattr(setup_module, "PROG_DIR", str(prog_dir))
    setup_module.create_test_data("FLU", ["BIOSAMPLE", "GENBANK"], str(tmp_path))
    metadata = pd.read_csv(tmp_path / "FLU_TEST_DATA" / "metadata.csv", dtype=str)
    assert metadata.columns.tolist() == [
        "bs-sample_name",
        "sample_name",
        "organism",
        "collection_date",
        "authors",
        "bioproject",
        "gb-sample_name",
        "sequence_name",
        "gb-extra",
    ]
    assert metadata.loc[0, "gb-extra"] == "keep-me"
    assert metadata.loc[0, "sequence_name"] == "seq1"
    assert metadata.columns.tolist().count("sample_name") == 1
    assert metadata.columns.tolist().count("organism") == 1
    assert metadata.columns.tolist().count("collection_date") == 1
    assert metadata.columns.tolist().count("authors") == 1
    assert metadata.columns.tolist().count("bioproject") == 1
    assert metadata.columns.tolist().count("bs-sample_name") == 1

def test_create_test_data__creates_combined_files_for_selected_databases(setup_module, tmp_path, monkeypatch):
    prog_dir = tmp_path / "seqsender"
    _make_test_data_tree(prog_dir, "FLU")
    monkeypatch.setattr(setup_module, "PROG_DIR", str(prog_dir))
    setup_module.create_test_data("FLU", ["BIOSAMPLE", "SRA", "GENBANK"], str(tmp_path))
    out_dir = tmp_path / "FLU_TEST_DATA"
    assert (out_dir / "config.yaml").is_file()
    assert (out_dir / "sequence.fasta").is_file()
    assert (out_dir / "raw_reads" / "fastq_1_R1.fastq.gz").is_file()
    metadata = pd.read_csv(out_dir / "metadata.csv", dtype=str)
    assert "bs-sample_name" in metadata.columns
    assert "sra-sample_name" in metadata.columns
    assert "gb-sample_name" in metadata.columns

def test_create_test_data__uses_exact_paths_read_csv_merge_copy_and_stdout(setup_module, tmp_path, monkeypatch, capsys):
    prog_dir = tmp_path / "seqsender"
    _make_test_data_tree(prog_dir, "FLU")
    monkeypatch.setattr(setup_module, "PROG_DIR", str(prog_dir))

    read_csv_calls: list[dict[str, Any]] = []
    original_read_csv = setup_module.pd.read_csv

    def fake_read_csv(*args, **kwargs):
        read_csv_calls.append({"args": args, "kwargs": kwargs})
        return original_read_csv(*args, **kwargs)

    merge_calls: list[dict[str, Any]] = []
    original_merge = setup_module.pd.merge

    def fake_merge(*args, **kwargs):
        merge_calls.append({"args": args, "kwargs": kwargs})
        return original_merge(*args, **kwargs)

    copy_calls: list[tuple[str, str]] = []
    original_copy = setup_module.shutil.copy

    def fake_copy(src, dst):
        copy_calls.append((str(src), str(dst)))
        return original_copy(src, dst)

    to_csv_calls: list[dict[str, Any]] = []
    original_to_csv = setup_module.pd.DataFrame.to_csv

    def fake_to_csv(self, path, *args, **kwargs):
        to_csv_calls.append({"columns": self.columns.tolist(), "path": str(path), "args": args, "kwargs": kwargs})
        return original_to_csv(self, path, *args, **kwargs)

    makedirs_calls: list[dict[str, Any]] = []
    original_makedirs = setup_module.os.makedirs

    def fake_makedirs(path, *args, **kwargs):
        makedirs_calls.append({"path": str(path), "args": args, "kwargs": kwargs})
        return original_makedirs(path, *args, **kwargs)

    monkeypatch.setattr(setup_module.pd, "read_csv", fake_read_csv)
    monkeypatch.setattr(setup_module.pd, "merge", fake_merge)
    monkeypatch.setattr(setup_module.shutil, "copy", fake_copy)
    monkeypatch.setattr(setup_module.pd.DataFrame, "to_csv", fake_to_csv)
    monkeypatch.setattr(setup_module.os, "makedirs", fake_makedirs)
    setup_module.create_test_data("FLU", ["BIOSAMPLE", "SRA", "GENBANK"], str(tmp_path))
    out_dir = tmp_path / "FLU_TEST_DATA"
    raw_reads = out_dir / "raw_reads"
    assert makedirs_calls == [{"path": str(out_dir), "args": (), "kwargs": {"exist_ok": True}}, {"path": str(raw_reads), "args": (), "kwargs": {"exist_ok": True}}]
    assert read_csv_calls == [
        {
            "args": (str(prog_dir / "test_data" / "FLU" / "flu_biosample_metadata.csv"),),
            "kwargs": {
                "header": 0,
                "dtype": str,
                "engine": "python",
                "encoding": "utf-8",
                "index_col": False,
                "na_filter": False,
            },
        },
        {
            "args": (str(prog_dir / "test_data" / "FLU" / "flu_sra_metadata.csv"),),
            "kwargs": {
                "header": 0,
                "dtype": str,
                "engine": "python",
                "encoding": "utf-8",
                "index_col": False,
                "na_filter": False,
            },
        },
        {
            "args": (str(prog_dir / "test_data" / "FLU" / "flu_genbank_metadata.csv"),),
            "kwargs": {
                "header": 0,
                "dtype": str,
                "engine": "python",
                "encoding": "utf-8",
                "index_col": False,
                "na_filter": False,
            },
        },
    ]
    assert len(merge_calls) == 2
    assert all(call["kwargs"] == {"how": "left", "left_index": True, "right_index": True} for call in merge_calls)
    assert to_csv_calls[-1]["path"] == str(out_dir / "metadata.csv")
    assert to_csv_calls[-1]["kwargs"] == {"index": False}
    assert "bs-sample_name" in to_csv_calls[-1]["columns"]
    assert "sra-sample_name" in to_csv_calls[-1]["columns"]
    assert "gb-sample_name" in to_csv_calls[-1]["columns"]
    assert "sample_name" in to_csv_calls[-1]["columns"]
    assert to_csv_calls[-1]["columns"].count("sample_name") == 1
    assert copy_calls == [
        (
            str(prog_dir / "test_data" / "FLU" / "flu_config.yaml"),
            str(out_dir / "config.yaml"),
        ),
        (
            str(prog_dir / "test_data" / "FLU" / "flu_sequence.fasta"),
            str(out_dir / "sequence.fasta"),
        ),
        (
            str(prog_dir / "test_data" / "FLU" / "flu_fastq_1_R1.fastq.gz"),
            str(raw_reads / "fastq_1_R1.fastq.gz"),
        ),
        (
            str(prog_dir / "test_data" / "FLU" / "flu_fastq_1_R2.fastq.gz"),
            str(raw_reads / "fastq_1_R2.fastq.gz"),
        ),
        (
            str(prog_dir / "test_data" / "FLU" / "flu_fastq_2_R1.fastq.gz"),
            str(raw_reads / "fastq_2_R1.fastq.gz"),
        ),
        (
            str(prog_dir / "test_data" / "FLU" / "flu_fastq_2_R2.fastq.gz"),
            str(raw_reads / "fastq_2_R2.fastq.gz"),
        ),
    ]
    assert capsys.readouterr().out == (
        "\n"
        "Generating submission test_data\n"
        f"Files are stored at: {out_dir}\n"
    )

#*******************************************************************************
#                        download_table2asn
#*******************************************************************************

def test_download_table2asn__linux_downloads_gzip_and_marks_executable(setup_module, tmp_path, monkeypatch):
    target = tmp_path / "table2asn"
    payload = gzip.compress(b"#!/bin/sh\necho table2asn\n")

    class FakeResponse:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def read(self):
            return payload

    def fake_system():
        return "Linux"

    def fake_urlopen(url):
        return FakeResponse()

    monkeypatch.setattr(setup_module.platform, "system", fake_system)
    monkeypatch.setattr(setup_module.urllib.request, "urlopen", fake_urlopen)
    setup_module.download_table2asn(str(target))
    assert target.read_bytes().startswith(b"#!/bin/sh")
    assert os.access(target, os.X_OK)

@pytest.mark.parametrize(
    ("system_name", "expected_url"),
    [
        ("Linux", "https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/linux64.table2asn.gz"),
        ("Darwin", "https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/mac.table2asn.gz"),
    ],
)
def test_download_table2asn__unix_downloads_exact_url_and_chmods_exact_mode(
    setup_module,
    tmp_path,
    monkeypatch,
    system_name: str,
    expected_url: str,
):
    target = tmp_path / "table2asn"
    payload = gzip.compress(b"binary")
    opened_urls: list[str] = []
    chmod_calls: list[tuple[str, int]] = []

    class FakeResponse:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def read(self):
            return payload

    def fake_system():
        return system_name

    def fake_urlopen(url):
        opened_urls.append(url)
        return FakeResponse()

    original_chmod = setup_module.os.chmod

    def fake_chmod(path, mode):
        chmod_calls.append((str(path), mode))
        return original_chmod(path, mode)

    monkeypatch.setattr(setup_module.platform, "system", fake_system)
    monkeypatch.setattr(setup_module.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(setup_module.os, "chmod", fake_chmod)

    setup_module.download_table2asn(str(target))

    stat_mode = setup_module.os.stat(target).st_mode
    assert target.read_bytes() == b"binary"
    assert opened_urls == [expected_url]
    assert chmod_calls == [(str(target), stat_mode | setup_module.stat.S_IXOTH | setup_module.stat.S_IRWXU)]


def test_download_table2asn__windows_uses_exact_url_and_extracts_to_target_dir(
    setup_module,
    tmp_path,
    monkeypatch,
):
    zip_bytes = io.BytesIO()
    with zipfile.ZipFile(zip_bytes, "w") as zf:
        zf.writestr("table2asn.exe", b"binary")

    opened_urls: list[str] = []

    class FakeResponse:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def read(self):
            return zip_bytes.getvalue()

    def fake_system():
        return "Windows"

    def fake_urlopen(url):
        opened_urls.append(url)
        return FakeResponse()

    monkeypatch.setattr(setup_module.platform, "system", fake_system)
    monkeypatch.setattr(setup_module.urllib.request, "urlopen", fake_urlopen)

    setup_module.download_table2asn(str(tmp_path))

    assert opened_urls == [
        "https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/win64.table2asn.zip"
    ]
    assert (tmp_path / "table2asn.exe").read_bytes() == b"binary"


def test_download_table2asn__unsupported_platform_prints_exact_error(
    setup_module,
    monkeypatch,
    capsys,
):
    def fake_system():
        return "Unsupported"

    monkeypatch.setattr(setup_module.platform, "system", fake_system)

    with pytest.raises(SystemExit) as excinfo:
        setup_module.download_table2asn("/tmp/table2asn")

    assert excinfo.value.code == 1
    assert capsys.readouterr().err == (
        "Error: Cannot identify correct system platform. Please download correct Table2asn version "
        "for system and place it in script directory.\n"
    )

def test_download_table2asn__download_error_prints_exact_stderr(setup_module, tmp_path, monkeypatch, capsys):
    def fake_system():
        return "Linux"

    def raise_error(url):
        raise RuntimeError("network down")

    monkeypatch.setattr(setup_module.platform, "system", fake_system)
    monkeypatch.setattr(setup_module.urllib.request, "urlopen", raise_error)

    with pytest.raises(SystemExit) as excinfo:
        setup_module.download_table2asn(str(tmp_path / "table2asn"))

    assert excinfo.value.code == 1
    assert capsys.readouterr().err == (
        "Downloading table2asn error\n"
        "network down\n"
    )

def test_download_table2asn__windows_extracts_zip(setup_module, tmp_path, monkeypatch):
    zip_bytes = io.BytesIO()
    with zipfile.ZipFile(zip_bytes, "w") as zf:
        zf.writestr("table2asn.exe", b"binary")
    zip_payload = zip_bytes.getvalue()

    class FakeResponse:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def read(self):
            return zip_payload

    def fake_system():
        return "Windows"

    def fake_urlopen(url):
        return FakeResponse()

    monkeypatch.setattr(setup_module.platform, "system", fake_system)
    monkeypatch.setattr(setup_module.urllib.request, "urlopen", fake_urlopen)

    setup_module.download_table2asn(str(tmp_path))

    assert (tmp_path / "table2asn.exe").read_bytes() == b"binary"

def test_download_table2asn__unsupported_platform_exits(setup_module, monkeypatch):
    def fake_system():
        return "Unsupported"

    monkeypatch.setattr(setup_module.platform, "system", fake_system)
    with pytest.raises(SystemExit) as excinfo:
        setup_module.download_table2asn("/tmp/table2asn")
    assert excinfo.value.code == 1

def test_download_table2asn__download_error_exits(setup_module, tmp_path, monkeypatch):
    def fake_system():
        return "Linux"

    monkeypatch.setattr(setup_module.platform, "system", fake_system)

    def raise_error(url):
        raise RuntimeError("network down")

    monkeypatch.setattr(setup_module.urllib.request, "urlopen", raise_error)
    with pytest.raises(SystemExit) as excinfo:
        setup_module.download_table2asn(str(tmp_path / "table2asn"))
    assert excinfo.value.code == 1

#*******************************************************************************
#                               download_xml
#*******************************************************************************

def test_download_xml__sets_utf8_opens_write_plus_and_replaces_nbsp(setup_module, tmp_path, monkeypatch):
    response = SimpleNamespace(text="Alpha\xa0Beta", encoding=None)
    get_calls: list[str] = []
    open_calls: list[tuple[str, str]] = []
    original_open = open

    def fake_get(url):
        get_calls.append(url)
        return response

    def fake_open(file, mode="r", *args, **kwargs):
        open_calls.append((str(file), mode))
        return original_open(file, mode, *args, **kwargs)

    monkeypatch.setattr(setup_module.requests, "get", fake_get)
    monkeypatch.setattr("builtins.open", fake_open)
    output_file = tmp_path / "package.xml"
    setup_module.download_xml("https://example.test/package.xml", str(output_file))
    assert get_calls == ["https://example.test/package.xml"]
    assert response.encoding == "UTF-8"
    assert open_calls == [(str(output_file), "w+")]
    assert output_file.read_text(encoding="utf-8") == "Alpha Beta"

#*******************************************************************************
#                      download_biosample_xml_list
#*******************************************************************************

def test_download_biosample_xml_list__downloads_each_non_generic_package(setup_module, tmp_path, monkeypatch):
    prog_dir = tmp_path / "seqsender"
    biosample_dir = prog_dir / "config" / "biosample"
    biosample_dir.mkdir(parents=True)
    monkeypatch.setattr(setup_module, "PROG_DIR", str(prog_dir))

    downloaded = []
    converted = []
    updated = []

    def fake_download_xml(xml_url: str, output_file: str) -> None:
        downloaded.append((xml_url, output_file))
        if output_file.endswith("biosample_package_list.xml"):
            Path(output_file).write_text(
                "<Packages>\n<Name>Generic.1.0</Name>\n<Name>Pathogen.cl.1.0</Name>\n</Packages>\n",
                encoding="utf-8",
            )
        else:
            Path(output_file).write_text("<BioSamplePackages />", encoding="utf-8")

    def fake_biosample_package_to_pandera_schema(xml_file, name):
        converted.append((Path(xml_file).name, name))

    def fake_sleep(seconds):
        return None

    def fake_update_all_schema_templates():
        updated.append(True)

    monkeypatch.setattr(setup_module, "download_xml", fake_download_xml)
    monkeypatch.setattr(setup_module,"biosample_package_to_pandera_schema", fake_biosample_package_to_pandera_schema)
    monkeypatch.setattr(setup_module.time, "sleep", fake_sleep)
    setup_module.tools.update_all_schema_templates = fake_update_all_schema_templates
    setup_module.download_biosample_xml_list()
    assert any("biosample_package_list.xml" in item[1] for item in downloaded)
    assert converted == [("Pathogen.cl.1.0.xml", "Pathogen.cl.1.0")]
    assert updated == [True]

def test_download_biosample_xml_list__uses_exact_urls_paths_prints_and_error_messages(setup_module, tmp_path, monkeypatch, capsys):
    prog_dir = tmp_path / "seqsender"
    biosample_dir = prog_dir / "config" / "biosample"
    biosample_dir.mkdir(parents=True)
    monkeypatch.setattr(setup_module, "PROG_DIR", str(prog_dir))
    download_calls: list[dict[str, str]] = []
    convert_calls: list[tuple[str, str]] = []
    sleep_calls: list[int] = []
    update_calls: list[bool] = []

    def fake_download_xml(xml_url: str, output_file: str) -> None:
        download_calls.append({"xml_url": xml_url, "output_file": output_file})
        if output_file.endswith("biosample_package_list.xml"):
            Path(output_file).write_text(
                "<Packages>\n"
                "<Name>Generic.1.0</Name>\n"
                "<Name>Pathogen.cl.1.0</Name>\n"
                "<Name>BrokenDownload.1.0</Name>\n"
                "<Name>BrokenConvert.1.0</Name>\n"
                "</Packages>\n",
                encoding="utf-8",
            )
        elif "BrokenDownload" in output_file:
            raise RuntimeError("download boom")
        else:
            Path(output_file).write_text("<BioSamplePackages />", encoding="utf-8")

    def fake_convert(xml_file, name):
        convert_calls.append((xml_file, name))
        if name == "BrokenConvert.1.0":
            raise RuntimeError("convert boom")

    def fake_sleep(seconds):
        sleep_calls.append(seconds)

    def fake_update_all_schema_templates():
        update_calls.append(True)

    monkeypatch.setattr(setup_module, "download_xml", fake_download_xml)
    monkeypatch.setattr(setup_module, "biosample_package_to_pandera_schema", fake_convert)
    monkeypatch.setattr(setup_module.time, "sleep", fake_sleep)
    setup_module.tools.update_all_schema_templates = fake_update_all_schema_templates
    setup_module.download_biosample_xml_list()
    assert download_calls == [
        {
            "xml_url": setup_module.BIOSAMPLE_HTML_PREFIX + setup_module.BIOSAMPLE_HTML_SUFFIX,
            "output_file": str(biosample_dir / "biosample_package_list.xml"),
        },
        {
            "xml_url": setup_module.BIOSAMPLE_HTML_PREFIX + "/Pathogen.cl.1.0" + setup_module.BIOSAMPLE_HTML_SUFFIX,
            "output_file": str(biosample_dir / "Pathogen.cl.1.0.xml"),
        },
        {
            "xml_url": setup_module.BIOSAMPLE_HTML_PREFIX + "/BrokenDownload.1.0" + setup_module.BIOSAMPLE_HTML_SUFFIX,
            "output_file": str(biosample_dir / "BrokenDownload.1.0.xml"),
        },
        {
            "xml_url": setup_module.BIOSAMPLE_HTML_PREFIX + "/BrokenConvert.1.0" + setup_module.BIOSAMPLE_HTML_SUFFIX,
            "output_file": str(biosample_dir / "BrokenConvert.1.0.xml"),
        },
    ]
    assert convert_calls == [
        (str(biosample_dir / "Pathogen.cl.1.0.xml"), "Pathogen.cl.1.0"),
        (str(biosample_dir / "BrokenDownload.1.0.xml"), "BrokenDownload.1.0"),
        (str(biosample_dir / "BrokenConvert.1.0.xml"), "BrokenConvert.1.0"),
    ]
    assert sleep_calls == [1, 1, 1]
    assert update_calls == [True]
    captured = capsys.readouterr()
    assert captured.out == (
        "Downloading Package: Pathogen.cl.1.0\n"
        "Downloading Package: BrokenDownload.1.0\n"
        "Downloading Package: BrokenConvert.1.0\n"
    )
    assert captured.err == (
        "Error: BioSample package BrokenDownload.1.0 failed to download.\n"
        "download boom\n"
        "Error: BioSample package BrokenConvert.1.0 failed to convert to schema.\n"
        "convert boom\n"
    )

#*******************************************************************************
#                     biosample_package_to_pandera_schema
#*******************************************************************************

def test_biosample_package_to_pandera_schema__writes_schema_and_deletes_xml(setup_module, tmp_path, monkeypatch):
    prog_dir = tmp_path / "seqsender"
    biosample_dir = prog_dir / "config" / "biosample"
    biosample_dir.mkdir(parents=True)
    monkeypatch.setattr(setup_module, "PROG_DIR", str(prog_dir))
    xml_file = tmp_path / "package.xml"
    xml_file.write_text(
        """
        <BioSamplePackages>
          <Package>
            <Attribute use="mandatory">
              <HarmonizedName>host</HarmonizedName>
              <Name>host</Name>
              <Description>Host description</Description>
            </Attribute>
            <Attribute use="optional">
              <HarmonizedName>host_sex</HarmonizedName>
              <Name>host sex</Name>
              <Description>Gender or physical sex of host</Description>
              <Format type="select"><Description>male | female</Description></Format>
            </Attribute>
            <Attribute use="either_one_mandatory" group_name="Organism">
              <HarmonizedName>strain</HarmonizedName>
              <Name>strain</Name>
              <Description>Strain name</Description>
            </Attribute>
            <Attribute use="either_one_mandatory" group_name="Organism">
              <HarmonizedName>isolate</HarmonizedName>
              <Name>isolate</Name>
              <Description>Isolate name</Description>
            </Attribute>
            <Attribute use="mandatory">
              <HarmonizedName>collection_date</HarmonizedName>
              <Name>collection date</Name>
              <Description>Reserved field skipped</Description>
            </Attribute>
          </Package>
        </BioSamplePackages>
        """,
        encoding="utf-8",
    )
    setup_module.biosample_package_to_pandera_schema(str(xml_file), "Pathogen.cl.1.0")
    schema_file = biosample_dir / "Pathogen_cl_1_0.py"
    text = schema_file.read_text(encoding="utf-8")
    assert not xml_file.exists()
    assert '"bs-host"' in text
    assert '"bs-host_sex"' in text
    assert "Biological sex" in text
    assert '"bs-collection_date"' not in text
    assert "At least one required: Group" in text
    assert "df[\"bs-strain\"].isnull() & df[\"bs-isolate\"].isnull()" in text

def test_biosample_package_to_pandera_schema__serializes_xml_with_exact_tostring_kwargs(setup_module, tmp_path, monkeypatch):
    prog_dir = tmp_path / "seqsender"
    biosample_dir = prog_dir / "config" / "biosample"
    biosample_dir.mkdir(parents=True)
    monkeypatch.setattr(setup_module, "PROG_DIR", str(prog_dir))
    xml_file = tmp_path / "package.xml"
    xml_file.write_text(
        """
        <BioSamplePackages>
          <Package>
            <Attribute use="mandatory">
              <HarmonizedName>host</HarmonizedName>
              <Name>host</Name>
              <Description>Host description</Description>
            </Attribute>
          </Package>
        </BioSamplePackages>
        """,
        encoding="utf-8",
    )
    tostring_calls: list[dict[str, Any]] = []
    original_tostring = setup_module.ET.tostring

    def fake_tostring(element, *args, **kwargs):
        tostring_calls.append({"tag": element.tag, "args": args, "kwargs": kwargs})
        return original_tostring(element, *args, **kwargs)

    monkeypatch.setattr(setup_module.ET, "tostring", fake_tostring)
    setup_module.biosample_package_to_pandera_schema(str(xml_file), "Pathogen.cl.1.0")
    assert tostring_calls == [{"tag": "BioSamplePackages", "args": (), "kwargs": {"encoding": "utf-8", "method": "xml"}}]

def test_biosample_package_to_pandera_schema__writes_exact_generated_schema_fragments(setup_module, tmp_path, monkeypatch):
    prog_dir = tmp_path / "seqsender"
    biosample_dir = prog_dir / "config" / "biosample"
    biosample_dir.mkdir(parents=True)
    monkeypatch.setattr(setup_module, "PROG_DIR", str(prog_dir))
    xml_file = tmp_path / "package.xml"
    xml_file.write_text(
        """
        <BioSamplePackages>
          <Package>
            <Attribute use="mandatory">
              <HarmonizedName>host</HarmonizedName>
              <Name>host</Name>
              <Description>Host "quoted"
              description</Description>
            </Attribute>
            <Attribute use="optional">
              <HarmonizedName>host_sex</HarmonizedName>
              <Name>host sex</Name>
              <Description>Gender or physical sex of host</Description>
              <Format type="select"><Description>male | female</Description></Format>
            </Attribute>
            <Attribute use="optional">
              <HarmonizedName>geo_loc_name</HarmonizedName>
              <Name>geographic location</Name>
              <Description>Location description</Description>
              <Format type="free_text"><Description>free text</Description></Format>
            </Attribute>
            <Attribute use="either_one_mandatory" group_name="Organism">
              <HarmonizedName>strain</HarmonizedName>
              <Name>strain</Name>
              <Description>Strain name</Description>
            </Attribute>
            <Attribute use="either_one_mandatory" group_name="Organism">
              <HarmonizedName>isolate</HarmonizedName>
              <Name>isolate</Name>
              <Description>Isolate name</Description>
            </Attribute>
            <Attribute use="mandatory">
              <HarmonizedName>collection_date</HarmonizedName>
              <Name>collection date</Name>
              <Description>Reserved field skipped</Description>
            </Attribute>
            <Attribute use="mandatory">
              <HarmonizedName>gender_restroom</HarmonizedName>
              <Name>gender restroom</Name>
              <Description>Reserved field skipped</Description>
            </Attribute>
          </Package>
        </BioSamplePackages>
        """,
        encoding="utf-8",
    )
    setup_module.biosample_package_to_pandera_schema(str(xml_file), "Pathogen.cl.1.0")
    schema_file = biosample_dir / "Pathogen_cl_1_0.py"
    text = schema_file.read_text(encoding="utf-8")
    lines = text.splitlines()
    assert not xml_file.exists()
    assert '"bs-collection_date"' not in text
    assert '"bs-gender_restroom"' not in text
    required_fragments = [
        '"bs-host": Column(',
        'dtype="object",',
        'checks=None,',
        'nullable=False,',
        'unique=False,',
        'coerce=False,',
        'required=True,',
        'Host \\"quoted\\"',
        'description",',
        'title="host",',
        '"bs-host_sex": Column(',
        r'checks=Check.str_matches(r\"(?i)(\W|^)(male|female)(\W|$)\"),',
        'nullable=True,',
        'required=False,',
        'description="Biological sex of host",',
        'title="host sex",',
        '"bs-geo_loc_name": Column(',
        '"bs-strain": Column(',
        'description="At least one required: Group \\"Organism\\". Strain name",',
        '"bs-isolate": Column(',
        'description="At least one required: Group \\"Organism\\". Isolate name",',
        'checks=[',
        'Check(lambda df: ~(df["bs-strain"].isnull() & df["bs-isolate"].isnull()), ignore_na = False),',
        'index=None,',
        'strict="filter",',
        'name="biosample_package_Pathogen.cl.1.0_schema",',
        'ordered=False,',
        'unique=None,',
        'report_duplicates="all",',
        'unique_column_names=True,',
        'add_missing_columns=False,',
        'title="BioSample package Pathogen.cl.1.0 schema",',
        'description="Schema validation for BioSample database using Pathogen.cl.1.0 package.",',
    ]
    for fragment in required_fragments:
        assert fragment in text

    assert text.count('dtype="object",') == 7
    assert text.count('nullable=False,') == 3
    assert text.count('unique=False,') >= 6
    assert text.count('coerce=False,') >= 7
    assert text.count('required=True,') == 5
    assert text.count('required=False,') >= 3

    assert 'DTYPE="OBJECT",' not in text
    assert 'checks=none,' not in text
    assert 'CHECKS=NONE,' not in text
    assert 'nullable=false,' not in text
    assert 'NULLABLE=FALSE,' not in text
    assert 'nullable=true,' not in text
    assert 'NULLABLE=TRUE,' not in text
    assert 'unique=false,' not in text
    assert 'UNIQUE=FALSE,' not in text
    assert 'coerce=false,' not in text
    assert 'COERCE=FALSE,' not in text
    assert 'required=true,' not in text
    assert 'REQUIRED=TRUE,' not in text
    assert 'required=false,' not in text
    assert 'REQUIRED=FALSE,' not in text

    assert 'XX' not in text
    assert "CHECKS=[" not in text
    assert 'strict="filter",' in text
    assert 'unique_column_names=True,' in text
    assert text.rstrip().endswith(")")

def test_biosample_package_to_pandera_schema__non_select_format_writes_checks_none(setup_module, tmp_path, monkeypatch):
    prog_dir = tmp_path / "seqsender"
    biosample_dir = prog_dir / "config" / "biosample"
    biosample_dir.mkdir(parents=True)
    monkeypatch.setattr(setup_module, "PROG_DIR", str(prog_dir))
    xml_file = tmp_path / "package.xml"
    xml_file.write_text(
        """
        <BioSamplePackages>
          <Package>
            <Attribute use="optional">
              <HarmonizedName>free_text_field</HarmonizedName>
              <Name>free text field</Name>
              <Description>Free text description</Description>
              <Format type="free_text"><Description>free text</Description></Format>
            </Attribute>
          </Package>
        </BioSamplePackages>
        """,
        encoding="utf-8",
    )
    setup_module.biosample_package_to_pandera_schema(str(xml_file), "Pathogen.cl.1.0")
    text = (biosample_dir / "Pathogen_cl_1_0.py").read_text(encoding="utf-8")
    assert '"bs-free_text_field": Column(' in text
    field_block = text.split('"bs-free_text_field": Column(', 1)[1].split('),', 1)[0]
    assert 'checks=None,' in field_block
    assert 'Check.str_matches' not in field_block
    assert 'required=False,' in field_block
    assert 'nullable=True,' in field_block

def test_biosample_package_to_pandera_schema__reserved_fields_continue_not_break(setup_module, tmp_path, monkeypatch):
    prog_dir = tmp_path / "seqsender"
    biosample_dir = prog_dir / "config" / "biosample"
    biosample_dir.mkdir(parents=True)
    monkeypatch.setattr(setup_module, "PROG_DIR", str(prog_dir))
    xml_file = tmp_path / "package.xml"
    xml_file.write_text(
        """
        <BioSamplePackages>
          <Package>
            <Attribute use="mandatory">
              <HarmonizedName>collection_date</HarmonizedName>
              <Name>collection date</Name>
              <Description>Reserved field skipped</Description>
            </Attribute>
            <Attribute use="mandatory">
              <HarmonizedName>gender_restroom</HarmonizedName>
              <Name>gender restroom</Name>
              <Description>Reserved field skipped</Description>
            </Attribute>
            <Attribute use="optional">
              <HarmonizedName>host_disease</HarmonizedName>
              <Name>host disease</Name>
              <Description>Disease description</Description>
            </Attribute>
          </Package>
        </BioSamplePackages>
        """,
        encoding="utf-8",
    )
    setup_module.biosample_package_to_pandera_schema(str(xml_file), "Pathogen.cl.1.0")
    text = (biosample_dir / "Pathogen_cl_1_0.py").read_text(encoding="utf-8")
    assert '"bs-collection_date"' not in text
    assert '"bs-gender_restroom"' not in text
    assert '"bs-host_disease": Column(' in text

#*******************************************************************************
#                      test_internet_connection
#*******************************************************************************

def test_test_internet_connection__ncbi_exact_stdout_command_and_subprocess_kwargs(setup_module, monkeypatch, capsys):
    request_urls: list[str] = []
    download_calls: list[str] = []
    run_calls: list[dict[str, Any]] = []
    ftp_events: list[Any] = []

    class FakeFTP:
        def connect(self, host, port, timeout):
            ftp_events.append(("connect", host, port, timeout))

        def quit(self):
            ftp_events.append(("quit",))

    def fake_get(url):
        request_urls.append(url)
        return SimpleNamespace(status_code=200)

    def fake_gethostbyname(host):
        return "1.2.3.4"

    def fake_ftp():
        return FakeFTP()

    def fake_download_table2asn(table2asn_dir):
        download_calls.append(table2asn_dir)

    def fake_run(command, stdout, stderr, cwd):
        run_calls.append({"command": command, "stdout": stdout, "stderr": stderr, "cwd": cwd})
        return SimpleNamespace(returncode=0, stdout=b"", stderr=b"")

    monkeypatch.setattr(setup_module.requests, "get", fake_get)
    monkeypatch.setattr(setup_module.socket, "gethostbyname", fake_gethostbyname)
    monkeypatch.setattr(setup_module.ftplib, "FTP", fake_ftp)
    monkeypatch.setattr(setup_module, "download_table2asn", fake_download_table2asn)
    monkeypatch.setattr(setup_module.subprocess, "run", fake_run)
    setup_module.test_internet_connection(["NCBI"])
    assert request_urls == ["http://www.google.com", "https://www.google.com", "https://www.ncbi.nlm.nih.gov", "https://submit.ncbi.nlm.nih.gov"]
    assert ftp_events == [("connect", setup_module.NCBI_FTP_HOST, 21, 10), ("quit",)]
    assert download_calls == ["/tmp/table2asn"]
    assert run_calls == [
        {
            "command": ["/tmp/table2asn", "-version-full-xml"],
            "stdout": setup_module.subprocess.PIPE,
            "stderr": setup_module.subprocess.PIPE,
            "cwd": os.path.join(os.path.dirname(os.path.abspath(setup_module.__file__))),
        }
    ]
    captured = capsys.readouterr()
    assert captured.err == ""
    assert captured.out == (
        "Checking network settings...\n"
        "Checking HTTP connection...\n"
        "HTTP 'http://www.google.com' connectivity test ok.\n"
        "Checking HTTPS connection...\n"
        "HTTPS 'https://www.google.com' connectivity test ok.\n"
        "Checking NCBI connection...\n"
        "NCBI 'https://www.ncbi.nlm.nih.gov' connectivity test ok.\n"
        "Checking NCBI API connection...\n"
        "NCBI API 'https://submit.ncbi.nlm.nih.gov' connectivity test ok.\n"
        "Checking DNS resolution for FTP site...\n"
        f"DNS resolution test ok. Able to reach ('{setup_module.NCBI_FTP_HOST} -> 1.2.3.4)\n"
        "Checking port status...\n"
        f"{setup_module.NCBI_FTP_HOST} open on port 21.\n"
        "Checking Table2asn functionality...\n"
        "Running Table2asn.\n"
        "No network connection issues detected.\n"
    )

@pytest.mark.parametrize("status_code", [200, 204, 301, 302])
def test_test_internet_connection__general_accepts_exact_success_status_codes(setup_module, monkeypatch, capsys, status_code):

    def fake_get(url):
        return SimpleNamespace(status_code=status_code)

    monkeypatch.setattr(setup_module.requests, "get", fake_get)
    setup_module.test_internet_connection(["GENERAL"])
    captured = capsys.readouterr()
    assert captured.err == ""
    assert "No network connection issues detected.\n" in captured.out

def test_test_internet_connection__general_success(setup_module, monkeypatch, capsys):
    def fake_get(url):
        return SimpleNamespace(status_code=200)

    monkeypatch.setattr(setup_module.requests, "get", fake_get)

    setup_module.test_internet_connection(["GENERAL"])

    captured = capsys.readouterr()
    assert "Checking network settings" in captured.out
    assert "No network connection issues detected" in captured.out
    assert captured.err == ""

def test_test_internet_connection__reports_bad_http_status(setup_module, monkeypatch, capsys):
    def fake_get(url):
        return SimpleNamespace(status_code=503)

    monkeypatch.setattr(setup_module.requests, "get", fake_get)

    setup_module.test_internet_connection(["GENERAL"])

    captured = capsys.readouterr()
    assert "Error code received:'503'" in captured.err


def test_test_internet_connection__ncbi_success_includes_dns_ftp_and_table2asn(setup_module, monkeypatch, capsys):
    ftp_events = []

    class FakeFTP:
        def connect(self, host, port, timeout):
            ftp_events.append((host, port, timeout))

        def quit(self):
            ftp_events.append("quit")

    def fake_get(url):
        return SimpleNamespace(status_code=200)

    def fake_gethostbyname(host):
        return "1.2.3.4"

    def fake_ftp():
        return FakeFTP()

    def fake_download_table2asn(table2asn_dir):
        return None

    run_calls = []

    def fake_run(command, stdout, stderr, cwd):
        run_calls.append({"command": command, "stdout": stdout, "stderr": stderr, "cwd": cwd})
        return SimpleNamespace(returncode=0, stdout=b"", stderr=b"")

    monkeypatch.setattr(setup_module.requests, "get", fake_get)
    monkeypatch.setattr(setup_module.socket, "gethostbyname", fake_gethostbyname)
    monkeypatch.setattr(setup_module.ftplib, "FTP", fake_ftp)
    monkeypatch.setattr(setup_module, "download_table2asn", fake_download_table2asn)
    monkeypatch.setattr(setup_module.subprocess, "run", fake_run)
    setup_module.test_internet_connection(["NCBI"])
    captured = capsys.readouterr()
    assert "DNS resolution test ok" in captured.out
    assert ftp_events == [(setup_module.NCBI_FTP_HOST, 21, 10), "quit"]
    assert "No network connection issues detected" in captured.out
    assert run_calls == [
        {
            "command": ["/tmp/table2asn", "-version-full-xml"],
            "stdout": setup_module.subprocess.PIPE,
            "stderr": setup_module.subprocess.PIPE,
            "cwd": os.path.join(os.path.dirname(os.path.abspath(setup_module.__file__))),
        }
    ]

def test_test_internet_connection__ncbi_collects_dns_ftp_and_table2asn_errors(setup_module, monkeypatch, capsys):
    class FakeFTP:
        def connect(self, host, port, timeout):
            raise OSError("blocked")

    def fake_get(url):
        return SimpleNamespace(status_code=200)

    def fake_gethostbyname(host):
        return ""

    def fake_ftp():
        return FakeFTP()

    def fake_download_table2asn(table2asn_dir):
        raise RuntimeError("download failed")

    monkeypatch.setattr(setup_module.requests, "get", fake_get)
    monkeypatch.setattr(setup_module.socket, "gethostbyname", fake_gethostbyname)
    monkeypatch.setattr(setup_module.ftplib, "FTP", fake_ftp)
    monkeypatch.setattr(setup_module, "download_table2asn", fake_download_table2asn)

    setup_module.test_internet_connection(["NCBI"])

    captured = capsys.readouterr()
    assert "Unable to resolve address" in captured.err
    assert "Unable to download latest version of Table2asn" in captured.err


def test_test_internet_connection__table2asn_nonzero_reports_stdout_stderr(setup_module, monkeypatch, capsys):
    class FakeFTP:
        def connect(self, host, port, timeout):
            return None

        def quit(self):
            return None

    def fake_get(url):
        return SimpleNamespace(status_code=200)

    def fake_gethostbyname(host):
        return "1.2.3.4"

    def fake_ftp():
        return FakeFTP()

    def fake_download_table2asn(table2asn_dir):
        return None

    def fake_run(*args, **kwargs):
        return SimpleNamespace(returncode=2, stdout=b"bad out", stderr=b"bad err")

    monkeypatch.setattr(setup_module.requests, "get", fake_get)
    monkeypatch.setattr(setup_module.socket, "gethostbyname", fake_gethostbyname)
    monkeypatch.setattr(setup_module.ftplib, "FTP", fake_ftp)
    monkeypatch.setattr(setup_module, "download_table2asn",fake_download_table2asn)
    monkeypatch.setattr(setup_module.subprocess, "run", fake_run)
    setup_module.test_internet_connection(["NCBI"])
    captured = capsys.readouterr()
    assert captured.err == (
        "Table2asn-Error\n"
        "b'bad out'\n"
        "b'bad err'\n"
    )

def test_test_internet_connection__http_exception_reports_exact_website_and_exception(setup_module, monkeypatch, capsys):
    def fake_get(url):
        raise RuntimeError("network boom")

    monkeypatch.setattr(setup_module.requests, "get", fake_get)
    setup_module.test_internet_connection(["GENERAL"])
    captured = capsys.readouterr()
    assert captured.out == (
        "Checking network settings...\n"
        "Checking HTTP connection...\n"
        "Checking HTTPS connection...\n"
    )
    assert captured.err == (
        "HTTP connectivity test failed for 'http://www.google.com'. Check possible firewall issues. \n"
        "Exception:network boom\n"
        "HTTPS connectivity test failed for 'https://www.google.com'. Check possible firewall issues. \n"
        "Exception:network boom\n"
    )
