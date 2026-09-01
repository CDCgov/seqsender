import importlib.util
import os
import sys
import types
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import Any

def _source_root() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        mutant_src = parent / "src"
        if parent.name == "mutants" and (mutant_src / "genbank_handler.py").exists():
            return mutant_src

        normal_src = parent / "src"
        if (normal_src / "genbank_handler.py").exists() and parent.name != "mutants":
            return normal_src

    raise RuntimeError("Could not find src/genbank_handler.py")

SOURCE_DIR = _source_root()
MODULE_PATH = SOURCE_DIR / "genbank_handler.py"
if str(SOURCE_DIR) not in sys.path:
    sys.path.insert(0, str(SOURCE_DIR))

"""Unit tests initially generated with ChatGPT v5.5 "Deep Research" and "Thinking"
then modified for use and clarity with validation and coverage testing via mutmut."""

#*******************************************************************************
#                 create test genbank_handler.py connections
#*******************************************************************************

@pytest.fixture()
def genbank_handler_module(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.syspath_prepend(str(SOURCE_DIR))
    src_pkg: Any = types.ModuleType("src")
    src_pkg.__path__ = [str(SOURCE_DIR)]
    monkeypatch.setitem(sys.modules, "src", src_pkg)

    def alias_src_module(name, module):
        monkeypatch.setitem(sys.modules, f"src.{name}", module)
        monkeypatch.setitem(sys.modules, name, module)
        setattr(src_pkg, name, module)

    settings_stub: Any = types.ModuleType("settings")
    settings_stub.NCBI_API_URL = "https://example.test/files/FILE_ID"
    settings_stub.GENBANK_REGEX_SRC = (
        r"^gb-sample_name$|^src-|^bioproject$|^organism$|^collection_date$"
    )
    settings_stub.GENBANK_REGEX_CMT = r"^gb-sample_name$|^cmt-"

    setup_stub: Any = types.ModuleType("setup")

    def download_table2asn(table2asn_dir):
        return None

    setup_stub.download_table2asn = download_table2asn

    file_handler_stub: Any = types.ModuleType("file_handler")
    file_handler_stub.saved_csvs = []
    file_handler_stub.saved_xmls = []
    file_handler_stub.copied_files = []
    file_handler_stub.create_fasta_calls = []

    def save_csv(df, file_path, file_name=None, sep=","):
        file_handler_stub.saved_csvs.append(
            {"df": df.copy(), "file_path": file_path, "file_name": file_name, "sep": sep}
        )
        Path(file_path).mkdir(parents=True, exist_ok=True)
        output = Path(file_path) / file_name if file_name else Path(file_path)
        df.to_csv(output, header=True, index=False, sep=sep)

    def save_xml(submission_xml, submission_dir):
        file_handler_stub.saved_xmls.append(
            {"submission_xml": submission_xml, "submission_dir": submission_dir}
        )
        Path(submission_dir).mkdir(parents=True, exist_ok=True)
        (Path(submission_dir) / "submission.xml").write_bytes(submission_xml)

    def create_fasta(database, metadata, submission_dir, config_dict):
        records = []
        for _, row in metadata.iterrows():
            records.append(
                SeqRecord(
                    Seq(str(row.get("fasta_sequence_orig", "ACGT"))),
                    id=str(row["gb-sample_name"]),
                    description="",
                )
            )
        Path(submission_dir).mkdir(parents=True, exist_ok=True)
        with open(Path(submission_dir) / "sequence.fsa", "w") as handle:
            SeqIO.write(records, handle, "fasta")
        file_handler_stub.create_fasta_calls.append({"database": database,"metadata": metadata.copy(),"submission_dir": submission_dir,"config_dict": config_dict})

    def copy_file(source, destination):
        file_handler_stub.copied_files.append((source, destination))
        Path(destination).write_text(Path(source).read_text() if Path(source).exists() else "")

    def load_csv(file_path, sep=","):
        return pd.read_csv(file_path, header=0, dtype=str, sep=sep, engine="python", na_filter=False)

    file_handler_stub.save_csv = save_csv
    file_handler_stub.save_xml = save_xml
    file_handler_stub.create_fasta = create_fasta
    file_handler_stub.copy_file = copy_file
    file_handler_stub.load_csv = load_csv

    ncbi_stub: Any = types.ModuleType("ncbi_handler")

    def process_report_header(report_file):
        return ({"SubmissionStatus": {}}, "SUBMITTED", "SUB123")

    ncbi_stub.process_report_header = process_report_header

    upload_log_stub: Any = types.ModuleType("upload_log")
    upload_log_stub.updated = []

    def update_submission_status_csv(**kwargs):
        upload_log_stub.updated.append(kwargs)

    upload_log_stub.update_submission_status_csv = update_submission_status_csv

    nameparser_stub: Any = types.ModuleType("nameparser")

    class HumanName:
        def __init__(self, value):
            self.original = value
            self.title = ""
            self.first = ""
            self.middle = ""
            self.last = ""
            self.suffix = ""
            parts = [part.strip() for part in str(value).split(",")]
            if len(parts) >= 2:
                self.last = parts[0]
                given = parts[1].split()
                self.first = given[0] if given else ""
                self.middle = " ".join(given[1:]) if len(given) > 1 else ""
                self.suffix = parts[2].strip() if len(parts) > 2 else ""
            else:
                words = str(value).split()
                self.first = words[0] if words else ""
                self.last = words[-1] if len(words) > 1 else ""
                self.middle = " ".join(words[1:-1]) if len(words) > 2 else ""

    nameparser_stub.HumanName = HumanName

    for module_name in [
        "genbank_handler",
        "src.genbank_handler",
        "settings",
        "src.settings",
        "setup",
        "src.setup",
        "file_handler",
        "src.file_handler",
        "ncbi_handler",
        "src.ncbi_handler",
        "upload_log",
        "src.upload_log",
        "nameparser",
    ]:
        sys.modules.pop(module_name, None)

    monkeypatch.setitem(sys.modules, "src", src_pkg)
    alias_src_module("settings", settings_stub)
    alias_src_module("setup", setup_stub)
    alias_src_module("file_handler", file_handler_stub)
    alias_src_module("ncbi_handler", ncbi_stub)
    alias_src_module("upload_log", upload_log_stub)
    monkeypatch.setitem(sys.modules, "nameparser", nameparser_stub)

    spec = importlib.util.spec_from_file_location("genbank_handler", MODULE_PATH)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules["genbank_handler"] = module
    setattr(src_pkg, "genbank_handler", module)
    spec.loader.exec_module(module)
    return module

#*******************************************************************************
#                      create fake ncbi config file
#*******************************************************************************

@pytest.fixture()
def ncbi_config():
    return {
        "Spuid_Namespace": "CDCNS",
        "Specified_Release_Date": "",
        "GenBank_Auto_Remove_Failed_Samples": False,
        "Publication_Title": "Default publication title",
        "Publication_Status": "Unpublished",
        "Description": {
            "Organization": {
                "Type": "center",
                "Role": "owner",
                "Name": "CDC",
                "Submitter": {
                    "Name": {"First": "Ada", "Last": "Lovelace"},
                    "Email": "ada@example.org",
                    "Alt_Email": "",
                },
                "Address": {
                    "Affil": "CDC",
                    "Div": "Division",
                    "Street": "1 Main St",
                    "City": "Atlanta",
                    "Sub": "GA",
                    "Country": "USA",
                    "Email": "lab@example.org",
                    "Phone": "404-555-1212",
                    "Postal_Code": 30329,
                },
            }
        },
    }

#*******************************************************************************
#                       create fake metadata file
#*******************************************************************************

@pytest.fixture()
def genbank_metadata():
    return pd.DataFrame(
        [
            {
                "organism": "Influenza A virus",
                "authors": "Lovelace, Ada; Hopper, Grace",
                "collection_date": "2024-01-02",
                "bioproject": "PRJNA123",
                "sequence_name": "seq1",
                "gb-sample_name": "GB001",
                "gb-title": "Custom portal title",
                "gb-comment": "Custom portal comment",
                "gb-fasta_definition_line_modifiers": "[country=USA]",
                "src-geo_loc_name": "USA: Georgia",
                "src-Host": "Homo sapiens",
                "src-Isolate": "iso1",
                "cmt-StructuredCommentPrefix": "Assembly-Data",
                "cmt-Assembly Method": "iVar v1",
                "cmt-StructuredCommentSuffix": "Assembly-Data",
                "fasta_sequence_orig": "ACGTACGT",
            }
        ]
    )


def parse_xml(xml_bytes):
    return ET.fromstring(xml_bytes.decode("utf-8"))

#*******************************************************************************
#                         create_submission_xml
#*******************************************************************************

def test_create_submission_xml__blank_title_and_comment_use_defaults(genbank_handler_module, ncbi_config, genbank_metadata):
    metadata = genbank_metadata.copy()
    metadata["gb-title"] = "   "
    metadata["gb-comment"] = "   "
    xml = genbank_handler_module.create_submission_xml(organism="FLU", submission_name="sub1", config_dict=ncbi_config, metadata=metadata)
    root = parse_xml(xml)
    assert root.findtext("./Description/Title") == "sub1-GB"
    assert root.findtext("./Description/Comment") == "GenBank Submission"
    assert root.find("./Action/AddFiles") is not None
    assert root.find("./Action/AddFiles").tag == "AddFiles"

def test_create_submission_xml__has_exact_required_xml_structure(genbank_handler_module, ncbi_config, genbank_metadata):
    xml = genbank_handler_module.create_submission_xml(organism="FLU", submission_name="sub1", config_dict=ncbi_config, metadata=genbank_metadata)
    assert xml.startswith(b"<?xml version='1.0' encoding='utf-8'?>\n")
    root = parse_xml(xml)
    assert root.tag == "Submission"

    description = root.find("./Description")
    assert description is not None
    assert description.findtext("Title") == "Custom portal title"
    assert description.findtext("Comment") == "Custom portal comment"

    organization = description.find("Organization")
    assert organization is not None
    assert organization.attrib == {"type": "center", "role": "owner"}
    assert organization.findtext("Name") == "CDC"

    action = root.find("./Action")
    assert action is not None

    addfiles = action.find("./AddFiles")
    assert addfiles is not None
    assert addfiles.attrib == {"target_db": "GenBank"}

    file_el = addfiles.find("./File")
    assert file_el is not None
    assert file_el.attrib == {"file_path": "sub1.zip"}
    assert file_el.findtext("DataType") == "genbank-submission-package"

    attrs = [(el.attrib, el.text) for el in addfiles.findall("./Attribute")]
    assert attrs == [({"name": "wizard"}, "BankIt_influenza_api"),({"name": "auto_remove_failed_seqs"}, "no")]

    identifier = addfiles.find("./Identifier")
    assert identifier is not None
    spuid = identifier.find("./SPUID")
    assert spuid is not None
    assert spuid.text == "sub1"
    assert spuid.attrib == {"spuid_namespace": "CDCNS"}

def test_create_submission_xml__flu_custom_title_comment_hold_and_auto_remove(genbank_handler_module, ncbi_config, genbank_metadata):
    ncbi_config = dict(ncbi_config)
    ncbi_config["Specified_Release_Date"] = "2030-01-01"
    ncbi_config["GenBank_Auto_Remove_Failed_Samples"] = True

    xml = genbank_handler_module.create_submission_xml(
        organism="FLU", submission_name="sub1", config_dict=ncbi_config, metadata=genbank_metadata
    )

    root = parse_xml(xml)
    assert root.findtext("./Description/Title") == "Custom portal title"
    assert root.findtext("./Description/Comment") == "Custom portal comment"
    assert root.find("./Description/Hold").attrib["release_date"] == "2030-01-01"
    assert root.find(".//AddFiles").attrib["target_db"] == "GenBank"
    attrs = {el.attrib["name"]: el.text for el in root.findall(".//Attribute")}
    assert attrs["wizard"] == "BankIt_influenza_api"
    assert attrs["auto_remove_failed_seqs"] == "yes"
    assert root.find(".//SPUID").text == "sub1"
    assert root.find(".//SPUID").attrib["spuid_namespace"] == "CDCNS"

@pytest.mark.parametrize(("organism", "expected_wizard"),
    [("COV", "BankIt_SARSCoV2_api"), ("OTHER", None), ("FLU", "BankIt_influenza_api")])
def test_create_submission_xml__default_values_and_wizard_variants(genbank_handler_module, ncbi_config, genbank_metadata, organism, expected_wizard):
    metadata = genbank_metadata.drop(columns=["gb-title", "gb-comment"])
    xml = genbank_handler_module.create_submission_xml(organism=organism, submission_name="sub1", config_dict=ncbi_config, metadata=metadata)

    root = parse_xml(xml)
    assert root.findtext("./Description/Title") == "sub1-GB"
    assert root.findtext("./Description/Comment") == "GenBank Submission"
    attrs = {el.attrib["name"]: el.text for el in root.findall(".//Attribute")}
    assert attrs["auto_remove_failed_seqs"] == "no"
    assert attrs.get("wizard") == expected_wizard

def test_create_submission_xml__serializes_with_exact_tostring_options(monkeypatch, genbank_handler_module, ncbi_config, genbank_metadata):
    original_tostring = genbank_handler_module.etree.tostring
    observed: list[dict[str, Any]] = []

    def fake_tostring(element, *args, **kwargs):
        observed.append(
            {
                "tag": element.tag,
                "args": args,
                "kwargs": kwargs,
            }
        )
        return original_tostring(element, *args, **kwargs)

    monkeypatch.setattr(genbank_handler_module.etree, "tostring", fake_tostring)
    genbank_handler_module.create_submission_xml(organism="FLU",submission_name="sub1",config_dict=ncbi_config,metadata=genbank_metadata)
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
#*******************************************************************************
#                              create_authorset
#*******************************************************************************
def normalize_newlines(text: str):
    return text.replace("\r\n", "\n")

def test_create_authorset__writes_middle_names_and_omits_blank_optional_contact_fields(tmp_path, genbank_handler_module, ncbi_config):
    ncbi_config["Description"]["Organization"]["Submitter"]["Alt_Email"] = ""
    ncbi_config["Description"]["Organization"]["Address"]["Phone"] = ""
    metadata = pd.DataFrame([{"authors": "Lovelace, Ada Byron; Hopper, Grace Murray"}])
    genbank_handler_module.create_authorset(config_dict=ncbi_config, metadata=metadata, submission_name="sub1", submission_dir=str(tmp_path), publication_title=None, publication_status=None)
    text = normalize_newlines((tmp_path / "authorset.sbt").read_text())
    assert 'email "lab@example.org",' in text
    assert 'phone "' not in text
    assert "ALT EMAIL:" not in text
    assert text.count('last "Lovelace"') == 3
    assert text.count('first "Ada"') == 3
    assert text.count('middle "Byron"') == 2
    assert text.count('last "Hopper"') == 2
    assert text.count('first "Grace"') == 2
    assert text.count('middle "Murray"') == 2
    assert 'title "Default publication title"' in text
    assert 'cit "Unpublished"' in text

def test_create_authorset__writes_exact_authorset_with_alt_email_and_phone(tmp_path,genbank_handler_module,ncbi_config,genbank_metadata):
    ncbi_config["Description"]["Organization"]["Submitter"]["Alt_Email"] = "alt@example.org"
    genbank_handler_module.create_authorset(config_dict=ncbi_config, metadata=genbank_metadata, submission_name="sub1", submission_dir=str(tmp_path), publication_title="Override title", publication_status="Published")
    text = normalize_newlines((tmp_path / "authorset.sbt").read_text())
    expected = """Submit-block ::= {
  contact {
    contact {
      name name {
        last "Lovelace",
        first "Ada"
      },
      affil std {
        affil "CDC",
        div "Division",
        city "Atlanta",
        sub "GA",
        country "USA",
        street "1 Main St",
        email "lab@example.org",
        phone "404-555-1212",
        postal-code "30329"
      }
    }
  },
  cit {
    authors {
      names std {
        {
          name name {
            last "Lovelace",
            first "Ada"
          }
        },
        {
          name name {
            last "Hopper",
            first "Grace"
          }
        }
      },
      affil std {
        affil "CDC",
        div "Division",
        city "Atlanta",
        sub "GA",
        country "USA",
        street "1 Main St",
        postal-code "30329"
      }
    }
  },
  subtype new
}
Seqdesc ::= pub {
  pub {
    gen {
      cit "Published",
      authors {
        names std {
          {
            name name {
              last "Lovelace",
              first "Ada"
            }
          },
          {
            name name {
              last "Hopper",
              first "Grace"
            }
          }
        }
      },
      title "Override title"
    }
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "AdditionalComment",
      data str "ALT EMAIL: alt@example.org"
    }
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "AdditionalComment",
      data str "Submission Title: sub1"
    }
  }
}
"""

    assert text == expected

def test_create_authorset__semicolon_delimiter_splits_exactly_without_warning(tmp_path, genbank_handler_module, ncbi_config, capsys):
    metadata = pd.DataFrame([{"authors": "Alpha, Ann; Beta, Bob; Gamma, Gail"}])
    genbank_handler_module.create_authorset(config_dict=ncbi_config, metadata=metadata, submission_name="sub1", submission_dir=str(tmp_path), publication_title=None, publication_status=None)
    captured = capsys.readouterr()
    assert captured.err == ""
    text = normalize_newlines((tmp_path / "authorset.sbt").read_text())
    assert text.count('last "Alpha"') == 2
    assert text.count('last "Beta"') == 2
    assert text.count('last "Gamma"') == 2
    assert 'last "Alpha; Beta; Gamma"' not in text

def test_create_authorset__comma_delimited_authors_warns_and_omits_empty_phone_alt_email(tmp_path, genbank_handler_module, ncbi_config, capsys):
    ncbi_config["Description"]["Organization"]["Address"]["Phone"] = ""
    metadata = pd.DataFrame([{"authors": "Lovelace, Ada, Hopper, Grace"}])
    genbank_handler_module.create_authorset(config_dict=ncbi_config, metadata=metadata, submission_name="sub2", submission_dir=str(tmp_path), publication_title=None, publication_status=None)
    captured = capsys.readouterr()
    assert captured.err == "Warning: Metadata 'authors' field is not using a semi-colon ';' to separate individual names. Comma was detected instead and is being used to separate individual names. If this causes issues, separate each name with a semi-colon instead.\n"
    text = normalize_newlines((tmp_path / "authorset.sbt").read_text())
    assert 'phone "' not in text
    assert "ALT EMAIL:" not in text
    assert 'title "Default publication title"' in text
    assert 'cit "Unpublished"' in text
    assert 'postal-code "30329"' in text
    assert 'last "Lovelace"' in text
    assert 'first "Ada"' in text
    assert 'first "Hopper"' in text
    assert 'first "Grace"' in text
    assert 'last "Lovelace, Ada, Hopper, Grace"' not in text

#*******************************************************************************
#                             create_files
#*******************************************************************************

def test_create_files__writes_authorset_source_comment_and_copies_gff(tmp_path, genbank_handler_module, ncbi_config, genbank_metadata):
    gff = tmp_path / "input.gff"
    gff.write_text("##gff-version 3\n")
    genbank_handler_module.create_files(organism="FLU", config_dict=ncbi_config, metadata=genbank_metadata, submission_name="sub1", submission_dir=str(tmp_path), gff_file=str(gff), publication_title=None, publication_status=None)
    assert (tmp_path / "authorset.sbt").exists()
    assert (tmp_path / "sequence.fsa").exists()
    source = pd.read_csv(tmp_path / "source.src", sep="\t", dtype=str)
    assert source.columns.tolist() == ["Sequence_ID","organism","Collection_date","BioProject","geo_loc_name","Host","Isolate"]
    assert source.to_dict("records") == [
        {
            "Sequence_ID": "GB001",
            "organism": "Influenza A virus",
            "Collection_date": "2024-01-02",
            "BioProject": "PRJNA123",
            "geo_loc_name": "USA: Georgia",
            "Host": "Homo sapiens",
            "Isolate": "iso1",
        }
    ]
    comment = pd.read_csv(tmp_path / "comment.cmt", sep="\t", dtype=str)
    assert comment.columns.tolist() == ["SeqID", "StructuredCommentPrefix", "Assembly Method", "StructuredCommentSuffix"]
    assert comment.to_dict("records") == [
        {
            "SeqID": "GB001",
            "StructuredCommentPrefix": "Assembly-Data",
            "Assembly Method": "iVar v1",
            "StructuredCommentSuffix": "Assembly-Data",
        }
    ]
    assert (tmp_path / "sub1.gff").exists()

def test_create_files__forwards_publication_overrides_and_preserves_comment_order(tmp_path, monkeypatch, genbank_handler_module, ncbi_config, genbank_metadata):
    calls: list[dict[str, Any]] = []

    def fake_create_authorset(**kwargs):
        calls.append(kwargs)
        Path(kwargs["submission_dir"]).mkdir(parents=True, exist_ok=True)
        Path(kwargs["submission_dir"], "authorset.sbt").write_text("authorset", encoding="utf-8")

    monkeypatch.setattr(genbank_handler_module, "create_authorset", fake_create_authorset)
    genbank_handler_module.create_files(organism="FLU", config_dict=ncbi_config, metadata=genbank_metadata, submission_name="sub1", submission_dir=str(tmp_path), gff_file=None, publication_title="Exact title", publication_status="Published")
    assert len(calls) == 1
    assert calls[0]["config_dict"] is ncbi_config
    assert calls[0]["submission_name"] == "sub1"
    assert calls[0]["submission_dir"] == str(tmp_path)
    assert calls[0]["publication_title"] == "Exact title"
    assert calls[0]["publication_status"] == "Published"
    assert calls[0]["metadata"].columns.tolist() == genbank_metadata.drop(columns=["gb-title", "gb-comment"]).columns.tolist()
    assert calls[0]["metadata"].to_dict("records") == genbank_metadata.drop(columns=["gb-title", "gb-comment"]).to_dict("records")
    saved_comment = [call for call in genbank_handler_module.file_handler.saved_csvs if call["file_name"] == "comment.cmt"][0]["df"]
    assert saved_comment.columns.tolist() == ["SeqID", "StructuredCommentPrefix", "Assembly Method", "StructuredCommentSuffix"]

def test_create_files__skips_comment_and_gff_when_absent(tmp_path, genbank_handler_module, ncbi_config, genbank_metadata):
    metadata = genbank_metadata.drop(columns=["cmt-StructuredCommentPrefix", "cmt-Assembly Method", "cmt-StructuredCommentSuffix"])

    genbank_handler_module.create_files(
        organism="FLU",
        config_dict=ncbi_config,
        metadata=metadata,
        submission_name="sub1",
        submission_dir=str(tmp_path),
        gff_file=None,
        publication_title=None,
        publication_status=None,
    )

    assert not (tmp_path / "comment.cmt").exists()
    assert not (tmp_path / "sub1.gff").exists()

def test_create_files__drops_portal_xml_columns_and_forwards_exact_arguments(tmp_path, genbank_handler_module, ncbi_config, genbank_metadata):
    genbank_handler_module.create_files(
        organism="FLU",
        config_dict=ncbi_config,
        metadata=genbank_metadata,
        submission_name="sub1",
        submission_dir=str(tmp_path),
        gff_file=None,
        publication_title="Exact title",
        publication_status="Published",
    )

    create_fasta_call = genbank_handler_module.file_handler.create_fasta_calls[0]
    assert create_fasta_call["database"] == "GENBANK"
    assert create_fasta_call["submission_dir"] == str(tmp_path)
    assert create_fasta_call["config_dict"] is ncbi_config

    forwarded_columns = create_fasta_call["metadata"].columns.tolist()
    assert "gb-title" not in forwarded_columns
    assert "gb-comment" not in forwarded_columns
    assert "gb-sample_name" in forwarded_columns
    assert "cmt-StructuredCommentPrefix" in forwarded_columns

#*******************************************************************************
#                               create_zip
#*******************************************************************************

def test_create_zip__writes_exact_paths_and_waits_for_zip(tmp_path, monkeypatch, genbank_handler_module):
    for name in ["authorset.sbt", "sequence.fsa", "source.src", "comment.cmt"]:
        (tmp_path / name).write_text(name)

    write_calls: list[tuple[str, str]] = []
    zip_init_calls: list[tuple[str, str]] = []
    isfile_calls: list[str] = []
    sleep_calls: list[int] = []

    class FakeZip:
        def __init__(self, path, mode):
            zip_init_calls.append((path, mode))

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def write(self, filename, arcname=None):
            write_calls.append((filename, arcname))

    def fake_isfile(path):
        isfile_calls.append(path)
        if path == os.path.join(str(tmp_path), "comment.cmt"):
            return True
        if path == os.path.join(str(tmp_path), "sub1.zip"):
            return len([p for p in isfile_calls if p == path]) > 1
        return False

    def fake_sleep(seconds):
        sleep_calls.append(seconds)

    monkeypatch.setattr(genbank_handler_module, "ZipFile", FakeZip)
    monkeypatch.setattr(genbank_handler_module.os.path, "isfile", fake_isfile)
    monkeypatch.setattr(genbank_handler_module.time, "sleep", fake_sleep)
    genbank_handler_module.create_zip("sub1", str(tmp_path))
    assert zip_init_calls == [(os.path.join(str(tmp_path), "sub1.zip"), "w")]
    assert write_calls == [
        (os.path.join(str(tmp_path), "authorset.sbt"), "authorset.sbt"),
        (os.path.join(str(tmp_path), "sequence.fsa"), "sequence.fsa"),
        (os.path.join(str(tmp_path), "source.src"), "source.src"),
        (os.path.join(str(tmp_path), "comment.cmt"), "comment.cmt"),
    ]
    assert isfile_calls == [
        os.path.join(str(tmp_path), "comment.cmt"),
        os.path.join(str(tmp_path), "sub1.zip"),
        os.path.join(str(tmp_path), "sub1.zip"),
    ]
    assert sleep_calls == [10]

def test_create_zip__includes_required_files_and_optional_comment(tmp_path, genbank_handler_module):
    for name in ["authorset.sbt", "sequence.fsa", "source.src", "comment.cmt"]:
        (tmp_path / name).write_text(name)

    genbank_handler_module.create_zip("sub1", str(tmp_path))

    with zipfile.ZipFile(tmp_path / "sub1.zip") as zf:
        assert sorted(zf.namelist()) == [
            "authorset.sbt",
            "comment.cmt",
            "sequence.fsa",
            "source.src",
        ]

def test_create_zip__omits_comment_when_absent(tmp_path, genbank_handler_module):
    for name in ["authorset.sbt", "sequence.fsa", "source.src"]:
        (tmp_path / name).write_text(name)

    genbank_handler_module.create_zip("sub1", str(tmp_path))

    with zipfile.ZipFile(tmp_path / "sub1.zip") as zf:
        assert sorted(zf.namelist()) == ["authorset.sbt", "sequence.fsa", "source.src"]

#*******************************************************************************
#                            create_table2asn
#*******************************************************************************

def test_create_table2asn__downloads_when_missing_runs_with_comment_and_gff(tmp_path, monkeypatch, genbank_handler_module, capsys):
    (tmp_path / "authorset.sbt").write_text("sbt")
    (tmp_path / "sequence.fsa").write_text(">GB001\nACGT\n")
    (tmp_path / "source.src").write_text("Sequence_ID\nGB001\n")
    (tmp_path / "comment.cmt").write_text("SeqID\nGB001\n")
    (tmp_path / "sub1.gff").write_text("##gff-version 3\n")
    (tmp_path / "sub1.val").write_text("valid\n")
    original_isfile = os.path.isfile

    def fake_isfile(path):
        if path == "/tmp/table2asn":
            return False
        return original_isfile(path)

    monkeypatch.setattr(genbank_handler_module.os.path, "isfile", fake_isfile)
    downloaded = []

    def fake_download_table2asn(table2asn_dir):
        downloaded.append(table2asn_dir)

    genbank_handler_module.setup.download_table2asn = fake_download_table2asn
    runs = []

    def fake_run(command, stdout, stderr, cwd):
        runs.append({"command": command, "stdout": stdout, "stderr": stderr, "cwd": cwd})
        return types.SimpleNamespace(returncode=0, stdout=b"", stderr=b"")

    validation_calls: list[str] = []

    def fake_check_table2asn_submission(validation_file):
        validation_calls.append(validation_file)
        return "VALIDATED"

    monkeypatch.setattr(genbank_handler_module, "check_table2asn_submission", fake_check_table2asn_submission)
    monkeypatch.setattr(genbank_handler_module.subprocess, "run", fake_run)

    result = genbank_handler_module.create_table2asn("sub1", str(tmp_path))
    assert result == "VALIDATED"
    assert validation_calls == [os.path.join(str(tmp_path), "sub1.val")]
    captured = capsys.readouterr()
    assert captured.out == (
        "Downloading Table2asn.\n"
        "Running Table2asn.\n"
        "Validating Table2asn submission.\n"
    )
    assert captured.err == ""
    expected_command = [
        "/tmp/table2asn","-V","vb","-a","s","-t",
        os.path.join(str(tmp_path), "authorset.sbt"),
        "-i",
        os.path.join(str(tmp_path), "sequence.fsa"),
        "-src-file",
        os.path.join(str(tmp_path), "source.src"),
        "-o",
        os.path.join(str(tmp_path), "sub1.sqn"),
        "-w",
        os.path.join(str(tmp_path), "comment.cmt"),
        "-f",
        os.path.join(str(tmp_path), "sub1.gff"),
    ]
    assert runs == [
        {
            "command": expected_command,
            "stdout": genbank_handler_module.subprocess.PIPE,
            "stderr": genbank_handler_module.subprocess.PIPE,
            "cwd": os.path.join(os.path.dirname(os.path.abspath(genbank_handler_module.__file__))),
        }
    ]

def test_create_table2asn__does_not_include_gff_when_only_cwd_has_matching_file(tmp_path, monkeypatch, genbank_handler_module):
    cwd_gff = Path("sub1.gff")
    cwd_gff.write_text("wrong", encoding="utf-8")
    (tmp_path / "sub1.val").write_text("valid\n", encoding="utf-8")

    def fake_isfile(path):
        if path == "/tmp/table2asn":
            return True
        if path == os.path.join(str(tmp_path), "comment.cmt"):
            return False
        if path == os.path.join(str(tmp_path), "sub1.gff"):
            return False
        if path == os.path.join(str(tmp_path), "sub1.val"):
            return True
        return os.path.isfile(path)

    runs: list[list[str]] = []

    def fake_run(command, stdout, stderr, cwd):
        runs.append(command)
        return types.SimpleNamespace(returncode=0, stdout=b"", stderr=b"")

    def fake_check_table2asn_submission(validation_file):
        return "VALIDATED"

    monkeypatch.setattr(genbank_handler_module.os.path, "isfile", fake_isfile)
    monkeypatch.setattr(genbank_handler_module.subprocess, "run", fake_run)
    monkeypatch.setattr(genbank_handler_module, "check_table2asn_submission", fake_check_table2asn_submission)

    try:
        assert genbank_handler_module.create_table2asn("sub1", str(tmp_path)) == "VALIDATED"
        assert "-f" not in runs[0]
    finally:
        cwd_gff.unlink(missing_ok=True)

def test_create_table2asn__subprocess_failure_exits(tmp_path, monkeypatch, genbank_handler_module, capsys):
    def fake_isfile(path):
        return True

    monkeypatch.setattr(genbank_handler_module.os.path, "isfile", fake_isfile)

    def fake_run(command, stdout, stderr, cwd):
        return types.SimpleNamespace(returncode=1, stdout=b"bad stdout", stderr=b"bad stderr")

    monkeypatch.setattr(genbank_handler_module.subprocess, "run", fake_run)

    with pytest.raises(SystemExit) as exc:
        genbank_handler_module.create_table2asn("sub1", str(tmp_path))

    captured = capsys.readouterr()
    assert exc.value.code == 1
    assert captured.out == (
        "Running Table2asn.\n"
        "b'bad stdout'\n"
    )
    assert captured.err == (
        "Table2asn-Error\n"
        "b'bad stderr'\n"
    )

#*******************************************************************************
#                      check_table2asn_submission
#*******************************************************************************

@pytest.mark.parametrize(("content", "expected"), [(None, "ERROR"), ("Error: bad\n", "ERROR"), ("No problems\n", "VALIDATED"), ("", "ERROR")])
def test_check_table2asn_submission__validation_file_states(tmp_path, genbank_handler_module, content, expected):
    validation_file = tmp_path / "sub1.val"
    if content is not None:
        validation_file.write_text(content)

    assert genbank_handler_module.check_table2asn_submission(str(validation_file)) == expected

def test_check_table2asn_submission__opens_validation_file_in_read_mode(tmp_path, monkeypatch, genbank_handler_module):
    validation_file = tmp_path / "sub1.val"
    validation_file.write_text("No problems\n", encoding="utf-8")

    original_open = open
    observed: list[tuple[str, str]] = []

    def fake_open(file, mode="r", *args, **kwargs):
        observed.append((str(file), mode))
        return original_open(file, mode, *args, **kwargs)

    monkeypatch.setattr("builtins.open", fake_open)

    assert genbank_handler_module.check_table2asn_submission(str(validation_file)) == "VALIDATED"
    assert observed == [(str(validation_file), "r")]

def test_check_table2asn_submission__error_file_prints_exact_stderr(tmp_path, genbank_handler_module, capsys):
    validation_file = tmp_path / "sub1.val"
    validation_file.write_text("Error: bad\n", encoding="utf-8")

    result = genbank_handler_module.check_table2asn_submission(str(validation_file))

    assert result == "ERROR"
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == (
        "Submission has errors after running Table2asn.\n"
        "Resolve issues labeled \"Error:\" in table2asn validation file.\n"
        f"Validation file: {validation_file}\n"
    )

#*******************************************************************************
#                         create_genbank_submission
#*******************************************************************************

def test_create_genbank_submission__forwards_exact_create_files_arguments(
    tmp_path, monkeypatch, genbank_handler_module, ncbi_config, genbank_metadata):
    observed: list[dict[str, Any]] = []

    def fake_create_files(**kwargs):
        observed.append(kwargs)

    def fake_create_zip(**kwargs):
        observed.append({"function": "zip", **kwargs})

    def fake_create_submission_xml(**kwargs):
        observed.append({"function": "xml", **kwargs})
        return b"<Submission />"

    monkeypatch.setattr(genbank_handler_module, "create_files", fake_create_files)
    monkeypatch.setattr(genbank_handler_module, "create_zip", fake_create_zip)
    monkeypatch.setattr(genbank_handler_module, "create_submission_xml", fake_create_submission_xml)

    genbank_handler_module.create_genbank_submission(
        organism="FLU",
        submission_name="sub1",
        submission_dir=str(tmp_path),
        config_dict=ncbi_config,
        metadata=genbank_metadata,
        gff_file="/tmp/input.gff",
        table2asn=False,
        publication_title="Exact publication title",
        publication_status="Published",
    )
    assert observed[0] == {
        "organism": "FLU",
        "submission_name": "sub1",
        "submission_dir": str(tmp_path),
        "config_dict": ncbi_config,
        "metadata": genbank_metadata,
        "gff_file": "/tmp/input.gff",
        "publication_title": "Exact publication title",
        "publication_status": "Published",
    }

@pytest.mark.parametrize("organism", ["FLU", "COV", "OTHER"])
@pytest.mark.parametrize("table2asn", [True, False])
def test_create_genbank_submission__check_ftp_vs_table2asn(
    tmp_path, monkeypatch, genbank_handler_module, ncbi_config, genbank_metadata, organism, table2asn
):
    calls = []

    def fake_create_files(**kwargs):
        calls.append("files")

    def fake_create_table2asn(**kwargs):
        calls.append("table2asn")

    def fake_create_zip(**kwargs):
        calls.append("zip")

    def fake_create_submission_xml(**kwargs):
        calls.append("xml")
        return b"<Submission />"

    monkeypatch.setattr(genbank_handler_module, "create_files", fake_create_files)
    monkeypatch.setattr(genbank_handler_module, "create_table2asn", fake_create_table2asn)
    monkeypatch.setattr(genbank_handler_module, "create_zip", fake_create_zip)
    monkeypatch.setattr(genbank_handler_module, "create_submission_xml", fake_create_submission_xml)

    genbank_handler_module.create_genbank_submission(
        organism=organism,
        submission_name="sub1",
        submission_dir=str(tmp_path),
        config_dict=ncbi_config,
        metadata=genbank_metadata,
        gff_file=None,
        table2asn=table2asn,
        publication_title=None,
        publication_status=None,
    )
    if table2asn or organism == "OTHER":
        assert calls == ["files", "table2asn"]
        assert genbank_handler_module.file_handler.saved_xmls == []
    else:
        assert calls == ["files", "zip", "xml"]
        assert genbank_handler_module.file_handler.saved_xmls[0]["submission_xml"] == b"<Submission />"

#*******************************************************************************
#                       accession_report_to_status_report
#*******************************************************************************

def test_accession_report_to_status_report__renames_and_updates_log(genbank_handler_module):
    report = pd.DataFrame(
        [{"Sequence ID": "GB001", "#Accession": "OP123", "Release Date": "2025-01-01"}]
    )

    genbank_handler_module.accession_report_to_status_report("/tmp/sub", report)

    update = genbank_handler_module.upload_log.updated[0]
    assert update["submission_dir"] == "/tmp/sub"
    assert update["update_database"] == "GENBANK"
    df = update["update_df"]
    assert df.to_dict("records") == [
        {
            "gb-sample_name": "GB001",
            "genbank_status": "PROCESSED",
            "genbank_accession": "OP123",
            "genbank_message": "2025-01-01",
        }
    ]

#*******************************************************************************
#                         process_genbank_report
#*******************************************************************************

def test_process_genbank_report__downloads_to_submission_dir_and_waits_exactly(tmp_path, monkeypatch, genbank_handler_module):
    def fake_process_report_header(report_file):
        return (
            {
                "SubmissionStatus": {
                    "Action": {
                        "@status": "processed-ok",
                        "Response": [
                            {
                                "File": [
                                    {"@file_path": "AccessionReport.tsv", "@file_id": "FILE123"},
                                ]
                            }
                        ],
                    }
                }
            },
            "PROCESSED",
            "SUB123",
        )

    genbank_handler_module.ncbi_handler.process_report_header = fake_process_report_header
    requests_calls: list[tuple[str, bool]] = []
    exists_calls: list[str] = []
    sleep_calls: list[int] = []
    load_calls: list[dict[str, Any]] = []

    def fake_get(url, allow_redirects):
        requests_calls.append((url, allow_redirects))
        return types.SimpleNamespace(content=b"Sequence ID\t#Accession\tRelease Date\nGB001\tOP123\t2025-01-01\n")

    original_exists = genbank_handler_module.os.path.exists

    def fake_exists(path):
        exists_calls.append(path)
        return original_exists(path)

    def fake_sleep(seconds):
        sleep_calls.append(seconds)

    def fake_load_csv(file_path, sep=","):
        load_calls.append({"file_path": file_path, "sep": sep})
        return pd.DataFrame(
            [{"Sequence ID": "GB001", "#Accession": "OP123", "Release Date": "2025-01-01"}]
        )

    monkeypatch.setattr(genbank_handler_module.requests, "get", fake_get)
    monkeypatch.setattr(genbank_handler_module.os.path, "exists", fake_exists)
    monkeypatch.setattr(genbank_handler_module.time, "sleep", fake_sleep)
    monkeypatch.setattr(genbank_handler_module.file_handler, "load_csv", fake_load_csv)
    assert genbank_handler_module.process_genbank_report("report.xml", str(tmp_path)) == ("PROCESSED", "SUB123")
    assert requests_calls == [("https://example.test/files/FILE123", True)]
    assert (tmp_path / "AccessionReport.tsv").read_bytes() == b"Sequence ID\t#Accession\tRelease Date\nGB001\tOP123\t2025-01-01\n"
    assert exists_calls == [os.path.join(str(tmp_path), "AccessionReport.tsv")]
    assert sleep_calls == []
    assert load_calls == [{"file_path": os.path.join(str(tmp_path), "AccessionReport.tsv"), "sep": "\t"}]

def test_process_genbank_report__downloads_accession_report_and_updates_status(tmp_path, monkeypatch, genbank_handler_module):
    def fake_process_report_header(report_file):
        return (
            {
                "SubmissionStatus": {
                    "Action": {
                        "@status": "processed-ok",
                        "Response": [
                            {"Other": "ignored"},
                            {"File": [{"@file_path": "AccessionReport.tsv", "@file_id": "FILE123"}]},
                        ],
                    }
                }
            },
            "PROCESSED",
            "SUB123",
        )
    requested = []

    genbank_handler_module.ncbi_handler.process_report_header = fake_process_report_header
    def fake_get(url, allow_redirects):
        requested.append((url, allow_redirects))
        return types.SimpleNamespace(
            content=b"Sequence ID\t#Accession\tRelease Date\nGB001\tOP123\t2025-01-01\n"
        )

    monkeypatch.setattr(genbank_handler_module.requests, "get", fake_get)

    status, submission_id = genbank_handler_module.process_genbank_report("report.xml", str(tmp_path))

    assert (status, submission_id) == ("PROCESSED", "SUB123")
    assert requested == [("https://example.test/files/FILE123", True)]
    assert (tmp_path / "AccessionReport.tsv").exists()
    assert genbank_handler_module.upload_log.updated[0]["update_database"] == "GENBANK"


def test_process_genbank_report__malformed_report_is_safely_ignored(genbank_handler_module):
    def fake_process_report_header(report_file):
        return (
            {"SubmissionStatus": {"Action": {"@status": "processed-ok"}}},
            "PROCESSING",
            "SUB123",
        )

    genbank_handler_module.ncbi_handler.process_report_header = fake_process_report_header
    assert genbank_handler_module.process_genbank_report("report.xml", "/tmp") == (
        "PROCESSING",
        "SUB123",
    )
    assert genbank_handler_module.upload_log.updated == []

#*******************************************************************************
#                         update_genbank_files
#*******************************************************************************

def write_fasta(path: Path, records):
    with open(path, "w") as handle:
        SeqIO.write(records, handle, "fasta")

def test_update_genbank_files__missing_source_file_exits(tmp_path, genbank_handler_module, capsys):
    parent = tmp_path
    submission_dir = parent / "GENBANK"
    submission_dir.mkdir()
    pd.DataFrame([{"gb-sample_name": "GB001"}]).to_csv(parent / "submission_status_report.csv", index=False)

    with pytest.raises(SystemExit) as exc:
        genbank_handler_module.update_genbank_files(
            {"BIOSAMPLE": True, "SRA": False, "GISAID": False},
            "COV",
            str(submission_dir),
            {"Add_Definition_Line_Accessions": True},
        )
    assert exc.value.code == 1
    assert capsys.readouterr().err == f"Error: submission source file does not exist at {submission_dir / 'source.src'}\n"

def test_update_genbank_files__adds_biosample_sra_to_source_and_fasta_definition_line(tmp_path, genbank_handler_module):
    submission_dir = tmp_path / "GENBANK"
    submission_dir.mkdir()
    pd.DataFrame(
        [
            {
                "gb-sample_name": "GB001",
                "biosample_accession": "SAMN1",
                "sra_accession": "SRR1",
            }
        ]
    ).to_csv(tmp_path / "submission_status_report.csv", index=False)
    pd.DataFrame([{"Sequence_ID": "GB001", "Collection_date": "2024-01-01"}]).to_csv(submission_dir / "source.src", sep="\t", index=False)
    write_fasta(submission_dir / "sequence.fsa", [SeqRecord(Seq("ACGT"), id="GB001", description="GB001 [BioSample=OLD] note")])

    genbank_handler_module.update_genbank_files({"BIOSAMPLE": True, "SRA": True, "GISAID": False}, "COV", str(submission_dir), {"Add_Definition_Line_Accessions": True})

    source = pd.read_csv(submission_dir / "source.src", sep="\t", dtype=str)
    assert source.loc[0, "BioSample"] == "SAMN1"
    assert source.loc[0, "SRA"] == "SRR1"
    fasta_text = (submission_dir / "sequence.fsa").read_text()
    assert "[BioSample=SAMN1]" in fasta_text
    assert "[SRA=SRR1]" in fasta_text
    assert "OLD" not in fasta_text

def test_update_genbank_files__uses_exact_left_merges_for_source_and_comment(tmp_path, monkeypatch, genbank_handler_module):
    submission_dir = tmp_path / "GENBANK"
    submission_dir.mkdir()
    pd.DataFrame(
        [
            {
                "gb-sample_name": "GB001",
                "biosample_accession": "SAMN1",
                "sra_accession": "SRR1",
                "gisaid_accession_epi_isl_id": "EPI_ISL_1",
                "gisaid_accession_epi_id": "EPI123",
            }
        ]
    ).to_csv(tmp_path / "submission_status_report.csv", index=False)
    pd.DataFrame([{"Sequence_ID": "GB001"}]).to_csv(submission_dir / "source.src", sep="\t", index=False)
    pd.DataFrame([{"SeqID": "GB001", "StructuredCommentPrefix": "Assembly-Data", "StructuredCommentSuffix": "Assembly-Data"}]).to_csv(submission_dir / "comment.cmt", sep="\t", index=False)
    original_merge = genbank_handler_module.pd.merge
    observed: list[dict[str, Any]] = []

    def fake_merge(left, right, *args, **kwargs):
        observed.append({"left_columns": left.columns.tolist(), "right_columns": right.columns.tolist(), "args": args, "kwargs": kwargs})
        return original_merge(left, right, *args, **kwargs)

    monkeypatch.setattr(genbank_handler_module.pd, "merge", fake_merge)
    genbank_handler_module.update_genbank_files({"BIOSAMPLE": True, "SRA": True, "GISAID": True}, "COV", str(submission_dir), {"Add_Definition_Line_Accessions": False})
    assert observed[0]["kwargs"] == {"how": "left", "on": "Sequence_ID"}
    assert observed[1]["kwargs"] == {"how": "left", "on": "SeqID"}

def test_update_genbank_files__does_not_link_biosample_when_linking_false_even_if_accession_exists(tmp_path, genbank_handler_module):
    submission_dir = tmp_path / "GENBANK"
    submission_dir.mkdir()
    pd.DataFrame([{"gb-sample_name": "GB001", "biosample_accession": "SAMN1"}]).to_csv(tmp_path / "submission_status_report.csv", index=False)
    pd.DataFrame([{"Sequence_ID": "GB001"}]).to_csv(submission_dir / "source.src", sep="\t", index=False)
    genbank_handler_module.update_genbank_files({"BIOSAMPLE": False, "SRA": False, "GISAID": False}, "COV", str(submission_dir), {"Add_Definition_Line_Accessions": False})
    assert genbank_handler_module.file_handler.saved_csvs == []

def test_update_genbank_files__does_not_link_sra_when_linking_false_even_if_accession_exists(tmp_path, genbank_handler_module):
    submission_dir = tmp_path / "GENBANK"
    submission_dir.mkdir()
    pd.DataFrame([{"gb-sample_name": "GB001", "sra_accession": "SRR1"}]).to_csv(tmp_path / "submission_status_report.csv", index=False)
    pd.DataFrame([{"Sequence_ID": "GB001"}]).to_csv(submission_dir / "source.src", sep="\t", index=False)
    genbank_handler_module.update_genbank_files({"BIOSAMPLE": False, "SRA": False, "GISAID": False}, "COV", str(submission_dir), {"Add_Definition_Line_Accessions": False})
    assert genbank_handler_module.file_handler.saved_csvs == []

def test_update_genbank_files__does_not_link_gisaid_when_linking_false_even_if_accessions_exist(tmp_path, genbank_handler_module):
    submission_dir = tmp_path / "GENBANK"
    submission_dir.mkdir()
    pd.DataFrame(
        [
            {
                "gb-sample_name": "GB001",
                "gisaid_accession_epi_isl_id": "EPI_ISL_1",
                "gisaid_accession_epi_id": "EPI123",
            }
        ]
    ).to_csv(tmp_path / "submission_status_report.csv", index=False)
    pd.DataFrame([{"Sequence_ID": "GB001"}]).to_csv(submission_dir / "source.src", sep="\t", index=False)
    genbank_handler_module.update_genbank_files({"BIOSAMPLE": False, "SRA": False, "GISAID": False}, "COV", str(submission_dir), {})
    assert genbank_handler_module.file_handler.saved_csvs == []
    assert not (submission_dir / "comment.cmt").exists()

def test_update_genbank_files__rewrites_fasta_definition_line_exactly(tmp_path, genbank_handler_module):
    submission_dir = tmp_path / "GENBANK"
    submission_dir.mkdir()
    pd.DataFrame([{"gb-sample_name": "GB001","biosample_accession": "SAMN1","sra_accession": "SRR1"}]).to_csv(tmp_path / "submission_status_report.csv", index=False)
    pd.DataFrame([{"Sequence_ID": "GB001"}]).to_csv(submission_dir / "source.src", sep="\t", index=False)
    write_fasta(submission_dir / "sequence.fsa", [SeqRecord(Seq("ACGT"), id="GB001", description="GB001 old desc [BioSample=OLD] [SRA=OLD]")])
    genbank_handler_module.update_genbank_files({"BIOSAMPLE": True, "SRA": True, "GISAID": False}, "COV", str(submission_dir), {"Add_Definition_Line_Accessions": True})
    records = list(SeqIO.parse(submission_dir / "sequence.fsa", "fasta"))
    assert len(records) == 1
    assert records[0].id == "GB001"
    assert records[0].description == "GB001 old desc [BioSample=SAMN1] [SRA=SRR1]"

def test_update_genbank_files__does_not_save_when_no_linking_accessions_available(tmp_path, genbank_handler_module):
    submission_dir = tmp_path / "GENBANK"
    submission_dir.mkdir()
    pd.DataFrame([{"gb-sample_name": "GB001"}]).to_csv(tmp_path / "submission_status_report.csv", index=False)
    pd.DataFrame([{"Sequence_ID": "GB001"}]).to_csv(submission_dir / "source.src", sep="\t", index=False)
    before_saved = list(genbank_handler_module.file_handler.saved_csvs)
    genbank_handler_module.update_genbank_files({"BIOSAMPLE": False, "SRA": False, "GISAID": False}, "COV", str(submission_dir), {"Add_Definition_Line_Accessions": True})
    assert genbank_handler_module.file_handler.saved_csvs == before_saved

def test_update_genbank_files__does_not_rewrite_fasta_when_config_disabled( tmp_path, genbank_handler_module):
    submission_dir = tmp_path / "GENBANK"
    submission_dir.mkdir()
    pd.DataFrame([{"gb-sample_name": "GB001", "biosample_accession": "SAMN1"}]).to_csv(tmp_path / "submission_status_report.csv", index=False)
    pd.DataFrame([{"Sequence_ID": "GB001"}]).to_csv(submission_dir / "source.src", sep="\t", index=False)
    fasta = submission_dir / "sequence.fsa"
    fasta.write_text(">GB001 old description\nACGT\n")
    genbank_handler_module.update_genbank_files({"BIOSAMPLE": True, "SRA": False, "GISAID": False}, "COV", str(submission_dir), {"Add_Definition_Line_Accessions": False})
    assert fasta.read_text() == ">GB001 old description\nACGT\n"

def test_update_genbank_files__adds_gisaid_accessions_to_existing_comment_file(tmp_path, genbank_handler_module):
    submission_dir = tmp_path / "GENBANK"
    submission_dir.mkdir()
    pd.DataFrame(
        [
            {
                "gb-sample_name": "GB001",
                "gisaid_accession_epi_isl_id": "EPI_ISL_1",
                "gisaid_accession_epi_id": "EPI123",
            }
        ]
    ).to_csv(tmp_path / "submission_status_report.csv", index=False)
    pd.DataFrame([{"Sequence_ID": "GB001"}]).to_csv(submission_dir / "source.src", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "SeqID": "GB001",
                "StructuredCommentPrefix": "Assembly-Data",
                "Assembly Method": "iVar",
                "StructuredCommentSuffix": "Assembly-Data",
            }
        ]
    ).to_csv(submission_dir / "comment.cmt", sep="\t", index=False)
    genbank_handler_module.update_genbank_files({"BIOSAMPLE": False, "SRA": False, "GISAID": True}, "COV", str(submission_dir), {})
    cmt = pd.read_csv(submission_dir / "comment.cmt", sep="\t", dtype=str)
    assert cmt.columns.tolist() == ["SeqID", "StructuredCommentPrefix", "EPI_ISOLATE_ID", "Assembly Method", "EPI_SEQUENCE_ID", "StructuredCommentSuffix"]
    assert cmt.loc[0, "EPI_ISOLATE_ID"] == "EPI_ISL_1"
    assert cmt.loc[0, "EPI_SEQUENCE_ID"] == "EPI123"

@pytest.mark.parametrize(("organism", "expected_prefix"), [("FLU", "FluData"), ("COV", "Assembly-Data")])
def test_update_genbank_files__creates_comment_file_for_gisaid_when_missing(tmp_path, genbank_handler_module, organism, expected_prefix):
    submission_dir = tmp_path / "GENBANK"
    submission_dir.mkdir()
    pd.DataFrame(
        [
            {
                "gb-sample_name": "GB001",
                "gisaid_accession_epi_isl_id": "EPI_ISL_1",
                "gisaid_accession_epi_id": "EPI123",
            }
        ]
    ).to_csv(tmp_path / "submission_status_report.csv", index=False)
    pd.DataFrame([{"Sequence_ID": "GB001"}]).to_csv(submission_dir / "source.src", sep="\t", index=False)
    genbank_handler_module.update_genbank_files({"BIOSAMPLE": False, "SRA": False, "GISAID": True}, organism, str(submission_dir), {})
    cmt = pd.read_csv(submission_dir / "comment.cmt", sep="\t", dtype=str)
    assert cmt.loc[0, "StructuredCommentPrefix"] == expected_prefix
    assert cmt.loc[0, "StructuredCommentSuffix"] == expected_prefix
    assert cmt.loc[0, "EPI_ISOLATE_ID"] == "EPI_ISL_1"
    assert cmt.loc[0, "EPI_SEQUENCE_ID"] == "EPI123"

def test_update_genbank_files__fasta_write_permission_error_exits(tmp_path, monkeypatch, genbank_handler_module, capsys):
    submission_dir = tmp_path / "GENBANK"
    submission_dir.mkdir()
    pd.DataFrame([{"gb-sample_name": "GB001", "biosample_accession": "SAMN1"}]).to_csv(tmp_path / "submission_status_report.csv", index=False)
    pd.DataFrame([{"Sequence_ID": "GB001"}]).to_csv(submission_dir / "source.src", sep="\t", index=False)
    write_fasta(submission_dir / "sequence.fsa", [SeqRecord(Seq("ACGT"), id="GB001", description="GB001")])
    real_open = open

    def fake_open(path, mode="r", *args, **kwargs):
        if Path(path).name == "sequence.fsa" and "w" in mode:
            raise PermissionError("PermissionError")
        return real_open(path, mode, *args, **kwargs)

    monkeypatch.setattr(genbank_handler_module, "open", fake_open, raising=False)
    with pytest.raises(SystemExit) as exc:
        genbank_handler_module.update_genbank_files({"BIOSAMPLE": True, "SRA": False, "GISAID": False}, "COV", str(submission_dir), {"Add_Definition_Line_Accessions": True})
    assert exc.value.code == 1
    captured = capsys.readouterr()
    assert captured.err == f"Error: Permission error when trying to save 'sequence.fsa' to path: {submission_dir}\nPermissionError\n"

def test_update_genbank_files__loads_existing_comment_file_with_exact_tsv_path(tmp_path, monkeypatch, genbank_handler_module):
    submission_dir = tmp_path / "GENBANK"
    submission_dir.mkdir()
    pd.DataFrame([{"gb-sample_name": "GB001", "gisaid_accession_epi_isl_id": "EPI_ISL_1"}]).to_csv(tmp_path / "submission_status_report.csv", index=False)
    pd.DataFrame([{"Sequence_ID": "GB001"}]).to_csv(submission_dir / "source.src", sep="\t", index=False)
    pd.DataFrame([{"SeqID": "GB001", "StructuredCommentPrefix": "Assembly-Data", "StructuredCommentSuffix": "Assembly-Data"}]).to_csv(submission_dir / "comment.cmt", sep="\t", index=False)
    original_load_csv = genbank_handler_module.file_handler.load_csv
    load_calls: list[dict[str, Any]] = []

    def fake_load_csv(file_path, sep=","):
        load_calls.append({"file_path": file_path, "sep": sep})
        return original_load_csv(file_path=file_path, sep=sep)

    monkeypatch.setattr(genbank_handler_module.file_handler, "load_csv", fake_load_csv)
    genbank_handler_module.update_genbank_files({"BIOSAMPLE": False, "SRA": False, "GISAID": True}, "COV", str(submission_dir), {})
    assert load_calls == [
        {"file_path": os.path.join(str(tmp_path), "submission_status_report.csv"), "sep": ","},
        {"file_path": os.path.join(str(submission_dir), "source.src"), "sep": "\t"},
        {"file_path": os.path.join(str(submission_dir), "comment.cmt"), "sep": "\t"},
        {"file_path": os.path.join(str(submission_dir), "comment.cmt"), "sep": "\t"},
    ]

def test_update_genbank_files__fasta_definition_omits_missing_or_blank_accessions(tmp_path, genbank_handler_module):
    submission_dir = tmp_path / "GENBANK"
    submission_dir.mkdir()
    pd.DataFrame([{"gb-sample_name": "GB001", "biosample_accession": "", "sra_accession": ""}]).to_csv(tmp_path / "submission_status_report.csv", index=False)
    pd.DataFrame([{"Sequence_ID": "GB001"}]).to_csv(submission_dir / "source.src", sep="\t", index=False)
    write_fasta(submission_dir / "sequence.fsa", [SeqRecord(Seq("ACGT"), id="GB001", description="GB001 original")])
    genbank_handler_module.update_genbank_files({"BIOSAMPLE": True, "SRA": True, "GISAID": False}, "COV", str(submission_dir), {"Add_Definition_Line_Accessions": True})
    records = list(SeqIO.parse(submission_dir / "sequence.fsa", "fasta"))
    assert records[0].description == "GB001 original"
    assert "None" not in records[0].description
    assert "XXXX" not in records[0].description
    assert "[BioSample=" not in records[0].description
    assert "[SRA=" not in records[0].description

def test_update_genbank_files__fasta_write_unexpected_error_exits_exactly(tmp_path, monkeypatch, genbank_handler_module, capsys):
    submission_dir = tmp_path / "GENBANK"
    submission_dir.mkdir()
    pd.DataFrame([{"gb-sample_name": "GB001", "biosample_accession": "SAMN1"}]).to_csv(tmp_path / "submission_status_report.csv", index=False)
    pd.DataFrame([{"Sequence_ID": "GB001"}]).to_csv(submission_dir / "source.src", sep="\t", index=False)
    write_fasta(submission_dir / "sequence.fsa", [SeqRecord(Seq("ACGT"), id="GB001", description="GB001")])
    real_open = open

    def fake_open(path, mode="r", *args, **kwargs):
        if Path(path).name == "sequence.fsa" and "w" in mode:
            raise RuntimeError("unexpected")
        return real_open(path, mode, *args, **kwargs)

    monkeypatch.setattr(genbank_handler_module, "open", fake_open, raising=False)
    with pytest.raises(SystemExit) as exc:
        genbank_handler_module.update_genbank_files({"BIOSAMPLE": True, "SRA": False, "GISAID": False}, "COV", str(submission_dir), {"Add_Definition_Line_Accessions": True})

    assert exc.value.code == 1
    captured = capsys.readouterr()
    assert captured.err == f"Error: An unexpected error occurred when trying to save 'sequence.fsa' to path: {submission_dir}\nunexpected\n"
