#!/usr/bin/env python3

"""
    Upload consensus sequences and metadata to GISAID's EpiCoV
    Copyright (C) 2022 Freunde von GISAID e.V.
"""


import argparse, getpass, sys, requests, secrets, hashlib, json, time, csv, os
from itertools import islice
import pandas as pd
from cli3.utils import __st__
from urllib3.exceptions import InsecureRequestWarning

# Suppress only the single warning from urllib3 needed.
requests.packages.urllib3.disable_warnings(category=InsecureRequestWarning)

# global configuration
proxy = None
GPS_SERVICE_URL_TEST = "https://gpsapi-test.epicov.org/epi3/gps_api"
GPS_SERVICE_URL_LOCAL = "http://localhost:8080/epi3/gps_api"
GPS_SERVICE_URL = "https://gpsapi.epicov.org/epi3/gps_api"
TAB="\t"
FRAMESHIFTS = {"catch_none": "all", # all frameshifts shall pass
               "catch_all": "none", # no frameshifts shall pass
               "catch_novel": "com"} # only common frameshifts shall pass

def get_service_url(client_type):
    """Get service url."""
    if client_type == "test":
        return GPS_SERVICE_URL_TEST
    elif client_type == "local":
        return GPS_SERVICE_URL_LOCAL
    else:
        return GPS_SERVICE_URL


def call_api(params, debug, client_type):
    """ call the GPS-API """
    global proxy
    body = json.dumps(params)
    try:
        r = requests.post(get_service_url(client_type),
                          data = body,
                          proxies = proxy,
                          verify = False)
        time.sleep(0.1)
        return r.json()
    except requests.exceptions.ProxyError:
        if debug:
            raise
        else:
            return {"rc": "proxy_error"}
    except KeyboardInterrupt:
        raise
    except:
        if debug:
            raise
        else:
            return {"rc": "connection_error"}

def sha512_hexdigest(inp):
    """Create sha512"""
    hasher = hashlib.sha512()
    hasher.update(inp)
    return hasher.hexdigest()

def request_authentication_token(ctx, client_id, user, password, debug, client_type):
    """
    Request authentication token from GISAID GPS-API.
    Used only for the `authenticate` sub-command (since >=v2.0.1)
    """
    password = password.encode("utf8")
    salt = secrets.token_urlsafe(64)
    hashed = salt + "/" + sha512_hexdigest((salt + sha512_hexdigest(password)).encode("ascii"))

    # resp is response
    resp = call_api({"cmd": "state/auth/get_token",
                     "api": {"version": 1},
                     "ctx": ctx,
                     "client_id": client_id,
                     "login": user,
                     "hash": hashed}, debug,
                     client_type)

    return resp["rc"], resp.get("auth_token"), resp.get("valid_until")


def write_authentication_token_file(ctx, filename, client_id, token, valid_until): #deleted unused 'debug' from args
    """
    Write the authentication token.
    """
    # It seems unusual behaviour that 'tokens' variable is available outside scope of
    # try-except block.  This issue is discussed here: https://bugs.python.org/issue40728
    # and is done because python needs to clean the exception variable outside the block
    # to avoid reference cycles.

    tokens = {}
    tokens[ctx] = {'client-id': client_id,
                   'token': token,
                   'expiry': valid_until,
                  }
    with open(filename, "w") as file:
        file.write(json.dumps(tokens))
    return tokens

def load_authentication_token(ctx, filename, debug):
    """Load the authentication token."""
    try:
        with open(filename, "r") as file:
            return json.loads(file.read())[ctx] # get the row out of json *.authtoken file
    except:
        if debug:
            raise
        else:
            return None

def parse_csv(filename):
    """Parse a csv file."""
    df = None
    with open(filename, "r") as input_handle:
        import pandas as pd
        df = pd.read_csv(input_handle)
        print(df)
    file = open(filename, "r", encoding = "utf-8-sig")
    return csv.DictReader(file, delimiter = ",", quotechar = "\"")

def parse_fasta(fn):
    """
    >>> import tempfile
    >>> with tempfile.NamedTemporaryFile() as tmp:
    ...    tmp.write(b'''>EPI_ISL_8747624786432768.1
    ...    ACGACCTAACTGAGGA
    ...    ACACCTAACTGAGAGA''')
    ...    tmp.seek(0)
    ...    sequences = {header: sequence for header, sequence in parse_fasta(tmp.name)}
    ...    print(sequences)
    67
    0
    {'EPI_ISL_8747624786432768.1': 'ACGACCTAACTGAGGAACACCTAACTGAGAGA'}

    """


    seq = None
    header = None
    for row in open(fn, "r", encoding = "utf-8-sig"):
        row = row.strip()
        if row:
            if row[0] == ">":
                if seq:
                    yield header, "".join(seq)
                header = row[1:].lstrip()
                seq = []
            else:
                if ">" in row:
                    raise ValueError("Found '>' in sequence-data")
                seq.append(row)
    if seq:
        yield header, "".join(seq)


def split_every(n, iterable):
    """Split iterable into n-sized chunks."""
    iterable = iter(iterable)
    yield from iter(lambda: list(islice(iterable, n)), [])


def handle(args):
    """
    Handle execution of the workflow:
        authenticate
        upload
    """
    if args.subparser_name == "upload" and args.proxy:
        proxy = {"https": args.proxy,
                 "http": args.proxy}
    client_type = 'live'

    logfile = []; log_dict = [];
    def log(code, *msg):
        log_dict.append({"code": code, "msg": " ".join(msg)})
        print(f"{code}: {', '.join(msg)}")

    failed_writer = None #can avoid this setup in later iterations of cli3 software
    failed_file = None
    def output_faulty(submission):
        if failed_writer:
            failed_writer.writerow(submission)

    # execute CLI-commands
    if args.subparser_name == "authenticate":
        from pathlib import PurePath, Path
        if args.force and Path(args.token).exists():
            print(f"Will overwrite '{Path(args.token)}'.")
        elif (args.force and not Path(args.token).exists()) or (not args.force and not Path(args.token).exists()):
            print(f"Will write authentication token to '{Path(args.token)}'")
        else:
            sys.exit(f"'{args.token}' already exists and force overwrite is off, so doing nothing.")
        if not args.username:
            args.username = getpass.getpass("Enter username: ")
        if not args.password:
            args.password = getpass.getpass("Enter password: ")
        if not args.client_id:
            args.client_id = getpass.getpass("Enter client-ID: ")

        if args.client_id.startswith("TEST-"):
            client_type = "test"
        elif args.client_id.startswith("LOCAL-"):
            client_type = "local"

        rc, auth_token, valid_until = request_authentication_token(ctx = args.database,
                                                                   client_id = args.client_id,
                                                                   user = args.username,
                                                                   password = args.password,
                                                                   debug = True,
                                                                   client_type = client_type)
        if rc == "ok":
            logfile.append(log("\nok", "user authenticated"))
            logfile.append(log("valid until", time.ctime(valid_until)))
        else:
            print("Response from server:", rc)

        # save it or report error
        if auth_token:
            write_authentication_token_file(args.database, args.token, args.client_id, auth_token, time.ctime(valid_until))
            return 0, logfile
        else:
            return 1, logfile

    elif args.subparser_name == "upload":

        # get the auth-token from the file or from the command line
        from pathlib import PurePath, Path
        if not Path(args.token).exists():
            sys.exit(f"{args.token} does not exist.  Create file using 'cli3 authenticate'.")

        auth_token = load_authentication_token(args.database, args.token, True)
        client_id = auth_token['client-id']
        auth_token = auth_token['token']

        if client_id.startswith("TEST-"):
            client_type = "test"
        elif client_id.startswith("LOCAL-"):
            client_type = "local"

        metadata = parse_csv(args.metadata)
##        metadata = pd.read_csv(metadata)

        # prepare a CSV-Write with the same header as the input-file for the failed uploads
        if args.failed: # todo: switch this to pandas, *** args.failed is always true
            failed_file = open(args.failed, "a")
            failed_writer = csv.DictWriter(failed_file,
                                           metadata.fieldnames,
                                           delimiter = ",",
                                           quotechar = "\"",
                                           extrasaction = "ignore")
            failed_writer.writeheader()

        # parse FASTA
        try:
            sequences = {header: sequence for header, sequence in parse_fasta(args.fasta)}
        except ValueError:
            log("fasta_format_error", "upload aborted")
            return 1, logfile

        # check for completeness
        submissions = []
        missing_data = False
        for row in metadata: #column fn does nothing to match with sequence filename, metadata can run without it.
            submission = dict(row)

            if "covv_virus_name" in submission:
                virus_name_key = "covv_virus_name"
                sequence_key = "covv_sequence"
            elif "pox_virus_name" in submission:
                virus_name_key = "pox_virus_name"
                sequence_key = "pox_sequence"
            else:
                log("missing_virus_name", f"{submission}")

            if submission[virus_name_key] not in sequences:
                log("missing_seq", f"{submission[virus_name_key]}")
                missing_data = True

            # name in metadata is key, name in fasta is value
            if submission["submitter"] != "Submitter": # a fix to skip demo row if present
                submission[sequence_key] = sequences.get(submission[virus_name_key])
                submission["_subm_confirmed"] = {"frameshift": FRAMESHIFTS[args.frameshift]}
                submission["_st"] = {"st": f"v {__st__}"}
                submissions.append(submission)

        submission_chunks = split_every(500, submissions)
        count = 0
        ok_count = 0
        for submission_chunk in submission_chunks: # to stop overload on curation, blocks of 500.
            # logon = open session
            resp = call_api({"cmd": "state/session/logon", #cmd is logon
                             "api": {"version": 1},
                             "ctx": args.database,
                             "client_id": client_id,
                             "auth_token": auth_token},
                            debug = True,
                            client_type = client_type)

            time.sleep(0.5) # give system time to save session

            if resp["rc"] == "auth_token_invalid":
                log("invalid_token", "invalid auth-token")
                return 1, logfile
            elif resp["rc"] != "ok":
                log("error", "unexpected response from GPS-API", repr(resp))
                return 1, logfile

            sid = resp["sid"]

            # do uploads
            for submission in submission_chunk:

                if "covv_virus_name" in submission:
                    virus_name_key = "covv_virus_name"
                    sequence_key = "covv_sequence"
                elif "pox_virus_name" in submission:
                    virus_name_key = "pox_virus_name"
                    sequence_key = "pox_sequence"
                else:
                    log("missing_virus_name", f"{submission}")

                submitter = submission["submitter"]
                resp = call_api({"cmd": "data/hcov-19/upload", #cmd is upload
                                 "api": {"version": 1},
                                 "sid": sid,
                                 "ctx": args.database,
                                 "data": submission,
                                 "submitter": submitter,
                                 }, debug = True, client_type = client_type)
                if resp["rc"] == "ok":
                    ok_count += 1
                    log("epi_isl_id", f"{submission.get(virus_name_key, '')}; {resp['accession_id']}")
                else:
                    output_faulty(submission)
                    if "validation" in resp:
                        log("validation_error", "%s; %s; %s" % (submission.get(virus_name_key, ""), resp["rc"], json.dumps(resp["validation"])))
                    else:
                        log("upload_error", "%s; %s" % (submission.get(virus_name_key, ""), resp["rc"])) #resp["rc"]
                count += 1

            # close session
            resp = call_api({"cmd": "state/session/logoff",
                             "api": {"version": 1},
                             "ctx": args.database,
                             "sid": sid},
                            debug = True,
                            client_type = client_type)

        log("upload_count", f"submissions uploaded: {ok_count}")
        log("failed_count", f"submissions failed: {count - ok_count}")

        # Save log file
        if args.log:
            with open(args.log, "wt") as f:
                f.write(json.dumps(log_dict, indent = 4))
                f.close

        # close the file for the failed uploads if needed
        if failed_file:
            failed_file.close()

        return 0, log_dict
