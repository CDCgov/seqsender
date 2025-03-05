
# Python Libraries
from Bio import Entrez
from settings import BATCH_SIZE, DEFAULT_RATE_LIMIT, API_RATE_LIMIT
import time

def load_entrez_credentials(config_dict: Dict[str, Any]) -> None:
    Entrez.email = config_dict[]
    Entrez.api_key = config_dict[]

# Return time value in seconds
def calculate_entrez_rate_limit(config_dict: Dict[str, Any]):
    rate_limit = DEFAULT_RATE_LIMIT
    if load_entrez_api_key(config_dict=config_dict):
        rate_limit = API_RATE_LIMIT
    return 1 / rate_limit

def create_entrez_query(database: str, query_list: List[str]) -> Tuple[str]:
    web_env = ""
    query_key = ""
    with Entrez.epost(db=database, id=query_list) as connection:
        epost_info = Entrez.read(connection)
    if "WebEnv" in epost_info and "QueryKey" in epost_info:
        web_env = epost_info["WebEnv"]
        query_key = epost_info["QueryKey"]
    else:
        logger.critical(f"Unable to connect to Entrez:\n{epost_info}")
        sys.exit(1)
    return web_env, query_key

def search_entrez(database: str, submission_dir: str, config_dict: Dict[str, Any], query_type: str, query_list: List[str]):
    web_env, query_key = create_entrez_query(query_list=query_list)
    xml_document = []
    # Limit requests
    delay = calculate_entrez_rate_limit(config_dict=config_dict)
    for index in range(0, len(query_list), BATCH_SIZE):
        with Entrez.efetch(db=database, query_key=query_key, WebEnv=web_env, rettype="xml", retstart=index, retmax=BATCH_SIZE) as connection:
            xml_document.append(connection.read())
    for index, entrez_query in enumerate(batch_query):
        connection = Entrez.efetch(db=database, id = entrez_query, rettype = "xml")
        xml_string = xmltodict(connection.read().decode("utf-8"))
        time.sleep(delay)

def search_entrez_for_linking_databases():
