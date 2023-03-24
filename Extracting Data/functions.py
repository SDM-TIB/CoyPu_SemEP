from random import random
import pandas as pd
import requests
import base64
import time
import functools
from os.path import join
from rdflib import Graph
from rdflib.plugins.sparql.processor import SPARQLResult
from SPARQLWrapper import SPARQLWrapper, JSON


def auth_oauth(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        url = args[0].url + "/auth/realms/cmem/protocol/openid-connect/token"
        payload = 'grant_type=client_credentials&client_id={}&client_secret={}'\
            .format(args[0].id_or_user, args[0].pass_or_secret)

        headers = {
                    'Content-Type': 'application/x-www-form-urlencoded'
                }

        response = requests.request("POST", url, headers=headers, data=payload)

        if response.status_code == 200:
            args[0].auth = 'Bearer ' + response.json()['access_token']
        
        return func(*args, **kwargs)
    return wrapper


def auth_basic(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        usr_pass = args[0].id_or_user + ':' + args[0].pass_or_secret
        args[0].auth =  "Basic {}".format(base64.b64encode(usr_pass.encode()).decode())
        return func(*args, **kwargs)
    return wrapper


def timer(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        time_pre = time.time()
        ret = func(*args, **kwargs)
        print('Total time taken by {} function: {}'.format(func.__name__,time.time()-time_pre))
        return ret
    return wrapper


@timer
def json_to_csv(json, save_path, filename, columns=None, ret=False):
    df = pd.json_normalize(json)
    # print(df.columns)
    print(df.memory_usage())
    print(df.info(memory_usage=True))
    df.to_csv(join(save_path, filename+'.csv'), encoding='utf-8', columns=columns, index=False)
    if ret:
        return df[columns]


def turtle_to_df(path, f_name):
    g = Graph()
    g.parse(path + f_name + ".ttl", format="ttl")
    print('triples:', len(g))
    g.serialize(destination=path+f_name + ".nt", format="nt", )
    df_g = pd.read_csv(f_name + '.nt', delimiter=" ", header=None)
    df_g.columns = ['s', 'p', 'o', 'dot']
    df_g = df_g[['s', 'p', 'o']]
    return df_g, g


def adding_prefix(df_g):
    df_g.replace('https://data.coypu.org/event/', 'event:', regex=True, inplace=True)
    df_g.replace('https://schema.coypu.org/global#', 'coy:', regex=True, inplace=True)
    df_g.replace('http://www.w3.org/1999/02/22-rdf-syntax-ns#', 'type:', regex=True, inplace=True)
    df_g.replace('https://data.coypu.org/country/', 'country:', regex=True, inplace=True)
    df_g.replace('http://www.w3.org/2001/XMLSchema#', 'xsd:', regex=True, inplace=True)
    df_g.replace('https://data.coypu.org/icews/', 'SPEvent:', regex=True, inplace=True)
    df_g.replace('https://data.coypu.org/event/acled/', 'conflict:', regex=True, inplace=True)
    df_g.replace('https://data.coypu.org/acled/', 'acled:', regex=True, inplace=True)
    df_g.replace('>', '', regex=True, inplace=True)
    df_g.replace('<', '', regex=True, inplace=True)
    return df_g


def filter_by_year(icews, year):
    df = icews.loc[icews.p == 'coy:hasDate']
    df['o'] = pd.to_datetime(df['o'])
    df = df.loc[(df.o < pd.Timestamp(year - 1, 12, 31)) | (df.o > pd.Timestamp(year, 12, 31))]

    icews = icews.loc[~icews.s.isin(df.s.tolist())]
    return icews


def filter_by_year_disaster(g, year, df_g):
    query = """PREFIX coy: <https://schema.coypu.org/global#>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
    select distinct ?s
    where {
        ?s coy:hasEventStatus ?EventStatus .
        ?s coy:hasStartTimestamp ?StartTimestamp .
        FILTER(?EventStatus in ("Ongoing"@en, "Current"@en, "Alert"@en)  &&
            (?StartTimestamp > xsd:dateTime(""" + '"' + str(year) + """-12-31T00:00:00Z")))
        }
        """
    qres = g.query(query)
    df = sparql_results_to_df(qres)

    query = """PREFIX coy: <https://schema.coypu.org/global#>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
    select distinct ?s
    where {
        ?s coy:hasEventStatus ?EventStatus .
        ?s coy:hasStartTimestamp ?StartTimestamp .
        ?s coy:hasEndTimestamp ?EndTimestamp .
        FILTER((?StartTimestamp > xsd:dateTime(""" + '"' + str(
        year) + """-12-31T00:00:00Z") && (?EndTimestamp > xsd:dateTime(""" + '"' + str(year) + """-12-31T00:00:00Z")))
                || (?EndTimestamp < xsd:dateTime(""" + '"' + str(year - 1) + """-12-31T00:00:00Z")))
        FILTER(?EventStatus="Past"@en)
        }
        """
    qres = g.query(query)
    df_1 = sparql_results_to_df(qres)
    df = pd.concat([df, df_1])
    df['s'] = '<' + df['s'].astype(str) + '>'
    df_g = df_g.loc[~df_g.s.isin(df.s.tolist())]
    return df_g


def filter_by_year_conflict(fuseki_infai, path, year, df_g):
    query = """PREFIX coy: <https://schema.coypu.org/global#>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
    select distinct ?s
    where {
        ?s a coy:Conflict .
        ?s coy:hasTimestamp ?Timestamp .
        FILTER((?Timestamp > xsd:date(""" + '"' + str(year) + """-12-31")) ||
                (?Timestamp < xsd:date(""" + '"' + str(year - 1) + """-12-31")))
        }
        """
    fuseki_infai.post_query(query)
    df = fuseki_infai.to_save(path, filename='conflict_year')
    df['s'] = '<' + df['s'].astype(str) + '>'
    df_g = df_g.loc[~df_g.s.isin(df.s.tolist())]
    return df_g


def sparql_results_to_df(results: SPARQLResult) -> pd.DataFrame:
    """
    Export results from an rdflib SPARQL query into a `pandas.DataFrame`,
    using Python types. See https://github.com/RDFLib/rdflib/issues/1179.
    """
    return pd.DataFrame(
        data=([None if x is None else x.toPython() for x in row] for row in results),
        columns=[str(x) for x in results.vars],
    )


def read_result_country(query, sparql):
    df = []
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    for r in results['results']['bindings']:
        row = {'s': r['s']['value'].replace(',', ''),
               'o': r['o']['value']}
        df.append(row)

    return pd.DataFrame.from_dict(df)
