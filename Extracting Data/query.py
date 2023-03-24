import requests
# from .credentials import *
from functions import *
import pandas as pd
from io import StringIO
print(__package__)
from os.path import join
from functools import reduce

class Query:
    def __init__(self, url, id_or_user=None, pass_or_secret=None, auth_type:str=None, is_fdq=False):
        self.url = url
        if not is_fdq:
            self.id_or_user = id_or_user
            self.pass_or_secret = pass_or_secret
        self.is_fdq = is_fdq
        self.auth_type = auth_type
        self.__connect_to()

    def __set_headers(self, accept_type='text/csv'):
        if self.auth_type:
            headers = {
                'Content-Type': 'application/sparql-query',
                'Accept': accept_type,
                'Authorization': self.auth}  # text/html
        else:
            headers = {'Content-Type': 'application/x-www-form-urlencoded'}
        return headers

    @auth_oauth
    def __set_req_params_os2(self, accept_type='text/csv'):
        self.url = self.url + "/dataplatform/proxy/default/sparql"
        self.headers = self.__set_headers(accept_type)

    @auth_basic
    def __set_req_params_basic(self, accept_type='text/csv'):
        self.headers = self.__set_headers(accept_type)

    def __set_req_params(self, accept_type='text/json'):
        self.headers = self.__set_headers(accept_type)

    def __connect_to(self, accept_type='text/csv'):
        if self.auth_type=='oauth':
            self.__set_req_params_os2()
        elif self.auth_type=='basic':
            self.__set_req_params_basic()
        elif self.is_fdq:
            self.__set_req_params()
        else:
            raise ValueError('wrong parameters or auth_type (oauth, basic ) is missing !!')

    @timer
    def post_query(self, query, is_service=False):
        try:
            if self.is_fdq and is_service:
                response = requests.request("POST", self.url, headers=self.headers, data="query="+query+"&sparql1_1=True")
            elif self.is_fdq:
                response = requests.request("POST", self.url, headers=self.headers, data="""query="""+query)
            else:
                response = requests.request("POST", self.url, headers=self.headers, data=query, stream=True, timeout=3600)

            if response.status_code == 200:
                print("Passed Query Description: {}".format(response.status_code))
                self.response = response
            else:
                print("Failed Query: {}".format(response.content))
        except Exception as e:
            print("Failed Query: {}".format(str(e)))

    @timer
    def to_save(self, save_path, filename:str):
        filename = '_'.join(filename.split(' '))
        if self.response.headers['Content-Type']=='text/csv; charset=utf-8':
            return pd.read_csv(StringIO(str(self.response.content, 'utf-8')))
            #with open(join(save_path, filename+'.csv'), 'wb') as fd:
                #for chunk in self.response.iter_content(chunk_size=128):
                    #fd.write(chunk)
        elif self.response.headers['Content-Type']=='text/turtle':
            with open(join(save_path, filename+'.ttl'), 'wb') as fd:
                for chunk in self.response.iter_content(chunk_size=1024):
                    fd.write(chunk)
        elif self.response.headers['Content-Type']=='application/json':
            with open(join(save_path, filename+'.json'), 'wb') as fd:
                for chunk in self.response.iter_content(chunk_size=128):
                    fd.write(chunk)
            json_to_csv(self.response.json()['results']['bindings'], save_path, \
                filename, columns=[var+'.value' for var in self.response.json()['head']['vars']])
        else:
            raise TypeError("Unknown response content-type")

        del self.response
        #return pd.read_csv(StringIO(str(self.response.content, 'utf-8')))


def main(client_url='', client_id='', client_secret='',
         query="""SELECT DISTINCT ?s ?o WHERE{?s a ?o.} LIMIT 10"""):

    cmemc_query = Query(client_url, client_id, client_secret)
    print(cmemc_query.get_response(query))


#    fuseki_query = Query(fuseki_endpoint, fuseki_user_infai, fuseki_pw_infai)
#    print (fuseki_query.get_response(query))


if __name__ == "__main__":

    query = """PREFIX  rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        SELECT DISTINCT ?s ?o WHERE{?s a ?o.} LIMIT 10"""

    main(client_url_tib, client_id_tib,
         client_secret_tib, query)
    main (client_url_imp, client_id_imp,
         client_secret_imp, query)
    main(client_url_infai, client_id_infai,
         client_secret_infai, query)
