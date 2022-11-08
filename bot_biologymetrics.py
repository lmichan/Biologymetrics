#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import Entrez
from Bio import Medline
from datetime import date, timedelta
from posixpath import isfile
import tweepy.client





def get_dates() -> tuple:
    '''
    Obtiene una tupla con las cadenas de la fecha de inicio y la fecha de fin
    de la semana anterior a la semana en la que se encuentra el ultimo lunes.
    Las cadenas de las fechas tienen el formato empleado para hacer la busqueda de Entrez.
    
    Returns
    -------
    tuple
        Las fechas de inicio y de fin de la semana anterior
        a la semana en la que se encuentra el ultimo lunes.
    '''
    
    #Se obtiene un objeto date que contiene la fecha actual
    current_date = date.today()
    
    #Se obtiene un objeto date que contiene la fecha del ultimo lunes
    last_monday_date = current_date - timedelta(days = current_date.weekday())
    
    #Se obtienen dos objetos date que contienen las fechas de inicio y de fin
    #de la semana anterior a la semana en la que se encuentra el ultimo lunes
    last_week_start_date = last_monday_date - timedelta(days=8)
    last_week_end_date = last_monday_date - timedelta(days=2)
    
    #Se obtienen las cadenas del ultimo par de fechas creadas
    #en el formato empleado para hacer la busqueda
    last_week_start_date_string = last_week_start_date.strftime('%Y/%m/%d')
    last_week_end_date_string = last_week_end_date.strftime('%Y/%m/%d')
    
    return last_week_start_date_string, last_week_end_date_string
    




def get_records(from_date: str, to_date: str) -> list:
    '''
    Obtiene la informacion de los registros creados en PubMed
    sobre publicaciones de bibliometria entre las fechas
    especificadas por ``from_date`` y ``to_date``.
    
    Parameters
    ----------
    from_date : str
        La fecha mas distante en la que se obtienen los registros.
    to_date : str
        La fecha mas cercana en la que se obtienen los registros.
    
    Returns
    -------
    list
        La informacion de los registros.
    
    '''
    
    #Si es lunes se borra el contenido de un archivo
    #que contiene los ids de los registros publicados
    if date.today().weekday() == 0:
        open('IDs_de_registros_publicados.txt', 'w').close()
    
    #Se lee de un archivo los ids de los registros que ya
    #se publicaron y estos se guardan en un conjunto
    if isfile('IDs_de_registros_publicados.txt'):
        with open('IDs_de_registros_publicados.txt') as f:
            posted_record_ids = set([x.rstrip() for x in f])
    else:
        posted_record_ids = set()
    
    #Se establece el parametro email de Entrez
    Entrez.email = "A.N.Other@example.com"
    
    #Se define el termino de busqueda
    query = '(altmetr*[Title] OR bibliometr*[Title] OR scientometr*[Title] OR bibliometrics[MeSH Terms]) \
            AND (\"' + from_date + '\"[Date - Create] : \"' + to_date + '\"[Date - Create])'
    
    #Se ejecute una búsqueda de Entrez y se obtiene un identificador a los resultados
    handle = Entrez.esearch(db="pubmed", term=query, retmax=500)
    
    #Se analiza un archivo XML de NCBI Entrez Utilities en objetos de python
    record = Entrez.read(handle)
    
    #Se crea una lista con los ids que se obtienen
    #con la busqueda y que no se han publicado ya
    id_list = [x for x in record["IdList"] if x not in posted_record_ids]
    
    #Se obtienen los resultados de Entrez que son devueltos como un identificador
    handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
    
    #Se leen los registros de Medline desde el identificador y se guardan en una lista
    records = list(Medline.parse(handle))
    
    return records





def tweet_records(records: list):
    '''
    Publica tuits, con la etiqueta #biologymetrics, que contienen los nombres
    de las publicaciones que estan en ``records`` y los enlaces electronicos 
    a sus registros en PubMed.
    
    Parameters
    ----------
    records : list
        La informacion de los registros que se publicaran.
    
    Returns
    -------
    None.
    
    '''
        
    #Se obtiene un agente de usuario utilizado al realizar solicitudes a la API de twitter
    client = tweepy.Client(
        consumer_key = 'consumer_key',
        consumer_secret = 'consumer_secret',
        access_token = 'access_token',
        access_token_secret = 'access_token_secret'
    )
    
    #Se asignan cadenas que contendran las publicaciones
    etiquetas = '#biologymetrics'
    protocol_and_host = 'https://pubmed.ncbi.nlm.nih.gov/'
    
    #Se crea un conjunto al que se agregara los ids de los registros publicados
    ids_registros_publicados = set()
    
    #Se publican los registros disponibles
    for record in records:
        try:
            client.create_tweet(
                text = etiquetas + '\n' + record.get("TI") + '\n' + protocol_and_host + record.get("PMID") + '/'
            )
            ids_registros_publicados.add(record.get("PMID"))
        except:
            pass
    
    #Se escriben los ids de los registros publicados en un archivo
    with open('IDs_de_registros_publicados.txt', 'a') as handle:
        handle.writelines(x + '\n' for x in ids_registros_publicados)





if __name__ == "__main__":
    
    #Se obtienen las cadenas de las fechas
    from_date, to_date = get_dates()
    
    #Se obtienen una lista con los registros que se publicaran
    records = get_records(from_date, to_date)
    
    #Se hacen los tuits si hay registros para publicar
    if records:
        tweet_records(records)






