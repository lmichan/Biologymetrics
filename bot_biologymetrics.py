#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import Entrez
from Bio import Medline
from datetime import date, timedelta
from mastodon import Mastodon
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





def get_pubmed_records(term_search: str) -> list:
    '''
    Obtiene la informacion de los registros creados en PubMed que se obtienen
    con ``term_search``.
    
    Parameters
    ----------
    term_search : str
        adena del termino de busqueda.
    
    Returns
    -------
    list
        La informacion de los registros.
    
    '''
        
    #Se establece el parametro email de Entrez
    Entrez.email = 'A.N.Other@example.com'
    
    #Se ejecute una búsqueda de Entrez y se obtiene un identificador a los resultados
    handle = Entrez.esearch(db='pubmed', term=term_search, retmax=500)
    
    #Se analiza un archivo XML de NCBI Entrez Utilities en objetos de python
    record = Entrez.read(handle)
    
    #Se crea una lista con los ids que se obtienen con la busqueda
    id_list = record['IdList']
    
    #Se obtienen los resultados de Entrez que son devueltos como un identificador
    handle = Entrez.efetch(db='pubmed', id=id_list, rettype='medline', retmode='text')
    
    #Se leen los registros de Medline desde el identificador y se guardan en una lista
    records = list(Medline.parse(handle))
    
    return records





def set_text(records: list) -> dict:
    '''
    Obtiene de una lista con la informacion de registros de PubMed, un diccionario
    cuyos valores son el texto que se publicara y cuyas llaves son los PMIDs asociados.

    Parameters
    ----------
    records : list
        La informacion de los registros de PubMed.

    Returns
    -------
    dict
        Diccionario cuyas llaves son PMIDs y cuyos valores son cadenas
        de la informacion que se publicara.

    '''
    
    hashtag = '#biologymetrics'
    protocol_and_host = 'https://pubmed.ncbi.nlm.nih.gov/'
    
    return {record.get('PMID'):
            hashtag + '\n' + record.get('TI', '?') + '\n' + protocol_and_host + record.get('PMID') + '/'
            for record in records}




def post_records(record_dict: dict):
    '''
    Publica en Mastodon la informacion de los valores del diccionario y
    escribe en los PMIDs asociados en un archivo.
    
    Parameters
    ----------
    records : dict
        La informacion de los registros que se publicaran y sus PMIDs sociados.
    
    Returns
    -------
    None.
    
    '''
    
    #Si es lunes se borra el contenido de un archivo que contiene
    #los ids de PubMed de los registros publicados en Mastodon
    if date.today().weekday() == 0:
        open('IDs_de_registros_publicados_en_mastodon.txt', 'w').close()
    
    #Se lee de un archivo los ids de los registros que ya
    #se publicaron y estos se guardan en un conjunto
    if isfile('IDs_de_registros_publicados_en_mastodon.txt'):
        with open('IDs_de_registros_publicados_en_mastodon.txt') as f:
            posted_record_ids = set([x.rstrip() for x in f])
    else:
        posted_record_ids = set()
    
    #Se obtiene un agente de usuario utilizado al realizar solicitudes a la API de Mastodon
    mastodon = Mastodon(
        access_token = 'access_token',
        api_base_url = 'https://botsin.space/'
    )
    
    #Se crea un conjunto al que se agregara los ids de los registros publicados
    ids = set()
    
    #Se publican los registros disponibles
    for pmid, info in record_dict.items():
        if pmid not in posted_record_ids:
            try:
                mastodon.status_post(info)
                ids.add(pmid)
            except:
                pass
    
    #Se escriben los ids de los registros publicados en un archivo
    with open('IDs_de_registros_publicados_en_mastodon.txt', 'a') as handle:
        handle.writelines(x + '\n' for x in ids)




def tweet_records(record_dict: dict):
    '''
    Publica tuits con la informacion de los valores del diccionario y
    escribe en los PMIDs asociados en un archivo.
    
    Parameters
    ----------
    records : dict
        La informacion de los registros que se publicaran y sus PMIDs sociados.
    
    Returns
    -------
    None.
    
    '''
    
    #Si es lunes se borra el contenido de un archivo que contiene 
    #los ids de PubMed de los registros publicados en Twitter
    if date.today().weekday() == 0:
        open('IDs_de_registros_publicados_en_twitter.txt', 'w').close()
    
    #Se lee de un archivo los ids de los registros que ya
    #se publicaron y estos se guardan en un conjunto
    if isfile('IDs_de_registros_publicados_en_twitter.txt'):
        with open('IDs_de_registros_publicados_en_twitter.txt') as f:
            posted_record_ids = set([x.rstrip() for x in f])
    else:
        posted_record_ids = set()
    
    #Se obtiene un agente de usuario utilizado al realizar solicitudes a la API de Twitter
    client = tweepy.Client(
        consumer_key = 'consumer_key',
        consumer_secret = 'consumer_secret',
        access_token = 'access_token',
        access_token_secret = 'access_token_secret'
    )
    
    #Se crea un conjunto al que se agregara los ids de los registros publicados
    ids = set()
    
    #Se publican los registros disponibles
    for pmid, info in record_dict.items():
        if pmid not in posted_record_ids:
            try:
                client.create_tweet(text=info)
                ids.add(pmid)
            except:
                pass
    
    #Se escriben los ids de los registros publicados en un archivo
    with open('IDs_de_registros_publicados_en_twitter.txt', 'a') as handle:
        handle.writelines(x + '\n' for x in ids)








if __name__ == "__main__":
    
    #Se obtienen las cadenas de las fechas
    from_date, to_date = get_dates()
    
    #Se define el termino de busqueda de PubMed
    pubmed_search_term = '(altmetr*[Title] OR bibliometr*[Title] OR scientometr*[Title] OR bibliometrics[MeSH Terms]) \
            AND (\"' + from_date + '\"[Date - Create] : \"' + to_date + '\"[Date - Create])'
    
    #Se obtienen una lista con los registros actuales de PubMed
    pubmed_records = get_pubmed_records(pubmed_search_term)
    
    #Se obtienen un diccionario con la informacion que se publicara
    record_dict = set_text(pubmed_records)
    
    #Se publica en Mastodon si hay registros para publicar
    post_records(record_dict)

    #Se publica en Twitter si hay registros para publicar
    tweet_records(record_dict)



