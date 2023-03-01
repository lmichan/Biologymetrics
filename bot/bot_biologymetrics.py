#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import Entrez
from Bio import Medline
from datetime import date, timedelta
from mastodon import Mastodon
from posixpath import isfile




def get_dates():
    '''
    Obtiene una tupla con las cadenas de la fecha de inicio y la fecha de fin
    de la semana anterior a la semana en la que se encuentra el ultimo lunes.
    Las cadenas de las fechas tienen el formato empleado para hacer una consulta en Entrez.
    
    Returns
    -------
    tuple
        Las fechas de inicio y de fin de la semana anterior
        a la semana en la que se encuentra el ultimo lunes.
    '''
    
    #Se obtiene un objeto datetime.date que contiene la fecha actual
    current_date = date.today()
    
    #Se obtiene un objeto datetime.date que contiene la fecha del ultimo lunes
    last_monday_date = current_date - timedelta(days = current_date.weekday())
    
    #Se obtienen dos objetos datetime.date que contienen las fechas de inicio y de fin
    #de la semana anterior a la semana en la que se encuentra el ultimo lunes
    last_week_start_date = last_monday_date - timedelta(days=8)
    last_week_end_date = last_monday_date - timedelta(days=2)
    
    #Se obtienen las cadenas del ultimo par de objetos datetime.date creados
    #en el formato empleado para hacer la busqueda
    last_week_start_date_string = last_week_start_date.strftime('%Y/%m/%d')
    last_week_end_date_string = last_week_end_date.strftime('%Y/%m/%d')
    
    return last_week_start_date_string, last_week_end_date_string





def get_pubmed_records(term_search):
    '''
    Devuelve la informacion de los registros de PubMed que se obtiene con ``term_search``.
    
    Parameters
    ----------
    term_search : str
        Cadena del termino de busqueda.
    
    Returns
    -------
    list
        La informacion de los registros.
    '''
        
    #Se establece el parametro email de Entrez
    Entrez.email = 'A.N.Other@example.com'
    
    #Se ejecute una busqueda de Entrez y se obtiene un identificador a los resultados
    handle = Entrez.esearch(db='pubmed', term=term_search, retmax=500)
    
    #Se analiza un archivo XML de NCBI Entrez Utilities
    record = Entrez.read(handle)
    
    #Se crea una lista con los ids que se obtienen con la busqueda
    id_list = record['IdList']
    
    #Se obtienen los resultados de Entrez que son devueltos como un identificador
    handle = Entrez.efetch(db='pubmed', id=id_list, rettype='medline', retmode='text')
    
    #Se leen los registros de Medline desde el identificador y se guardan en una lista
    records = list(Medline.parse(handle))
    
    return records





def set_text(records):
    '''
    Devuelve un diccionario cuyos valores son cadenas con la informacion de ´´records´´
    con un formato para su publicacion en Mastodon y cuyas llaves son los PMIDs asociados.
    
    Parameters
    ----------
    records : list
        Informacion obtenida con la funcion get_pubmed_records de registros de PubMed.
    
    Returns
    -------
    dict(str, str)
        Diccionario cuyos valores son cadenas con la informacion de ´´records´´
        con un formato para su publicacion en Mastodon y cuyas llaves son los PMIDs asociados.
    '''
    
    record_dict = dict()
    protocol_and_host = 'https://pubmed.ncbi.nlm.nih.gov/'
    
    for record in records:
        labels = '#biologymetrics'
        if 'bibliometr' in record.get('TI', '').lower():
            labels += ' #bibliometrics'
        elif 'MH' in record and (('Bibliometrics' in record['MH']) or ('*Bibliometrics' in record['MH'])):
            labels += ' #bibliometrics'
        if 'altmetr' in record.get('TI', '').lower():
            labels += ' #altmetrics'
        if 'scientometr' in record.get('TI', '').lower():
            labels += ' #scientometrics'
        if 'PMID' in record:
            record_dict[record.get('PMID')] = labels + '\n' + record.get('TI', '?') + '\n' + protocol_and_host + record.get('PMID')
    
    return record_dict





def post_records(record_dict):
    '''
    Publica en Mastodon la informacion de los valores de ´´record_dict´´´
    que no se ha publicado ya y actualiza un archivo con los PMIDs
    de la informacion que se publique.
    
    Parameters
    ----------
    records : dict(str, str)
        Diccionario cuyos valores son cadenas con la informacion de registros de PubMed
        con un formato para su publicacion en Mastodon y cuyas llaves son los PMIDs asociados.
    
    Returns
    -------
    None.
    '''
    
    #Se lee de un archivo los ids de los registros que ya
    #se publicaron y se guardan en un conjunto
    path_pmids = 'path/Posted record PMIDs.txt'
    if isfile(path_pmids):
        with open(path_pmids) as f:
            previously_posted_record_ids = set([x.rstrip() for x in f])
    else:
        previously_posted_record_ids = set()
    
    #Se obtiene un agente de usuario para realizar solicitudes a la API de Mastodon
    mastodon = Mastodon(
        access_token = 'access_token',
        api_base_url = 'https://botsin.space/'
    )
    
    #Se crea un conjunto al que se agregara los ids de los registros que se publicaran
    currently_posted_record_ids = set()
    
    #Se publican los registros en Mastodon y se actualiza
    #un archivo con los ids de la informacion publicada
    try:
        for pmid in record_dict.keys() - previously_posted_record_ids:
            mastodon.status_post(record_dict[pmid])
            currently_posted_record_ids.add(pmid)
    finally:
        posted_record_updated_ids = currently_posted_record_ids | (record_dict.keys() & previously_posted_record_ids)
        with open(path_pmids, 'w') as handle:
            handle.writelines(pmid + '\n' for pmid in posted_record_updated_ids)







#%%

if __name__ == "__main__":
    
    #Se obtienen las cadenas de las fechas
    from_date, to_date = get_dates()
    
    #Se define el termino de busqueda de PubMed
    pubmed_search_term = '(altmetr*[Title] OR bibliometr*[Title] OR scientometr*[Title] OR bibliometrics[MeSH Terms]) \
            AND (\"' + from_date + '\"[Date - Create] : \"' + to_date + '\"[Date - Create])'
    
    #Se obtienen una lista con los registros actuales de PubMed
    pubmed_records = get_pubmed_records(pubmed_search_term)
    
    #Se obtienen un diccionario con la informacion de registros de PubMed
    record_dict = set_text(pubmed_records)
    
    #Se publica en Mastodon si hay informacion en ´´record_dict´´ sin publicar
    post_records(record_dict)



