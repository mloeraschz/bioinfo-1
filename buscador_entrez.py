# -*- encoding: utf-8 -*-


from Bio import Entrez
import sys
import time
Entrez.email = "bio.mloera@gmail.com"

## Buscador para obtener IdList con usehistory en ON

search_handle = Entrez.esearch(db="gene",term="hexokinase AND Aspergillus[Orgn]", usehistory="y")
search_results = Entrez.read(search_handle)

search_handle.close()
gi_list = search_results["IdList"]
count = int(search_results["Count"])
#assert count == len(gi_list)
webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"]
print(webenv,query_key)

## Definimos el tamaño del batch con el que se hará .efetch()
start = 0
batch_size = 3

## En este archivo se van a guardar las secuencias que se están bajando

out_handle = open("ejercicios/hexokinase_aspergillus.fasta", "a")

## Empieza el loop para hacer .efetch() a todos los elementos del IdList

for start in range(0,len(gi_list),batch_size):
    time.sleep(3)
    fetch_handle = Entrez.efetch(db="gene", rettype="fasta", retmode="text",
                                         retstart=start, retmax=batch_size,
                                         webenv=webenv, query_key=query_key)
        
    data = fetch_handle.read()
    fetch_handle.close()
    out_handle.write(data)
out_handle.close()
