### BIOPYTHON
###
### "El proyecto Biopython es una asociación internacional de desarrolladores
## de herramientas para biología molecular computacional" 
## (http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc2)

## Entre las herramientas que pueden encontrarse en Python están:
## 1. Parsers (interpretadores) de estructuras de datos usados en bioinformática
## 2. Acceso a bases de datos bionformáticas 
## 3. Interfaces para programas comunes en bioinformática
## 4. Estructuras de datos (clases) re-usables para secuencias y métodos
## asociados a ellas

## TEMA 1. PARSERS DE ESTRUCTURAS DE DATOS COMUNES EN BIOINFORMÁTICA

## Un parser, o interpretador,  permite acceder y usar la información que
## se encuentra codificada en una estructura de datos.

## Entre las estructuras de datos que Biopython puede interpretar, se encuentran:
## 1. Blast output (XML)
## 2. Clustalw
## 3. Fasta
## 4. GenBank
## 5. PubMed and Medline
## 6. ExPASy files, like Enzyme and Prosite
## 7. SCOP, including ‘dom’ and ‘lin’ files
## 8. UniGene
## 9. SwissProt


## 1.1 PARSER DE ARCHIVOS FASTA

from Bio import SeqIO # Primero hay que importar Biopython y la clase SeqIO

## Cuando se trata de un archivo con una sola secuencia, puede usarse:
## Bio.SeqIO.read()
## En cambio, cuando se trata de un archivo multifasta, debe usarse:
## Bio.SeqIO.parse()

for elemento in SeqIO.parse('dna.fasta','fasta'):

    print(elemento.id) # El ID del elemento
    print(elemento.description)
    #print(elemento.seq) # La secuencia del elemento
    #print(repr(elemento.seq)) # Descripción breve de la secuencia
    #print(len(elemento)) # Longitud de la secuencia

## RETO 1: Cambia el nombre de las secuencias del archivo "dna.fasta" para
## que solamente quede el identificador EMB.
##
## RETO 2: Cambia el nombre de las secuencias del archivo "dna.fasta" para
## que solamente quede el nombre de la especie. (Recordatorio de Python básico)
##
## RETO 3: Encuentra las posiciones de todos los codones de inicio y de stop
## en las secuencias del archivo "dna.fasta"


## 1.2 PARSER DE ARCHIVOS GENBANK

elemento = SeqIO.read('NC_005816.gb','genbank')

for feature in elemento.features:
    print(feature.type, feature.qualifiers.get('db_xref'), feature.location
                    )
    

## RETO 3: extrae los "qualifiers" de solamente los CDS del archivo 'NC_005816.gb'
## usando ".qualifiers.get('db_xref')"


## 1.3 PARSING DE DATOS GENBANK DESDE EL INTERNET

from Bio import Entrez
Entrez.email = 'bio.mloera@gmail.com'
handle = Entrez.efetch(db='nucleotide',rettype='gb',retmode='text',id='6273291,6273290,6273289')
for elemento in SeqIO.parse(handle,'gb'):
    print("%s %s..." % (elemento.id, elemento.description[:50]))
    print("Sequence length %i, %i features, from: %s"
          % (len(elemento), len(elemento.features), elemento.annotations["source"]))
handle.close()

## 1.4 PARSING DE DATOS DE SWISSPROT DESDE EL INTERNET

from Bio import ExPASy
from Bio import SeqIO
handle = ExPASy.get_sprot_raw("Q197F7")
seq_record = SeqIO.read(handle, "swiss")
handle.close()
print(seq_record.id)
print(seq_record.name)
print(seq_record.description)
print(seq_record.seq)
print(seq_record.annotations["keywords"])
