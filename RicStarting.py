from Bio.PDB import *
from Bio import pairwise2
import sys
import os
import argparse
import re

parser = argparse.ArgumentParser(description = """This program analyze pdb binary chain interactions provided, calculate and reconstructs the polipeptide complex
represented by the provided files. 
""")

parser.add_argument('-i', '--input',
                                        dest = "infile",
                                        action = "store",
                                        default = os.getcwd(), #current input default is complete current directory, we need to discuss what default will be set.
                                        help = """You need to provide the path to the directory where the binary interaction pdb files are stored, 
                                                  they will be all analized""") # I think we could use default = None, and guide the user to provide only a path. Discussion.

parser.add_argument('-o', '--output',
                                        dest = "outfile",
                                        action = "store",
                                        default = None,
                                        help = "Enter the name of output file")

args = parser.parse_args()

###########################################################################################################################

#PROPIO#########################################################

# En la investigación vi que exiten por lo menos tres instancias de interés para analizar los pdb:
# Structure, Chain y sequence. Structure está representada por la totalidad del pdb file, chain se
# refiere al nombre de la cadena, en los binary interaction files habrá solo dos cadenas. Sequence 
# es como su nombre lo dice la sequencia. En este caso de aminoácidos.

# Basados en el texto anterior, estas dos líneas siguientes analizan los archivos pdb y almacenan la 
# información de la estructura en general en objetos de python.
structure1 = PDBParser().get_structure('test', "chainA.pdb")
structure2 = PDBParser().get_structure('test', "chainB.pdb")

# obtener la secuencia. El texto en la función es propio.
def get_sequence(structure):
    """returns the sequence of the structure instance passed
    through the method. It uses CaPPBuilder from Bio.pdb module
    so it identifies each alpha Carbon of the peptide chain to
    build the sequence"""

    ppb = CaPPBuilder()
    for pp in ppb.build_peptides(structure):
        seq = pp.get_sequence()

    return(seq)

# En ésta función se hace un alineamiento global. después se calcula el sequence simlarity (SS)
# con la fórmula   SS = score/len(minSeq)   en donde score es el valor del score obtenido en la 
# alineación pareada. OJO !!! - éste cálculo es erroneo, se está construyendo esta función, estoy
# investigando la fórmula necesarioa para calcular el SS de forma correcta. Se tiene que penalizar
# de forma diferente el cambio de aminoácidos según las propiedades fisicoquímicas de los camvbios de
# los residuos
def compare_seqs(seq1, seq2):
    """Returns the percentage sequence similarity between
    two given sequences. """

    alignment = pairwise2.align.globalxx(seq1, seq2)
    score = alignment[0][2]
    seq_sim = score/min(len(seq1), len(seq2))

    return(seq_sim)


#TESTING#################################################
######################
## éstas son funciones en construcción, no funcionan aún, y están siendo exploradas. Las funciones de 
## éstas vienen del paquete Bio.pdb

# def get_name_structure(pdb_files):
#         """Parses the names of the input files and returns the name of the files without .pdb extension"""
#         p = re.compile('(.*).pdb') #using regular expression, impor re needed
#         m = p.match(pdb_files)
#
#         return m.group(1)
#
# def get_pdb_info(pdb_files):
#     p = PDBParser(PERMISSIVE=1, QUIET=True)
#     error_files = []
#     for f in pdb_files:
#         try:
#             str_id = get_name_structure(f)
#             print(str_id)
#             str = p.get_structure(str_id,f)
#         except:
#             error_files.append(f)
#             continue
#
#         chains = []
#         sequences = []
#         chains_to_remove = []
#
#         for chain in str.get_chains():
#             seq = get_sequence(chain)
#             if seq is None:
#                 chains_to_remove.append(chain.id)
#
#             else:
#                 sequences.append(seq)
#                 chains.append(chain)
######################TESTING###############################3

if  __name__=="__main__":

    # corriendo las dos funciones. ambas funcionan, se extraen las secuencias de las dos cadenas 
    # (vienen de archivos independientes, hay que desarrollar la forma de hacerlo desde un solo archivo binario)
    # tras extraer las secuencias, se hace el alineamiento pareado y se ejecuta la formula. Estoy trabajando 
    # en obtener la fórmula correcta.
    
    seq1 = get_sequence(structure1)
    seq2 = get_sequence(structure2)

    compare_seqs(seq1, seq2)
