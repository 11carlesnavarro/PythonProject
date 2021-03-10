from Bio.PDB import *
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist # Used for applying subtitution matrix in sequence comparison
from Bio.pairwise2 import format_alignment # DEBUGING Used for printing the alignment. DEBUGING purposes
from math import log # For use in defining the gap function for sequence comparison parameters
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
  
 def gap_opening_function(x, y): # x is gap position in seq, y is gap length
    """ This function handles the penalties of gap openings and extensions in
    compare_seqsBLOSUM method. change the values on the return funtions for changing
    the penalty values. """

    if y == 0: # No gap
        return 0
    elif y == 1: #Gap open penalty
        return -2
    return - (2 + y/4.0 + log(y)/2.0)

def gap_extension_function(x, y): # x is gap position in seq, y is gap length
    """ This function handles the penalties of gap openings and extensions in
    compare_seqsBLOSUM method. change the values on the return funtions for changing
    the penalty values. """

    if y == 0: # No gap
        return 0
    elif y == 1: #Gap open penalty
        return -2
    return - (2 + y/4.0 + log(y)/2.0)

# En ésta función se hace un alineamiento global. después se calcula el sequence simlarity (SS)
# con la fórmula   SS = score/len(maxSeq)   en donde score es el valor del score obtenido en la 
# alineación pareada. se utiliza la matriz de substitución BLOSUM50, pero puede ser elegida alguna
# otra. Los penalties de gap opening y gap extension pueden ser cambiados en las funciones de:
# gap_opening_function y gap_extension_function

# OJO!!, aún falta confirmar si aunque tengamos matriz de substitución y penalties, es pertinente
# dividir el score entre la secuencia de máxima longitud para obtener el porcentaje de similitud.
def compare_seqsBLOSUM50(seq1, seq2):
    """This method returns the score of a pairwise sequence comparison. Expects
    two sequences. The matrix substitution is blosum50, and penalizes gap openings
    and extensions with values set in gap_function. default penalties are :

    gap opening: -2 # I am exploring which penalties are best for gap opening and extension.
    gap extension: -2 """

    matrix = matlist.blosum50
    alignment = pairwise2.align.globaldc(seq1, seq2, matrix, gap_opening_function, gap_extension_function)
    score = alignment[0][2]
    seq_sim = score/max(len(seq1), len(seq2))

    return seq_sim

if  __name__=="__main__":

    # corriendo las funciones. corren sin error, pero con un warning que estoy trabajando. se extraen las secuencias de las dos cadenas 
    # (vienen de archivos independientes, hay que desarrollar la forma de hacerlo desde un solo archivo binario)
    # tras extraer las secuencias, se hace el alineamiento pareado y se ejecuta la formula. Estoy trabajando 
    # en obtener la fórmula correcta.
    
    seq1 = get_sequence(structure1)
    seq2 = get_sequence(structure2)

    compare_seqs(seq1, seq2)
