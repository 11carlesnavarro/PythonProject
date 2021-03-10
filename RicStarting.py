from Bio.PDB import *
from Bio import pairwise2
# from Bio.SubsMat import MatrixInfo as matlist # Used for applying subtitution matrix in sequence comparison
from Bio.pairwise2 import format_alignment # DEBUGING Used for printing the alignment. DEBUGING purposes
from math import log # For use in defining the gap function for sequence comparison parameters
import sys
import os
import argparse
import re

arser = argparse.ArgumentParser(description = """This program analyze pdb binary chain interactions provided, calculate and reconstructs the polipeptide complex
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
                                        help = """The output file will be named after the name provided after this flag.
                                                This option is mandatory, if no output filename is provided, no analysis
                                                will be run and this message will be displayed. """)

parser.add_argument('-v', '--verbose',
                                        dest = "verbose",
                                        action = "store_true",
                                        default = False,
                                        help = """This option prints the execution of the program as it is the actions are being.
                                                performed. It is used without options as it activates the verbose mode on.
                                                It is adviseable to run the analysis with this option in case any exception is raised.""")

args = parser.parse_args()

###########################################################################################################################

# Por ahora la ejecución del programa funciona bien, pero toma solo éstas dos secuencias. Debe cambiarse para operar sobre 
# todos los archivos del directorio y desde los archivos pareados, analizándolos sin separar.
structure1 = PDBParser().get_structure('test', "chainA.pdb")
structure2 = PDBParser().get_structure('test', "chainC.pdb")

# regresa la secuencia de un archivo pdb. Debe trabajar sobre un solo archivo pareado.
def get_sequence(structure):
    """returns the sequence of the structure instance passed
       through the method. It uses CaPPBuilder from Bio.pdb module
       so it identifies each alpha Carbon of the peptide chain to
       build the sequence"""

    ppb = CaPPBuilder()
    for pp in ppb.build_peptides(structure):
        seq = pp.get_sequence()

    return(seq)

# Funciona bien, debe recibir dos secuencias separadas. El cálculo de la identidad de secuencias es correcto.
def compare_seqs(seq1, seq2):
    """This method aligns two given sequences, calculates percentage identity and
       tests if this identity is equal or higher than 95% or lower. If the given
       sequences show 95% identity or more, method returns True. otherwise it returns
       False."""
    # value 1 is for counting identical matches. Therefore the alignment score is
    # equal to counting identical matches. percentage identity is then calculated
    # as the ratio between the score and the length of the alignment. this percentage
    # identity is the BLAST definition.

    # definition of BLAST percentage identity is taken from:
    # Heng Li, 2021, at: https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
    align = pairwise2.align.globalmx(seq1, seq2, 1, 0)
    algn_len = len(align[0][1])
    match_score = align[0][2]

    identity = match_score/algn_len

    # Checking for identity, true = homodimer. false = heterodimer
    if identity >= .95:
        return True
    else:
        return False


if __name__=="__main__":

    seq1 = get_sequence(structure1)
    seq2 = get_sequence(structure2)

    compare = compare_seqs(seq1, seq2)
    print(compare)
