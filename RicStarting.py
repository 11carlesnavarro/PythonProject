from Bio.PDB import *
from Bio import pairwise2
import sys
import os
import argparse
import re

parser = argparse.ArgumentParser(description="""This program reconstructs the macrocomplex from protein interactions pdb files""")

parser.add_argument('-i', '--input',
                                        dest = "infile",
                                        action = "store",
                                        default = os.getcwd(),
                                        help = "Enter a list of interaction files, a directory or the current directory will be selected")

parser.add_argument('-o', '--output',
                                        dest = "outfile",
                                        action = "store",
                                        default = None,
                                        help = "Enter the name of output file")

args = parser.parse_args()

###########################################################################################################################

#PROPIO#########################################################
#if  __name__=="__main__":
structure1 = PDBParser().get_structure('test', "chainA.pdb")
structure2 = PDBParser().get_structure('test', "chainB.pdb")

def get_sequence(structure):
    """returns the sequence of the structure instance passed
    through the method. It uses CaPPBuilder from Bio.pdb module
    so it identifies each alpha Carbon of the peptide chain to
    build the sequence"""

    ppb = CaPPBuilder()
    for pp in ppb.build_peptides(structure):
        seq = pp.get_sequence()

    return(seq)

def compare_seqs(seq1, seq2):
    """Returns the percentage sequence similarity between
    two given sequences. """

    alignment = pairwise2.align.globalxx(seq1, seq2)
    score = alignment[0][2]
    seq_sim = score/min(len(seq1), len(seq2))

    return(seq_sim)


#TESTING#################################################
######################
#
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

    seq1 = get_sequence(structure1)
    seq2 = get_sequence(structure2)

    compare_seqs(seq1, seq2)
