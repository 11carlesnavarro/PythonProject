from Bio.PDB import *
from Bio import pairwise2
# from Bio.SubsMat import MatrixInfo as matlist # Used for applying subtitution matrix in sequence comparison
from Bio.pairwise2 import format_alignment # DEBUGING Used for printing the alignment. DEBUGING purposes
from math import log # For use in defining the gap function for sequence comparison parameters
import sys
import os
import argparse
import re


###################################################################################################################################################

parser = argparse.ArgumentParser(description = """This program analyzes pdb binary chain interactions provided, calculate and reconstructs the polipeptide complex
represented by the provided files.
""")

parser.add_argument('-i', '--input',
                                        dest = "inPath",
                                        action = "store",
                                        default = os.getcwd(), # Default is current directory
                                        help = """Provide the complete path to the pdb files you
                                                  want to analyze after -i flag. If no path is provided,
                                                  the program will assume that current directory is
                                                  going to be searched for pdb files to analize.""")


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

###################################################################################################################################################


def get_files(input):
    """ This method reads the directory that is provided in the command line argument.
        Default directory to read is the current directory. It searches for files with
        .pdb extension and returns them.
    """

    path = input
    extension = re.compile('.pdb$')

    pdb_files = [file for file in os.listdir(path) if extension.search(file) is not None]
    os.chdir(path)

    return pdb_files

def get_file_prefix(pdb_files):
        """Tests for match of regex containing .pdb extension, trims it and returns
           the file prefix.
        """

        p = re.compile('(.*).pdb$')
        m = p.search(pdb_files)

        return m.group(1)


def compare_seqs(seq1, seq2):
    """This method aligns two given sequences, calculates percentage identity and
       tests if this identity is equal or higher than 95% or lower. If the given
       sequences show 95% identity or more, method returns True. otherwise it returns
       False.
    """
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


def get_structures(pdb_files):
    """This method parses the providad pdb files and extracts chains, and sequences trom them.
       The sequences are extracting using the alpha carbons provided in the pdb file. Input is
       a list of pdb files.
    """
    ppb = CaPPBuilder()
    parser = PDBParser(QUIET=True)

    for file in pdb_files:

        id = get_file_prefix(file)
        structure = parser.get_structure(id, file)
        # print(structure)
        sequences = []
        chains = []

        for chain in structure.get_chains():

            chains.append(chain)

        for pp in ppb.build_peptides(structure):

            seq = pp.get_sequence()
            sequences.append(seq)

        print(sequences) #debugging purposes


if __name__=="__main__":

    files = get_files(args.inPath)
    print(get_structures(files))
