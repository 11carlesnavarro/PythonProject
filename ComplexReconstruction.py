
import arguments as arguments
import sys
import os
import argparse
import re
import shutil
from Bio.PDB import *
from Bio import pairwise2
# from Bio.SubsMat import MatrixInfo as matlist # Used for applying subtitution matrix in sequence comparison
from Bio.pairwise2 import format_alignment # DEBUGING Used for printing the alignment. DEBUGING purposes



########## GET FILE PREFIXES ###########

def get_file_prefix(pdb_files):
        """Tests for match of regex containing .pdb extension, trims it and returns
           the file prefix.
        """

        p = re.compile('(.*).pdb$')
        m = p.search(pdb_files)

        return m.group(1)


########### COMPARE SEQUENCES ###########

def compare_seqs(seq1, seq2):
    """
       This method aligns two given sequences, calculates percentage identity and
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


######### GET THE SEQUENCES OF EACH FILE ###########

def get_structures(pdb_files, path = arguments.args.inPath):
    """
       This method parses the provided pdb files and extracts and yields the structures 
       from them.
    """

    parser = PDBParser(QUIET=True)

    # Main directory 
    main_dir = os.getcwd()

    # Change dir to the dir of the files / input dir
    os.chdir(path)


    for file in pdb_files:

        id = get_file_prefix(file)
        yield parser.get_structure(id, file)
        #print(structure)

        # Set the path again to main

    os.chdir(main_dir)


def get_chains_structure(structure):
        """ 
        Returns a list of chains given a structure

        """
        chains = []

        for model in structure:
            for chain in model:
                chains.append(chain)
        return chains

def get_sequences(structure):

    """
    Returns a list of sequences given a structure
    """

    ppb = PPBuilder()
    sequences = []

    for pp in ppb.build_peptides(structure):

        seq = pp.get_sequence()
        sequences.append(seq)

    return sequences #debugging purposes
    
def get_sequences_string(chain):
    """
    Given a chain it returns its sequence as string

    """

    ppb = PPBuilder()
    for pp in ppb.build_peptides(chain):
        return pp.get_sequence()



if __name__=="__main__":


    files = arguments.get_files(arguments.args.inPath)
    files_dir = arguments.args.inPath


    print(files[0])
    print(get_file_prefix(files[0]))
    
    get_structures(files, files_dir)

    print(arguments.args)   
    print(os.getcwd())



    for structure in get_structures(files):
        chains = get_chains_structure(structure)
        n = 0
        for chain in chains:
            print(chains[n])
            print(get_sequences_string(chain))
            n += 1

    for structure in get_structures(files):
        print(structure)
        print(get_sequences(structure))

#######################################################################

