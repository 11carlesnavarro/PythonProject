
import sys
import os
import argparse
import re

parser = argparse.ArgumentParser(description = """This program analyzes pdb binary chain interactions provided, 
calculate and reconstructs the polipeptide complex represented by the provided files. """)

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

######### IncorrectInputFile Subclass ###########

class IncorrectInputFile(AttributeError):
    """ Exception when a letter that not belongs to the sequence alphabet is found"""
    __module__ = 'builtins'

    def __init__(self, pdbfile):
        self.pdbfile = pdbfile
    
    def __str__(self):
        """ Raise an exception if the file is not in the correct format"""

        return "The file %s is not appropiate" %(self.pdbfile)


####### Functions ########

def get_files(input):
    """ 
        This method reads the directory that is provided in the command line argument.
        Default directory to read is the current directory. It searches for files with
        .pdb extension and returns them.

    """
    # Print the execution of the program if -v is set
    if args.verbose:
        print("Reading pdb files...")

    # Reading pdb files and saving into a list

    path = input

    # Regular expression to check if the input file name is correct
    input_file = re.compile(r"^(?P<name>[a-zA-Z0-9]+)(\_)(?P<chain1>[a-zA-Z])(\_)(?P<chain2>[a-zA-Z])(.pdb$)")

    # List of pdb files
    pdb_files = []
    for f in os.listdir(path):
        if input_file.match(f):
            pdb_files.append(f)
        else:
            raise IncorrectInputFile(f)

    # Change dir to the dir of the files / input dir
    os.chdir(path)

    return pdb_files

def get_file_prefix(pdb_files):
        """
           Tests for match of regex containing .pdb extension, trims it and returns
           the file prefix.
        """

        p = re.compile('(.*).pdb$')
        m = p.search(pdb_files)

        return m.group(1)

if __name__ == "__main__":


    files = get_files(args.inPath)
    print(files[0])
    print(get_file_prefix(files[0]))
    
    #print(os.path.basename(args.inPath))
