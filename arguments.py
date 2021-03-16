
import sys
import os
import argparse
import re
import shutil


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

parser.add_argument('-f', '--force',
                                        dest = "force",
                                        action = "store_true",
                                        default = False,
                                        help = """If this option is False and the output directory already exists before the application is
                                                   executed, exit the program execution and warn the user that the directory already exists. 
                                                   If it is True, then the program can continue running and overwrite all the contents of the
                                                   output directory""")

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

# READ FILES

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

    

    return pdb_files


# OUTPUT DIRECTORY 

def save_output():
    """
        If the output directory does not exists, create it. 
    """

    # Verbose
    if args.verbose:
        print("Writing the output...")

    subfolder_names = ["Structures", "Analysis"]

    # If the argument force is not selected create the output directory if it does not exists 

    if not args.force:
        #if not os.path.exists('FinalComplex'):
        try:
            for subfolder_name in subfolder_names:
                os.makedirs(os.path.join('FinalComplex', subfolder_name))
        except OSError as err:
            raise err
    
    # If the argument force is selected create the output directory if it does not exists and, if exists, override it
    elif args.force:
        if not os.path.exists('FinalComplex'):
            for subfolder_name in subfolder_names:
                os.makedirs(os.path.join('FinalComplex', subfolder_name))
        else:
            #Remove the directory
            shutil.rmtree('FinalComplex')
            for subfolder_name in subfolder_names:
                os.makedirs(os.path.join('FinalComplex', subfolder_name))



if __name__ == "__main__":


    files = get_files(args.inPath)
    print(files[0])
    print(get_file_prefix(files[0]))
    
    #print(os.path.basename(args.inPath))
    save_output()