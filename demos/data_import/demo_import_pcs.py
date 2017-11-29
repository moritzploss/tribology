"""

Short demonstration of how to use data import methods to import data files from
PCS Instruments test rigs. This script is meant to be run from a terminal and
can be used to import all files in the current working directory and---depending
on parser arguments---all its child directories.

"""


import argparse
import os

import tribology as tr


def __parse_args():
    """

    Parse all parser arguments that are provided when the script is running in
    a terminal.

    Returns
    -------
    args: Namespace
        The parsed parser arguments.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--force', action="store_true", default=False,
                        help='overwrite existing database files during import')
    parser.add_argument('-e', '--extension', action="store", default='txt',
                        help='specify file extension. default is "txt"')
    parser.add_argument('-d', '--delimiter', action="store", default='\t',
                        help='specify column delimiter. default is tab (\\t)')
    parser.add_argument('-m', '--mark', action="store", default='.',
                        help='specify decimal mark for numeric data. default is'
                             ' dot (.)')
    parser.add_argument('-o', '--outformat', action="store", default='npz',
                        help='specify output database format. default is "npz"'
                             ' for numpy database. use "mat" for matlab '
                             ' database format.')
    parser.add_argument('-r', '--recursive', action="store_true", default=False,
                        help='recursively walk through all sub-directories of'
                             ' current working directory')
    parser.add_argument('-p', '--pcs', action="store_true", default=True,
                        help='indicate if files are pcs files.')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    """
    
    Use the code:`import_dir` method of the tribology package to import all 
    data files based on the parser arguments provided.
    
    """
    ARGS = __parse_args()
    tr.import_dir(os.getcwd(), in_ext=ARGS.extension, recursive=ARGS.recursive,
                  force=ARGS.force, deli=ARGS.delimiter, dec_mark=ARGS.mark,
                  out_ext=ARGS.outformat, out_dir='', print_stat=True,
                  pcs=ARGS.pcs)
