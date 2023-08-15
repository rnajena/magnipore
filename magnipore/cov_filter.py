#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from magnipore.__init__ import __version__
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import os
import pandas as pd

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description='A small script to filter the .magnipore output for a given coverage threshold. Please provide the path to the .magnipore file and a coverage threshold. All compared positions, where at least one sample has a coverage below your given threshold will be filtered out. The remaining positions are written to a file with the name of the provided magnipore file and the suffix \"_c10\": e.g. <given_magnipore_file>_c10.magnipore',
        prog='Magnipore',
    )
    parser.add_argument('magnipore', type=str, help='.magnipore file to filter for a given coverage threshold')
    parser.add_argument('coverage', type=int, help='Coverage threshold to filter for. Results, where at least one sample has a coverage below the given threshold are filtered out. Results, where both samples have a coverage equal or higher than the threshold remain.')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s' + f' {__version__}')
    return parser.parse_args()

def filter_coverage(magnipore_file : str, cov_threshold : int) -> None:
    magnipore = pd.read_csv(magnipore_file, sep='\t')
    outfile = os.path.splitext(magnipore_file)[0] + '_c10.magnipore'
    filtered_magnipore = magnipore[(magnipore['n_reads_1'] >= cov_threshold) & (magnipore['n_reads_2'] >= cov_threshold)]
    filtered_magnipore.to_csv(outfile, sep='\t')

def main() -> None:
    args = parse()
    filter_coverage(args.magnipore, args.coverage)

if __name__ == '__main__':
    main()