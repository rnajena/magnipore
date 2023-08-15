#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import os
import pandas as pd

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('magnipore', type=str, help='.magnipore file to filter for a given coverage threshold')
    parser.add_argument('coverage', type=int, help='Coverage threshold to filter for. Results, where at least one sample has a coverage below the given threshold are filtered out. Results, where both samples have a coverage equal or higher than the threshold remain.')
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