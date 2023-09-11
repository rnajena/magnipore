#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import pandas as pd
import numpy as np

from magnipore.__init__ import __version__

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description='A small script that compares the Magnipore output file with a given validationset in form of a table in a CSV file. Your table should contain the reference id/name, that was also used during the Magnipore analysis, and the validation positions on the given reference. The script will then check, if Magnipore found significant positions in a kmer range of the positions in the validation table (CSV). The kmer range depends on the used pore during sequencing. You can specify the used pore with the --pore parameter. You should also think about the coverage threshold. This script will by default filter out positions from the Magnipore output, where at least one sample has a coverage less than 10 reads.',
        prog='magnicheck',
    )
    parser.add_argument('magnipore', type=str, help='Magnipore file with called differential modifications.')
    parser.add_argument('magnipore_poscol', type=str, help='Which position to validate from the Magnipore output.')
    parser.add_argument('refid', type=str, help='Reference id or name - This id must match between the magnipore file and the your validation table (CSV).')
    parser.add_argument('eval_csv', type=str, help='CSV file containing the validation table. The table contains the reference id and position of validated (ground truth) modifications.')
    parser.add_argument('refcol', type=str, help='Column containing the reference ids.')
    parser.add_argument('poscol', type=str, help='Column containing the validated positions.')
    parser.add_argument('outfile', type=str, help='Name of the output file.')
    parser.add_argument('--coverage', type=int, default=10, help='Coverage filter to apply to the Magnipore output.')
    parser.add_argument('--valid_sep', type=str, default=',', help='Separation character in your CSV file.')
    parser.add_argument('--pore', type=str, default='r9', choices=['r9', 'r10'], help='Which pore you used during sequencing.')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s' + f' {__version__}')
    return parser.parse_args()

def main() -> None:
    args = parse()

    eval = pd.read_csv(args.eval_csv, sep=args.valid_sep)
    eval = eval[eval[args.refcol] == args.refid]
    magnipore_pos = pd.read_csv(args.magnipore, sep='\t')
    print(f'Number of called positions by magnipore: {len(magnipore_pos.index)}')
    magnipore_pos = magnipore_pos[(magnipore_pos['n_reads_1'] >= args.coverage) & (magnipore_pos['n_reads_2'] >= args.coverage)]
    print(f'Number of called positions after coverage filtering with {args.coverage} by magnipore: {len(magnipore_pos.index)}')
    
    if args.pore == 'r9':
        r = 5
    elif args.pore == 'r10':
        r = 7

    m = 0
    f = 0

    with open(args.outfile, 'w') as w:
        w.write('strand,found_positions\n')

        for strand in '+-':
            eval_pos = eval[eval['strand'] == strand][args.poscol].to_numpy()
            called_pos = magnipore_pos[magnipore_pos['strand'] == strand][args.magnipore_poscol].to_numpy()

        # found_pos = np.intersect1d(eval_pos, magnipore_pos, assume_unique=True)
        # print(f'directly found positions: {len(found_pos)}/{len(eval_pos)}')
        # print('fraction:', len(found_pos)/len(eval_pos))

            magniporeRange = []
            for pos in called_pos:
                for i in range(pos-r//2, pos+r//2+1):
                    magniporeRange.append(i)
            magniporeRange = np.unique(magniporeRange)

            within_kmer = np.intersect1d(eval_pos, magniporeRange, assume_unique=True)

            print(f'Found positions within {r}mer range {len(within_kmer)}/{len(eval_pos)} for strand: {strand}')
            print('Fraction:', len(within_kmer)/len(eval_pos))

            f += len(within_kmer)
            m += len(eval_pos)

            for pos in within_kmer:
                w.write(f'{strand},{pos}\n')

    print('Overall')
    print(f'Found positions within {r}mer range {f}/{m}')

if __name__ == '__main__':
    main()