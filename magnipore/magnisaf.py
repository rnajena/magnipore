#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import pandas as pd

SAMPLE : int

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description='Checks if magnipore positions are within annotated regions of the provided saf file.',
        prog='magnisaf',
    )
    parser.add_argument('saf', type=str, help='.saf file containing chromosome, gene, region and strand')
    parser.add_argument('magnipore', type=str, help='.magnipore file')
    parser.add_argument('outfile', type=str, metavar='TSV', help='TSV style output table')
    parser.add_argument('--sample', default=1, type=int, choices=[1, 2], help='Which sample to compare in magnipore file.')
    return parser.parse_args()

def readFile(file : str, sep : str) -> pd.DataFrame:
    return pd.read_csv(file, sep=sep)

def f_match(saf_row : pd.Series, magni_row : pd.Series) -> bool:
    if not saf_row['Chr'] == magni_row[f'ref_{SAMPLE}']:
        return False
    if not saf_row['Strand'] == magni_row['strand']:
        return False
    return True

def getGenes(saf : pd.DataFrame, magnipore : pd.DataFrame) -> pd.DataFrame:
    ret = pd.DataFrame(columns=['Chr', 'Strand', 'GeneID', 'Start', 'End', 'Magnipore'])
    ret = ret.astype({'Chr':str, 'Strand':str, 'GeneID':str, 'Start':int, 'End':int, 'Magnipore':int})
    print('Start comparing files')
    for idx, magni_entry in magnipore.iterrows():
        for _, saf_entry in saf.loc[(saf['Start']<=magni_entry[f'pos_{SAMPLE}']) & (magni_entry[f'pos_{SAMPLE}']<=saf['End'])].iterrows():
            if (idx + 1) % 100 == 0:
                print(f'Processing magnipore line {idx + 1}/{len(magnipore.index)}', end = '\r')
            if f_match(saf_entry, magni_entry):
                # new_entry = pd.DataFrame(saf_entry._append(pd.Series({'Magnipore': magni_entry[f'pos_{SAMPLE}']}))).transpose()
                new_entry = pd.DataFrame({
                    'Chr':[saf_entry['Chr']],
                    'Strand':[saf_entry['Strand']],
                    'GeneID':[saf_entry['GeneID']],
                    'Start':[saf_entry['Start']],
                    'End':[saf_entry['End']],
                    'Magnipore':[magni_entry[f'pos_{SAMPLE}']],
                })
                ret = pd.concat((ret, new_entry), ignore_index=True)
    print(f'Processed magnipore line {idx+1}/{len(magnipore.index)}     ')
    ret = ret.set_index('GeneID')
    return ret

def main() -> None:
    args = parse()
    global SAMPLE
    SAMPLE = args.sample
    saf = readFile(args.saf, '\t')
    magnipore = readFile(args.magnipore, '\t')
    output = getGenes(saf, magnipore)
    output.to_csv(args.outfile, sep='\t', index=False)

if __name__ == '__main__':
    main()