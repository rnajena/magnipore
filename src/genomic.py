#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import pandas as pd
from tqdm import tqdm
from src.__init__ import __version_str__

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description='Checks if magnipore positions are within annotated regions of the provided annotation file.',
        prog='magnipore genomic',
    )
    parser.add_argument('annot', type=str, help='.saf, .gff, or .gtf file containing chromosome, gene, region and strand')
    parser.add_argument('magnipore', type=str, help='.magnipore file')
    parser.add_argument('outfile', type=str, help='Output file will be in .tsv format.')
    parser.add_argument('--sample', default=1, type=int, choices=[1, 2], help='Which sample to compare in magnipore file.')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s' + f' {__version_str__}')
    return parser.parse_args()

import pandas as pd

def read_annot(file: str) -> pd.DataFrame:
    """
    Reads an annotation file (SAF, GTF, or GFF) and returns a pandas DataFrame.

    Parameters:
        file (str): Path to the annotation file.

    Returns:
        pd.DataFrame: Parsed annotation data.
    """
    if file.endswith('.saf'):
        # SAF files are tab-separated with 5 required columns: GeneID, Chr, Start, End, Strand
        df = pd.read_csv(file, sep="\t", comment="#")
        required_columns = {"GeneID", "Chr", "Start", "End", "Strand"}
        if not required_columns.issubset(df.columns):
            raise ValueError(f"SAF file must contain columns: {required_columns}")
    
    elif file.endswith('.gtf') or file.endswith('.gff'):
        # GTF/GFF files are tab-separated and have a standardized structure
        gtf_columns = [
            "seqname", "source", "feature", "start", "end", "score",
            "strand", "frame", "attribute"
        ]
        df = pd.read_csv(file, sep="\t", comment="#", names=gtf_columns)

        # Extract GeneID from the attribute column
        def extract_gene_id(attribute):
            """
            Extracts the GeneID from the attribute string of a GTF or GFF file entry.

            Parameters:
                attribute (str): The attribute field of a GTF/GFF file line, containing
                                key-value pairs separated by semicolons.

            Returns:
                str: The extracted GeneID. Returns "Unknown" if no GeneID is found.
            """
            for item in attribute.split(";"):
                if "gene_id" in item or "ID" in item:  # Handles both GTF and GFF formats
                    return item.split('"')[1] if '"' in item else item.split("=")[1]
            return "Unknown"

        df["GeneID"] = df["attribute"].apply(extract_gene_id)
        df = df[["GeneID", "seqname", "start", "end", "strand"]]
        df.columns = ["GeneID", "Chr", "Start", "End", "Strand"]  # Standardize column names
    
    else:
        raise ValueError(f"Unknown file type: {file}")

    return df


def read_magnipore(file : str) -> pd.DataFrame:
    """
    Reads a magnipore file and returns a pandas DataFrame.

    Parameters:
        file (str): Path to the magnipore file.

    Returns:
        pd.DataFrame: Parsed magnipore data.
    """
    return pd.read_csv(file, sep='\t')

def f_match(genomic_row : pd.Series, magni_row : pd.Series, sample : str) -> bool:
    """
    Checks if a magnipore row matches an annotation row on Chr, Strand, and position.

    Parameters:
        genomic_row (pd.Series): Row from the annotation DataFrame.
        magni_row (pd.Series): Row from the magnipore DataFrame.
        sample (str): Sample number to check in magnipore file. Must be either '1' or '2'.

    Returns:
        bool: True if the rows match, False otherwise.
    """
    if not genomic_row['Chr'] == magni_row[f'ref_{sample}']:
        return False
    if not genomic_row['Strand'] == magni_row['strand']:
        return False
    return True

def get_genes(annot: pd.DataFrame, magnipore: pd.DataFrame, sample: str) -> pd.DataFrame:
    """Find genes in the annotation that overlap with magnipore positions."""
    
    if annot.empty or magnipore.empty:
        return pd.DataFrame()
    if sample not in ["1", "2"]:
        raise ValueError("Sample must be either 1 or 2.")
    if not all(col in magnipore.columns for col in [f'pos_{sample}', f'base_{sample}', f'motif_{sample}']):
        raise ValueError(f"Magnipore file must contain columns: pos_{sample}, base_{sample}, motif_{sample}")
    
    # Create an interval index for fast lookups
    geneIntervals = pd.IntervalIndex.from_arrays(annot["Start"], annot["End"], closed="both")

    results = []  # Store results as list (avoids slow DataFrame appends)

    print("Start comparing files")
    
    for _, magni_entry in tqdm(magnipore.iterrows(), total=len(magnipore), desc="Processing Magnipore entries", unit=" entry"):
        pos = magni_entry[f'pos_{sample}']

        # Use IntervalIndex `.contains(pos)` for fast lookups
        overlapping_indices = geneIntervals.contains(pos)
        overlapping_genes = annot.loc[overlapping_indices]

        for _, genomic_entry in overlapping_genes.iterrows():
            
            if f_match(genomic_entry, magni_entry, sample):  # Assuming f_match() is an external function
                results.append({
                    "Chr": genomic_entry["Chr"],
                    "Strand": genomic_entry["Strand"],
                    "GeneID": genomic_entry["GeneID"],
                    "Start": genomic_entry["Start"],
                    "End": genomic_entry["End"],
                    "Magnipore": pos,
                    "Geneposition": pos - genomic_entry["Start"],
                    "Base" : magni_entry[f'base_{sample}'],
                    "Motif" : magni_entry[f'motif_{sample}'],
                })
    
    print(f"Processed all {len(magnipore)} lines.     ")

    # Convert list to DataFrame in one step
    return pd.DataFrame(results, columns=["Chr", "Strand", "GeneID", "Start", "End", "Magnipore", "Geneposition", "Base", "Motif"])

def main() -> None:
    args = parse()
    
    annot = read_annot(args.annot)    
    magnipore = read_magnipore(args.magnipore)
    output = get_genes(annot, magnipore, args.sample)
    output.to_csv(args.outfile, sep='\t', index=False)

if __name__ == '__main__':
    main()