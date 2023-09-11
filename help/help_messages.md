# Magnipore

```
usage: magnipore [-h] [--guppy_bin GUPPY_BIN] [--guppy_model GUPPY_MODEL] [--guppy_device GUPPY_DEVICE] [-b1 FASTQ] [-b2 FASTQ] [-s1 TXT] [-s2 TXT] [-d] [-t THREADS] [-fr] [-mx {map-ont,splice,ava-ont}] [-mk MINIMAP2K] [--timeit] [-rna] [-r10] [-km KMER_MODEL] [-v]
                raw_data_first_sample reference_first_sample label_first_sample raw_data_sec_sample reference_sec_sample label_sec_sample working_dir

Required tools: see github https://github.com/JannesSP/magnipore

positional arguments:
  raw_data_first_sample
                        Parent directory of FAST5 files of first sample, can also be a single SLOW5 or BLOW5 file of first sample, that contains all reads, ifFASTQs are provided
  reference_first_sample
                        reference FASTA file of first sample, POSITIVE (+) or FORWARD strand, ATTENTION: can only contain a single sequence
  label_first_sample    Name of the sample or pipeline run
  raw_data_sec_sample   Parent directory of FAST5 files of second sample, can also be SLOW5 or BLOW5 file of second sample, that contains all reads, if FASTQs
                        are provided
  reference_sec_sample  reference FASTA file of second sample, POSITIVE (+) or FORWARD strand, ATTENTION: can only contain a single sequence
  label_sec_sample      Name of the sample or pipeline run
  working_dir           Path to write all output files

optional arguments:
  -h, --help            show this help message and exit
  --guppy_bin GUPPY_BIN
                        Guppy binary (default: None)
  --guppy_model GUPPY_MODEL
                        Guppy model used for basecalling (default: None)
  --guppy_device GUPPY_DEVICE
                        Use the GPU to basecall "cuda:0" to use the GPU with ID 0 (default: cuda:0)
  -b1 FASTQ, --basecalls_first_sample FASTQ
                        Path to existing basecalls of first sample. Basecalls must be in one single file. (default: None)
  -b2 FASTQ, --basecalls_sec_sample FASTQ
                        Path to existing basecalls of second sample. Basecalls must be in one single file. (default: None)
  -s1 TXT, --sequencing_summary_first_sample TXT
                        Use, when sequencing summary is not next to your FASTQ file. Path to existing sequencing summary file of second sample. (default:
                        None)
  -s2 TXT, --sequencing_summary_sec_sample TXT
                        Use, when sequencing summary is not next to your FASTQ file. Path to existing sequencing summary file of first sample. (default: None)
  -d, --calculate_data_density
                        Will calculate data density after building the models. Will increase runtime! (default: False)
  -t THREADS, --threads THREADS
                        Number of threads to use (default: 1)
  -fr, --force_rebuild  Run commands regardless if files are already present (default: False)
  -mx {map-ont,splice,ava-ont}, --minimap2x {map-ont,splice,ava-ont}
                        -x parameter for minimap2 (default: map-ont)
  -mk MINIMAP2K, --minimap2k MINIMAP2K
                        -k parameter for minimap2 (default: 14)
  --timeit              Measure and print time used by submodules (default: False)
  -rna, --rna           Use when data is rna (default: False)
  -r10, --r10           Use when data is from R10.4.1 flowcell (default: False)
  -km KMER_MODEL, --kmer_model KMER_MODEL
                        custom kmer model file for f5c eventalign (default: None)
  -v, --version         show program's version number and exit
```

# Nanosherlock

```
usage: nanosherlock [-h] [--guppy_bin GUPPY_BIN] [--guppy_model GUPPY_MODEL] [--guppy_device GUPPY_DEVICE] [--basecalls FASTQ] [--sequencing_summary TXT] [--calculate_data_density] [-t THREADS] [-fr] [-mx {map-ont,splice,ava-ont}] [-mk MINIMAP2K] [--timeit] [--max_lines MAX_LINES] [-rna] [-r10] [-km KMER_MODEL] [-e ERROR_PREFIX] [-v] raw_data reference working_dir sample_label

Required tools: see github https://github.com/JannesSP/magnipore

positional arguments:
  raw_data              Parent directory of FAST5 files, can also be a direct path to a single SLOW5 or BLOW5 file, that contains all reads, if FASTQs are provided
  reference             reference FASTA file, POSITIVE (+) or FORWARD strand
  working_dir           Path to write all output files
  sample_label          Name of the sample or pipeline run

optional arguments:
  -h, --help            show this help message and exit
  --guppy_bin GUPPY_BIN
                        Guppy binary (default: None)
  --guppy_model GUPPY_MODEL
                        Guppy model for basecalling (default: None)
  --guppy_device GUPPY_DEVICE
                        Use the gpu to basecall with cuda:0 (default: None)
  --basecalls FASTQ     Path to existing basecalls. Basecalls must be in one single file. (default: None)
  --sequencing_summary TXT
                        Use, when sequencing summary is not next to your FASTQ file. Path to existing sequencing summary file of sample. (default: None)
  --calculate_data_density
                        Will calculate data density after building the models. Will increase runtime! (default: False)
  -t THREADS, --threads THREADS
                        Number of threads to use (default: None)
  -fr, --force_rebuild  Run commands regardless if files are already present (default: False)
  -mx {map-ont,splice,ava-ont}, --minimap2x {map-ont,splice,ava-ont}
                        -x parameter for minimap2 (default: map-ont)
  -mk MINIMAP2K, --minimap2k MINIMAP2K
                        -k parameter for minimap2 (default: 14)
  --timeit              Measure and print time used by submodules (default: False)
  --max_lines MAX_LINES
                        Only process first given number of lines from segmentation eventalign (default: None)
  -rna                  Use when data is rna (default: False)
  -r10                  Use when data is from R10.4.1 flowcell (default: False)
  -km KMER_MODEL, --kmer_model KMER_MODEL
                        custom kmer model file for f5c eventalign (default: None)
  -e ERROR_PREFIX, --error_prefix ERROR_PREFIX
  -v, --version         show program's version number and exit
```

# Magniplot

```
usage: magniplot [-h] [-f FONTSIZE] [-t THREADS] [-nl NUM_LINES] [-ml MAX_LINES] [-c COVERAGE] [-s SEED] [-v]
                 magnipore outdir label_first_sample label_sec_sample

Creating plots (MeDAS, etc) for a given .magnipore file. If the number of entries in the .magnipore file exceeds `max_lines`, this script will randomly sample from the .magnipore file to create the plot.

positional arguments:
  magnipore             Magnipore-style output
  outdir                Path to write plots.
  label_first_sample    Name of the sample or pipeline run
  label_sec_sample      Name of the sample or pipeline run

optional arguments:
  -h, --help            show this help message and exit
  -f FONTSIZE, --fontsize FONTSIZE
                        Fontsize for plots (default: 18)
  -t THREADS, --threads THREADS
                        Number of processes to use to create plots (default: 1)
  -nl NUM_LINES, --num_lines NUM_LINES
                        Providing the number of lines in file speeds up the process. (default: None)
  -ml MAX_LINES, --max_lines MAX_LINES
                        Plot max this number of entries. Do not set this score too high, as it increases runtime and memory usage. If you have data with a low coverage, many entries/lines in the .magnipore file could be NANs. These are filtered out. Increase this number have more get more entries with data. (default: 1500000)
  -c COVERAGE, --coverage COVERAGE
                        Coverage cutoff threshold for the plots. (default: 10)
  -s SEED, --seed SEED  Set a random seed to reproduce the same image. (default: None)
  -v, --version         show program's version number and exit
```

# Magnicheck

```
usage: magnicheck [-h] [--coverage COVERAGE] [--valid_sep VALID_SEP] [--pore {r9,r10}] [-v] magnipore magnipore_poscol refid eval_csv refcol poscol outfile

A small script that compares the Magnipore output file with a given validationset in form of a table in a CSV file. Your table should contain the reference
id/name, that was also used during the Magnipore analysis, and the validation positions on the given reference. The script will then check, if Magnipore found
significant positions in a kmer range of the positions in the validation table (CSV). The kmer range depends on the used pore during sequencing. You can
specify the used pore with the --pore parameter. You should also think about the coverage threshold. This script will by default filter out positions from the
Magnipore output, where at least one sample has a coverage less than 10 reads.

positional arguments:
  magnipore             Magnipore file with called differential modifications.
  magnipore_poscol      Which position to validate from the Magnipore output.
  refid                 Reference id or name - This id must match between the magnipore file and the your validation table (CSV).
  eval_csv              CSV file containing the validation table. The table contains the reference id and position of validated (ground truth) modifications.
  refcol                Column containing the reference ids.
  poscol                Column containing the validated positions.
  outfile               Name of the output file.

optional arguments:
  -h, --help            show this help message and exit
  --coverage COVERAGE   Coverage filter to apply to the Magnipore output. (default: 10)
  --valid_sep VALID_SEP
                        Separation character in your CSV file. (default: ,)
  --pore {r9,r10}       Which pore you used during sequencing. (default: r9)
  -v, --version         show program's version number and exit
```

# Magnifilter

```
usage: magnifilter [-h] [-c COVERAGE] [-v] magnipore

A small script to filter the .magnipore output for a given coverage threshold. Please provide the path to the .magnipore file and a coverage threshold. All
compared positions, where at least one sample has a coverage below your given threshold will be filtered out. The remaining positions are written to a file
with the name of the provided magnipore file and the suffix "_c10": e.g. <given_magnipore_file>_c10.magnipore

positional arguments:
  magnipore             .magnipore file to filter for a given coverage threshold

optional arguments:
  -h, --help            show this help message and exit
  -c COVERAGE, --coverage COVERAGE
                        Coverage threshold to filter for. Results, where at least one sample has a coverage below the given threshold are filtered out. Results, where both samples have a coverage equal or higher than the threshold remain. (default: 10)
  -v, --version         show program's version number and exit
```