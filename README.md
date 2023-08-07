# ![](figures/magnipore_logo.png)

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-teal.svg)](https://www.gnu.org/licenses/gpl-3.0)![conda](https://img.shields.io/badge/Uses-conda-green.svg)

[![Conda package](https://anaconda.org/jannessp/magnipore/badges/version.svg)](https://anaconda.org/jannessp/magnipore) ![Conda](https://img.shields.io/conda/dn/jannessp/magnipore)
[![Conda package](https://anaconda.org/jannessp/magnipore/badges/latest_release_date.svg)](https://anaconda.org/jannessp/magnipore) [![Conda package](https://anaconda.org/jannessp/magnipore/badges/platforms.svg)](https://anaconda.org/jannessp/magnipore)

[![DOI](https://zenodo.org/badge/545997776.svg)](https://zenodo.org/badge/latestdoi/545997776)

[![Twitter Follow](https://img.shields.io/twitter/follow/Ja_Spangenberg)](https://twitter.com/Ja_Spangenberg)


>If you find a bug, please add it to the issues on GitHub with a detailed description.
___
1.  [Installation via Conda](#installation-via-conda)
2.  [Description](#description)
3.  [Dependencies](#dependencies)
4.  [Workflow](#workflow)
5.  [Usage](#usage)
6.  [Output File Description](#output-file-description)
7.  [Error Codes Explanation](#error-codes-explanation)
___
## Installation via Conda

To install Magnipore we recommend to use Conda:
Magnipore is available for **linux-64 and osx-64**.

```bash
conda install mamba
mamba create -n magnipore -c jannessp magnipore
conda activate magnipore
```

*Alternatively you can create a conda environment using the [conda_env.yml](conda.recipe/conda_env.yml) and mamba.*
```bash
conda install mamba
mamba env create -f conda/conda_env.yml
conda activate magnipore
git clone https://github.com/JannesSP/magnipore.git
```

If you want to basecall your ONT data you also need a Guppy version from [Oxford Nanopore Technologies](https://community.nanoporetech.com).

---

## Description

Magnipore is a tool written in python3 to analyze and pair-wise compare sequencing samples from Oxford Nanopore Technologies (ONT) sequencing.

Magnipore compares two ONT samples on a signal level to find differential signals between them in single base resolution.
Such differences are caused by mutations or modifications.
Magnipore classifies these differences and provides the user with a position-wise comparison.

---

## Dependencies

Magnipore depends on/requires other tools to preprocess and analyze the data.

<details><summary>Click here to see dependencies</summary>

Conda Dependencies
- python (>=3.8,<3.11)
- h5py>=3.7
- biopython>=1.80
- mafft>=7.508
- matplotlib>=3.6.2
- numpy>=1.23
- scipy>=1.9
- minimap2>=2.24
- pandas>=1.5
- seaborn>=0.12
- psutil>=5.9
- hdf5plugin>=3.3.1
- ont_vbz_hdf_plugin>=1.0.1
- pytest>=7.1
- gzip>=1.12
- read5>=1.1.6
- f5c>=1.2
- read5>=1.2.0

</details>

---

## Workflow

### Input

For each sample in the comparison, Magnipore takes:
- (FASTA) exactly ONE reference sequence
- (FAST5) the raw sequencing data from ONT
- (optinal FASTQ) optionally basecalls, if you do not have the guppy binary or do not want to basecall the raw ONT data (again).

### Output

- Magnipore file (TSV)
  - all compared positions
  - classified into mutation and potential modification
  - with the TD score
  - with the Kullback-Leibler divergence
  - with a bayesian p-Value
- reference sequence alignment file
- stockholm file (significant positions are marked)
- multiple plots about the data of the samples like

---

## Usage

If you are not using the conda package replace "magnipore" by "python3 magnipore.py".

### Without basecalling:

<details><summary>Click here to see command:</summary>

```bash
magnipore raw_data_first_sample reference_first_sample label_first_sample raw_data_sec_sample reference_sec_sample label_sec_sample working_dir --basecalls_first_sample basecalls_first_sample --basecalls_sec_sample basecalls_sec_sample
```
</details>

### With basecalling

<details><summary>Click here to see command:</summary>

```bash
magnipore raw_data_first_sample reference_first_sample label_first_sample raw_data_sec_sample reference_sec_sample label_sec_sample working_dir --guppy_bin PATH --guppy_model PATH
```
</details>

### Using a single sequencing run with demultiplexed FASTQs

<details><summary>Click here to see command:</summary>

- basecalls_first_sample/basecalls_sec_sample containing the demultiplexed FASTQs
 - *label_first_sample.fastq* contains only those reads of the first condition
 - *label_sec_sample.fastq* contains only those reads of the second condition
- be sure that the *sequencing_summary.txt* is next to your FASTQ files, otherwise provide them using
 - -s1, --sequencing_summary_first_sample
 - -s2, --sequencing_summary_sec_sample

```bash
magnipore --basecalls_first_sample basecalls_first_sample --basecalls_sec_sample basecalls_sec_sample raw_data_first_sample reference_first_sample label_first_sample raw_data_sec_sample reference_sec_sample label_sec_sample working_dir
```
</details>

### Using the same reference sequence

Using the same reference sequence for both samples results in no reported mutations. Magnipore will only report potential modifications in this case. If you assume there are mutations between the samples, try to provide different reference sequences containing these mutations.

### Help Message

<details><summary>Click here to see help message:</summary>

```bash
usage: Magnipore [-h] [--guppy_bin GUPPY_BIN] [--guppy_model GUPPY_MODEL] [--guppy_device GUPPY_DEVICE] [-b1 FASTQ] [-b2 FASTQ] [-s1 TXT] [-s2 TXT] [-d] [-t THREADS] [-fr]
                 [-mx {map-ont,splice,ava-ont}] [-mk MINIMAP2K] [--timeit] [--rna] [-v]
                 raw_data_first_sample reference_first_sample label_first_sample raw_data_sec_sample reference_sec_sample label_sec_sample working_dir

Required tools: see github https://github.com/JannesSP/magnipore

positional arguments:
  raw_data_first_sample
                        Parent directory of FAST5 files of first sample, can also be a single SLOW5 or BLOW5 file of first sample, that contains all reads, if FASTQs are
                        provided
  reference_first_sample
                        reference FASTA file of first sample, POSITIVE (+) or FORWARD strand, ATTENTION: can only contain a single sequence
  label_first_sample    Name of the sample or pipeline run
  raw_data_sec_sample   Parent directory of FAST5 files of second sample, can also be SLOW5 or BLOW5 file of second sample, that contains all reads, if FASTQs are provided
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
                        Use, when sequencing summary is not next to your FASTQ file. Path to existing sequencing summary file of second sample. (default: None)
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
  -rna                  Use when data is rna (default: False)
  -r10                  Use when data is from R10.4.1 flowcell (default: False)
  -km KMER_MODEL, --kmer_model KMER_MODEL
                        custom kmer model file for f5c eventalign (default: None)
  -v, --version         show program's version number and exit
```
</details>

#### required arguments:
use either the basecalling arguments or provide basecalls
- basecalling arguments:
    - guppy_bin : Path to guppy binary
    - guppy_model : Path to guppy model used for basecalling
    - (optional) guppy_device : Device used for basecalling (cpu or gpu cuda:0)
- provided basecalls (FASTQ)
    - basecalls_first_sample : Path
    - basecalls_sec_sample : Path

For optional arguments see magnipore.py --help. Includes small number of mapping parameters and the option to skip basecalling.

## Output File Description

<details><summary>Click here to see overview:</summary>
The .magnipore file is a TSV containing the following columns.

- strand : on which strand the comparison took place
- td_score : threshold distance score for the signal comparison
- kl_divergence : kullback leibler divergence for the signal comparison
- bayesian_p : p-value for the signal comparison
- signal_type : classification into "mod" for modification and "mut" for mutation
- ref_1 : contig name of sample 1
- pos_1 : position in contig of sample 1
- base_1 : base at the position of sample 1
- motif_1 : motif around the base at the position of sample 1
- signal_mean_1 : mean of the signal distribution at the position of sample 1
- signal_std_1 : standard deviation of the signal distribution at the position of sample 1
- n_datapoints_1 : number of data points that formed the signal distribution
- contained_datapoints_1 : number of data points withtin 3 standard deviations around the mean
- n_segments_1 : number of segments from nanopolish eventalign that formed the signal distribution
- contained_segments_1 : number of segments within 3 standard deviations around the mean
- n_reads_1 : number of reads (coverage) that formed the signal distribution

same for second sample:
- ref_2, pos_2, base_2, motif_2, signal_mean_2, signal_std_2, n_datapoints_2, contained_datapoints_2, n_segments_2, contained_segments_2, n_reads_2
</details>

## Error Codes Explanation

<details><summary>Click here to see error codes:</summary>

- 11: Concatenating both reference files failed
- 12: Building mafft alignment failed
- 13: Running nanosherlock of the first sample failed
- 14: Running nanosherlock of the second sample failed
- 15: Number of provided reference sequences is not equal 1 or 2
---
Errors of first sample:
- 119: Cannot basecall .slow5/.blow5 with guppy
- 120: Could not find raw data or unknown file format
- 121: Guppy basecalling failed
- 122: minimap2 mapping failed
- 123: Samtools indexing failed
- 124: f5c index failed
- 125: f5c eventalign failed
- 126: Could not find provided fastq files
- 127: Could not find provided sequencing summary file
---
Errors of second sample
- 219: Cannot basecall .slow5/.blow5 with guppy
- 220: Could not find raw data or unknown file format
- 221: Guppy basecalling failed
- 222: minimap2 mapping failed
- 223: Samtools indexing failed
- 224: f5c index failed
- 225: f5c eventalign failed
- 226: Could not find provided fastq files
- 227: Could not find provided sequencing summary file

### If Subscript Nanosherlock is Executed Separately

The -e parameter of nanosherlock specifies the leading number of the error code. Default is 0.
- 019: Cannot basecall .slow5/.blow5 with guppy
- 020: Could not find raw data or unknown file format 
- 021: Guppy basecalling failed
- 022: minimap2 mapping failed
- 023: Samtools indexing failed
- 024: f5c index failed
- 025: f5c eventalign failed
- 026: Could not find provided fastq files
- 027: Could not find provided sequencing summary file
  </details>