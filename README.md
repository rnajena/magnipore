![](figures/magnipore_logo.png)

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-teal.svg)](https://www.gnu.org/licenses/gpl-3.0)![conda](https://img.shields.io/badge/Uses-conda-green.svg)

[![Conda package](https://anaconda.org/jannessp/magnipore/badges/version.svg)](https://anaconda.org/jannessp/magnipore) ![Conda](https://img.shields.io/conda/dn/jannessp/magnipore)
[![Conda package](https://anaconda.org/jannessp/magnipore/badges/latest_release_date.svg)](https://anaconda.org/jannessp/magnipore) [![Conda package](https://anaconda.org/jannessp/magnipore/badges/platforms.svg)](https://anaconda.org/jannessp/magnipore)

[![DOI](https://zenodo.org/badge/545997776.svg)](https://zenodo.org/badge/latestdoi/545997776)

[![Twitter Follow](https://img.shields.io/twitter/follow/Ja_Spangenberg)](https://twitter.com/Ja_Spangenberg)


>If you find a bug, please add it to the issues on GitHub with a detailed description.
---
- [Installation via Conda](#installation-via-conda)
- [Description](#description)
- [Dependencies](#dependencies)
- [Usage](#usage)
  - [Using the same reference sequence](#using-the-same-reference-sequence)
- [Output](#output)
- [Output File Description](#output-file-description)
  - [File Structure](#file-structure)
  - [Example Output](#example-output)
- [Error Codes Explanation](#error-codes-explanation)
---
# Installation via Conda

To install Magnipore we recommend to use Conda:
Magnipore is available for **linux-64 and osx-64**.

```bash
conda create -n magnipore jannessp::magnipore
conda activate magnipore
```

If you want to basecall your ONT data you also need a Guppy version from [Oxford Nanopore Technologies](https://community.nanoporetech.com).

---

# Description

Magnipore is a tool written in python3 to analyze and pair-wise compare sequencing samples from Oxford Nanopore Technologies (ONT) sequencing.

Magnipore compares two ONT samples on a signal level to find differential signals between them in single base resolution.
Such differences are caused by mutations or modifications.
Magnipore classifies these differences and provides the user with a position-wise comparison.

---

# Dependencies

Magnipore depends on/requires other tools to preprocess and analyze the data.

Conda Dependencies
- h5py >= 3.7
- biopython >= 1.80
- matplotlib >= 3.6
- numpy >= 1.21
- scipy >= 1.9
- pandas >= 1.5
- seaborn >= 0.12
- psutil >= 6.0
- pytest >= 7.1
- read5_ont >= 1.2.7
- tqdm >= 4.0

---

# Usage

If you are not using the conda package replace "magnipore" by "python3 magnipore.py".

```bash
magnipore <raw_1> <raw_2> <basecalls_1> <basecalls_2> <uncalled4_1> <uncalled4_2> <alignment> <outdir>
```

## Using the same reference sequence

Using the same reference sequence for both samples results in no reported mutations. Magnipore will only report potential modifications in this case. If you assume there are mutations between the samples, try to provide different reference sequences containing these mutations.

# Output

- Magnipore file (TSV)
  - all compared positions
  - classified into mutation and potential modification
  - with the TD score
  - with the Kullback-Leibler divergence
  - with a bayesian p-Value
- stockholm file (significant positions are marked)
- multiple plots about the data of the samples like
  - MeDAS (mean deviation average standard deviation), shows the distribution of TD scores
  - Kullback-Leibler divergence distibution
  - TD score distribution

# Output File Description

The `.magnipore` file is a tab-separated values (TSV) file containing the results of signal comparisons between two DNA samples. It provides key statistics and classifications that help distinguish between modifications and mutations.

## File Structure
Each row in the `.magnipore` file represents a single position where a comparison was made between the two samples. The table below describes the columns in the output file.

| Column Name                 | Description |
|-----------------------------|-------------|
| **strand**                  | The DNA strand (`+` or `-`) on which the comparison took place. |
| **td_score**                | Threshold distance score for the signal comparison. |
| **kl_divergence**           | Kullback-Leibler divergence for the signal comparison. |
| **bayesian_p**              | P-value from Bayesian analysis for the signal comparison. |
| **signal_type**             | Classification of the signal: `mod` (modification) or `mut` (mutation). |
| **ref_1**                   | Contig name of sample 1. |
| **pos_1**                   | Position in the contig for sample 1 (0-based). |
| **base_1**                  | Nucleotide base at the position in sample 1. |
| **motif_1**                 | DNA motif surrounding the position in sample 1. |
| **signal_mean_1**           | Mean of the signal distribution at this position in sample 1. |
| **signal_std_1**            | Standard deviation of the signal distribution at this position in sample 1. |
| **n_datapoints_1**          | Number of data points used to form the signal distribution. |
| **contained_datapoints_1**  | Number of data points within 3 standard deviations of the mean. |
| **n_segments_1**            | Number of segments from Nanopolish eventalign used in the signal distribution. |
| **contained_segments_1**    | Number of segments within 3 standard deviations of the mean. |
| **n_reads_1**               | Number of reads (coverage) used to form the signal distribution. |
| **ref_2, pos_2, base_2, motif_2, signal_mean_2, signal_std_2, n_datapoints_2, contained_datapoints_2, n_segments_2, contained_segments_2, n_reads_2** | The same fields as above, but for sample 2. |

## Example Output
```
strand  td_score    kl_divergence  bayesian_p  signal_type  ref_1        pos_1  base_1  motif_1  signal_mean_1  signal_std_1  n_datapoints_1  contained_datapoints_1  n_segments_1  contained_segments_1  n_reads_1  ref_2        pos_2  base_2  motif_2  signal_mean_2  signal_std_2  n_datapoints_2  contained_datapoints_2  n_segments_2  contained_segments_2  n_reads_2
+       1.02158245  2.43555934     0.56475101  mod          NC_000913.3  8630   A       TCAAA    -0.51221108    0.47883821    2970            2970                    56            56                    56         NC_000913.3  8630   A       TCAAA    -0.12400217    0.28117663    1500            1489                    50            48                    50
+       1.2774802   3.29518479     0.48385991  mod          NC_000913.3  49969  A       CAATC    0.45179024     0.52977556    4822            4775                    49            46                    49         NC_000913.3  49969  A       CAATC    0.97852969     0.29487824    1869            1853                    47            42                    47
```
This structured format ensures clarity and makes it easier to interpret results at a glance.****

# Error Codes Explanation

- 1: Unknown Pore Type, check --help to see which pore types are supported
- 2: Number of provided reference sequences is not equal 1 or 2
- 3: Error in multiprocessing red building
- 4: Error in multiprocessing magnipore signal comparison
- 5: Error in magniplot
