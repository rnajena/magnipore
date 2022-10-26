# Magnipore

... still in development, meaning bugs will occur and not everything is tested.
This is a first draft of the magnipore repository for public use.
This tools will be published soon.
If you find a bug please add it to the issues on github with a detailed desciption. :)

## Description

Magnipore is used to compare two ONT samples with each other on the signal level to find differential signals between these samples.
These differences occur majorly due to mutations and modifications.
Magnipore tries to classify these differences and provides the user with a positionwise comparison on reference base resolution together with quality values like a bayesian p-value or the coverage.

## Usage

### Conda

>Magnipore conda coming soon

For now create a conda environment using the [conda.yml](conda_env/conda.yml).
If you want to basecall your ONT data you also need a Guppy version.

Required tools in environment: guppy, minimap2, ococo, medaka, nanopolish, h5py, samtools, scipy and mafft

#### Simplest use case:

    magnipore.py path_to_fast5_first_sample path_to_reference_first_sample first_sample_label path_to_fast5_sec_sample path_to_reference_sec_sample sec_sample_label working_dir --guppy_bin PATH --guppy_model PATH

#### positional arguments:
- path_to_fast5_first_sample : FAST5 directory of first sample
- path_to_reference_first_sample : reference FASTA file of first sample
- first_sample_label : Name of the sample or pipeline run
- path_to_fast5_sec_sample : FAST5 file of second sample
- path_to_reference_sec_sample : reference FASTA file of second sample
- sec_sample_label : Name of the sample or pipeline run
- working_dir : Path to write all output files

#### required arguments:
use either the basecalling arguments or provide basecalls
- basecalling arguments:
    - guppy_bin : Path to guppy binary
    - guppy_model : Path to guppy model used for basecalling
    - (optional) guppy_device : Device used for basecalling (cpu or gpu cuda:0)
- provided basecalls (FASTQ)
    - path_to_first_basecalls : Path
    - path_to_sec_basecalls : Path

For optional arguments see magnipore.py --help. Includes small number of mapping parameters and the option to skip basecalling.

## Output

The .magnipore file is a tsv containing the following columns.

- strand : on which strand the comparison took place
- td_score : threshold distance score for the signal comparison
- kl_divergence : kullback leibler divergence for the signal comparison
- bayesian_p : p-value for the signal comparison
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