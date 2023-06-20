#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import datetime
import os
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from pathlib import Path
import h5py
import numpy as np
from scipy import stats
from Bio import SeqIO
import pandas as pd
from time import perf_counter_ns
import os
import psutil

from magnipore.__init__ import __version__
from magnipore.Logger import Logger
from magnipore.Helper import ANSI, REDENCODER, STRANDDECODER, sizeof_fmt
import magnipore.OnlineMeanVar as omv

PROCESS = psutil.Process(os.getpid())
LOGGER : Logger = None
TIMEIT = False
ERROR_PREFIX : str = '0'

def initLogger(file = None) -> None:
    global LOGGER
    LOGGER = Logger(file)

def readFast5(fast5_path : str, sequencing_summary : str = None) -> dict:
    '''
    Returns a dictionary mapping the readid to the corresponding fast5 file containing the read.
    '''
    LOGGER.printLog(f'Loading readids from fast5 files in {fast5_path} ... ')
    readid2fast5 = {}

    if sequencing_summary is None:
        # recursively loop over the fast5 files in the given path
        fast5list = Path(fast5_path).rglob('*.fast5')

        for fidx, fast5 in enumerate(fast5list):
            LOGGER.printLog(f'Indexing fast5 file {fidx + 1}\r', newline_after=False)

            # because path is object not string
            fast5 = str(fast5)
            fast5_h5 = h5py.File(fast5, 'r')

            for readid in fast5_h5:
                readid = readid.split('read_')[1]
                readid2fast5[readid] = fast5
            fast5_h5.close()

    else:
        with open(sequencing_summary, 'r') as seqsum:
            seqsum.readline()

            for ridx, line in enumerate(seqsum):
                if (ridx + 1) % 1000 == 0:
                    LOGGER.printLog(f'Indexing read {ridx + 1}\r', newline_after=False)
            
                filename, read_id = line.strip().split('\t')[:2]
                readid2fast5[read_id] = os.path.join(fast5_path, filename)
    
    LOGGER.printLog(f'Indexed {ridx + 1} reads')
    return readid2fast5

def readNanoSum(nanoSum_path : str):
    '''
    Read summary file with read_indices and read_ids.\n
    Return a dictionary containing the read index from nanopolish as keys and the read id as values.
    '''
    read_index2ID = {}
    # LOGGER.printLog(f'Start loading ids from {nanoSum_path} ...')
    with open(nanoSum_path, 'r') as summary:
        next(summary) # skip header

        for lidx, line in enumerate(summary):
            if (lidx + 1) % 1000 == 0:
                LOGGER.printLog(f'Loading id {lidx + 1}\r', newline_after=False)
            r_index, r_ID = line.strip().split()[:2]
            read_index2ID[int(r_index)] = r_ID

    LOGGER.printLog(f'Found {len(read_index2ID)} readids')
    return read_index2ID

def parse() -> Namespace:

    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description='Required tools: see github https://github.com/JannesSP/magnipore',
        prog='Nanosherlock',
        ) 
    
    parser.add_argument('path_to_fast5', type = str, help='FAST5 file')
    parser.add_argument('path_to_reference', type = str, help='reference FASTA file, POSITIVE (+) or FORWARD strand')
    parser.add_argument('working_dir', type = str, help='Path to write all output files')
    parser.add_argument('sample_label', type = str, help='Name of the sample or pipeline run')
    
    parser.add_argument('--guppy_bin', type = str, default = None, help='Guppy binary')
    parser.add_argument('--guppy_model', type = str, default = None, help='Guppy model for basecalling')
    parser.add_argument('--guppy_device', type=str, default=None, help='Use the gpu to basecall with cuda:0')

    parser.add_argument('--path_to_basecalls', metavar='FASTQ_DIR', default = None, type = str, help = 'Path to existing basecalls and sequencing summary file. Basecalls must be in one single file with the name <sample_label>.fastq')
    parser.add_argument('--calculate_data_density', action = 'store_true', default = False, help = 'Will calculate data density after building the models. Will increase runtime!')
    parser.add_argument('-t', '--threads', type=int, help='Number of threads to use')
    parser.add_argument('-f5', '--fast5_out', action = 'store_true', help='Guppy generates FAST5 output (workspace folder)')
    parser.add_argument('-fr', '--force_rebuild', action = 'store_true', help='Run commands regardless if files are already present')
    parser.add_argument('-mx', '--minimap2x', default = 'splice', choices = ['map-ont', 'splice', 'ava-ont'], help = '-x parameter for minimap2')
    parser.add_argument('-mk', '--minimap2k', default = 14, type = int, help = '-k parameter for minimap2')
    parser.add_argument('--timeit', default = False, action = 'store_true', help = 'Measure and print time used by submodules')
    parser.add_argument('--max_lines', default=None, type=int, help='Only process first given number of lines from nanopolish eventalign')

    parser.add_argument('-e', '--error_prefix', type=str, default=None)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s' + f' {__version__}')

    return parser.parse_args()

def guppy_basecalling(guppy_bin : str, guppy_model : str,  guppy_device : str, path_to_fast5 : str, working_dir : str, sample_label : str, fast5_out : bool, force_rebuild : bool, sequencing_summary : str) -> str:

    basecalls_path = os.path.join(working_dir, 'basecalls', sample_label)
    basecalls = os.path.join(basecalls_path, sample_label + ".fastq")

    if sequencing_summary is None:
        sequencing_summary = os.path.join(basecalls_path, 'sequencing_summary.txt')
    else:
        assert os.path.exists(sequencing_summary)

    if os.path.exists(sequencing_summary) and os.path.exists(basecalls) and not force_rebuild:
        LOGGER.printLog(f'basecalls and sequencing summary already exist:\n-\t{basecalls}\n-\t{sequencing_summary}')
        return basecalls, sequencing_summary, force_rebuild
    
    else:
        force_rebuild = True

    if not os.path.exists(basecalls_path):
        # LOGGER.printLog(f'Creating sub-directory {basecalls_path}')
        os.makedirs(basecalls_path)
    
    # LOGGER.printLog(f'Run guppy basecalling on {path_to_fast5}')
    command = f'{guppy_bin} -i {path_to_fast5} -s {basecalls_path} --disable_pings --disable_qscore_filtering --calib_detect -c {guppy_model}'

    if fast5_out:
        command += ' --fast5_out'
    
    if guppy_device is not None:
        command += f' --device {guppy_device}'

    LOGGER.printLog(f'guppy command: {ANSI.GREEN}{command}{ANSI.END}')
    
    if TIMEIT:
        start = perf_counter_ns()
    ret = os.system(command)
    if ret != 0:
        LOGGER.error(f'Error in guppy basecalling with error code {ret}', error_type=ERROR_PREFIX+'21')
    if TIMEIT:
        end = perf_counter_ns()
        LOGGER.printLog(f'{ANSI.YELLOW}TIMED: Guppy basecalling took {pd.to_timedelta(end-start)}, {end - start} nanoseconds{ANSI.END}')

    # concatenate all fastq files into one fastq file
    os.system(f'cat {os.path.join(basecalls_path, "*.fastq")} > {os.path.join(basecalls_path, "*.tmpfastq")}')
    os.system(f'rm {os.path.join(basecalls_path, "*.fastq")}')
    os.system(f'mv {os.path.join(basecalls_path, "*.tmpfastq")} {basecalls}')
        
    # move basecalled fast5s to a seperate directory
    if fast5_out:
        basecalled_fast5_path = os.path.join(working_dir, "basecalled_fast5", sample_label)
        # LOGGER.printLog(f'Move basecalled fast5 files to {basecalled_fast5_path}')
        os.makedirs(basecalled_fast5_path)
        os.system(f'mv {os.path.join(basecalls_path, "workspace")} {basecalled_fast5_path}')
        
    return basecalls, sequencing_summary, force_rebuild

def minimap(path_to_reference : str, path_to_basecalls : str, working_dir : str, sample_label : str, threads : int, force_rebuild : bool, mx : str, mk : int) -> str:
    
    minimap_path = os.path.join(working_dir, 'minimap2', sample_label)
    bam_path = os.path.join(minimap_path, sample_label) + '.bam'

    if os.path.exists(bam_path) and not force_rebuild:
        LOGGER.printLog(f'Bam alignment/mapping file already exists\n-\t{bam_path}')
        return bam_path, force_rebuild
    else:
        force_rebuild = True
    
    if not os.path.exists(minimap_path):
        # LOGGER.printLog(f'Creating sub-directory {minimap_path}')
        os.makedirs(minimap_path)
        
    # LOGGER.printLog(f'Run minimiap2 on basecalls: {path_to_basecalls} and reference: {path_to_reference}')
    # splice because of sarscov genome/transciptome functionality
    command = f'minimap2 -a -x {mx} -k{mk} -t {threads} {path_to_reference} {path_to_basecalls} | samtools view -hbF4 | samtools sort > {bam_path}'
    LOGGER.printLog(f'minimap2 command: {ANSI.GREEN}{command}{ANSI.END}')
    
    if TIMEIT:
        start = perf_counter_ns()
    ret = os.system(command)
    if ret != 0:
        LOGGER.error(f'Error in minimap2 run with error code {ret}', error_type=ERROR_PREFIX+'22')
    if TIMEIT:
        end = perf_counter_ns()
        LOGGER.printLog(f'{ANSI.YELLOW}TIMED: minimap2 took {pd.to_timedelta(end-start)}, {end - start} nanoseconds{ANSI.END}')

    return bam_path, force_rebuild

def nanopolish(path_to_fast5 : str, path_to_sequencing_summary : str, path_to_basecalls : str, path_to_reference : str, alignment_bam : str, working_dir : str, sample_label : str, threads : int, force_rebuild : bool) -> str:
    
    nanopolish_path = os.path.join(working_dir, 'nanopolish', sample_label)
    summary_csv = os.path.join(nanopolish_path, "eventalign_summary.csv")
    result_csv = os.path.join(nanopolish_path, "eventalign_result.csv")

    if os.path.exists(summary_csv) and os.path.exists(result_csv) and not force_rebuild:
        LOGGER.printLog(f'nanopolish files already exist:\n-\t{summary_csv}\n-\t{result_csv}')
        return summary_csv, result_csv, force_rebuild

    else:
        force_rebuild = True

    if not os.path.exists(nanopolish_path):
        # LOGGER.printLog(f'Creating sub-directory {nanopolish_path}')
        os.makedirs(nanopolish_path)

    if not os.path.exists(alignment_bam + '.bai') or force_rebuild:

        # samtools indexing
        # LOGGER.printLog(f'Indexing {alignment_bam} with samtools')
        command = f'samtools index {alignment_bam}'
        LOGGER.printLog(f'samtools index command: {ANSI.GREEN}{command}{ANSI.END}')
    
        if TIMEIT:
            start = perf_counter_ns()
        ret = os.system(command)
        if ret != 0:
            LOGGER.error(f'Error in samtools indexing for nanopolish with error code {ret}', error_type=ERROR_PREFIX+'23')
        if TIMEIT:
            end = perf_counter_ns()
            LOGGER.printLog(f'{ANSI.YELLOW}TIMED: samtools indexing took {pd.to_timedelta(end-start)}, {end - start} nanoseconds{ANSI.END}')

    else:
        LOGGER.printLog(f"Using existing {os.path.exists(alignment_bam + '.bai')}")

    if not os.path.exists(path_to_basecalls + '.index') or force_rebuild:

        # LOGGER.printLog(f'Run nanopolish indexing')
        command = f'nanopolish index -d {path_to_fast5} -s {path_to_sequencing_summary} {path_to_basecalls}'
        LOGGER.printLog(f'nanopolish indexing command: {ANSI.GREEN}{command}{ANSI.END}')

        if TIMEIT:
            start = perf_counter_ns()
        ret = os.system(command)
        if ret != 0:
            LOGGER.error(f'Error in nanopolish indexing with error code {ret}', error_type=ERROR_PREFIX+'24')
        if TIMEIT:
            end = perf_counter_ns()
            LOGGER.printLog(f'{ANSI.YELLOW}TIMED: nanopolish indexing took {pd.to_timedelta(end-start)}, {end - start} nanoseconds{ANSI.END}')

    else:
        LOGGER.printLog(f"Using existing {os.path.exists(path_to_basecalls + '.index')}")

    
    # LOGGER.printLog(f'Run nanopolish eventalign')
    command = f'nanopolish eventalign --reads {path_to_basecalls} --bam {alignment_bam} --genome {path_to_reference} --summary={summary_csv} --scale-events --signal-index -t {threads} > {result_csv}'
    LOGGER.printLog(f'Nanopolish eventalign command: {ANSI.GREEN}{command}{ANSI.END}')
    
    if TIMEIT:
        start = perf_counter_ns()
    ret = os.system(command)
    if ret != 0:
        LOGGER.error(f'Error in nanopolish eventalign with error code {ret}', error_type=ERROR_PREFIX+'25')
    if TIMEIT:
        end = perf_counter_ns()
        LOGGER.printLog(f'{ANSI.YELLOW}TIMED: nanopolish eventalign took {pd.to_timedelta(end-start)}, {end - start} nanoseconds{ANSI.END}')

    return summary_csv, result_csv, force_rebuild

def read_line(line, event : dict):
    
    line = line.strip().split()

    event['contig'] = line[0]
    # add 2 to collect the signals for the middle base of the 5mer event
    # position describes the first base of the event kmer in the reference
    event['position'] = int(line[1]) + 2
    event['ref_kmer'] = line[2]
    event['read_index'] = int(line[3])
    # line_dic['strand'] = line[4]
    # line_dic['event_index'] = int(line[5])
    # line_dic['event_mean'] = float(line[6])
    # line_dic['event_std'] = float(line[7])
    # line_dic['event_len'] = float(line[8])
    event['model_kmer'] = line[9]
    # line_dic['model_mean'] = float(line[10])
    # line_dic['model_std'] = float(line[11])
    # line_dic['standardized_level'] = float(line[12])
    event['start_idx'] = int(line[13])
    event['end_idx'] = int(line[14])

def getZNormSignal(fast5_read : h5py.Group, mode : str = 'median') -> np.ndarray:
    '''
    Z normalise the nanopore signal by shift and scale.
    
    Parameters
    ----------
    fast5_read : h5py.Group
        Contain the FAST5 read
    mode : str
        mode to use for Z normalisation
        median: normalises median -> 0 and mad -> 1, not affected by outliers but takes longer
        mean: normalises mean -> 0 and stdev -> 1, faster but affected by outliers

    Returns
    -------
    normalised signal as numpy array
    '''
    dacs = fast5_read['Raw/Signal'][:]
    r_offset = fast5_read['channel_id'].attrs['offset']
    r_range = fast5_read['channel_id'].attrs['range']
    r_digit = fast5_read['channel_id'].attrs['digitisation']
    
    # https://github.com/nanoporetech/taiyaki/blob/master/docs/FILE_FORMATS.md
    # Dacs is the pA signal transformed to integers
    # using offset, range and digit to transform it back
    read_pA = (dacs + r_offset) * (r_range / r_digit)
    
    if mode == 'median':
        # https://github.com/nanoporetech/taiyaki/blob/9672c66d654ad8101c4faf481646493082e090c8/taiyaki/maths.py#L8
        # shift and scale are calculated as described by nanoporetech, see the links above
        r_shift = np.median(read_pA)
        r_scale = 1.4826 * np.median(np.abs(read_pA - r_shift))
        return (read_pA - r_shift) / r_scale
    elif mode == 'mean':
        return (read_pA - read_pA.mean()) / read_pA.std()

def aggregate_events(nanopolish_result_csv : str, nanopolish_summary_csv: str, path_to_fast5 : str, path_to_reference : str, working_dir : str, sample_label : str, force_rebuild : bool, sequencing_summary : str, calculate_data_density : bool, max_lines : int) -> str:
    
    if not os.path.exists(os.path.join(working_dir, 'magnipore', sample_label)):
        os.makedirs(os.path.join(working_dir, 'magnipore', sample_label))

    # write distributions into reference event distribution file
    red_file = os.path.join(working_dir, 'magnipore', sample_label, f'{sample_label}.red')

    if os.path.exists(red_file) and not force_rebuild:

        return red_file, force_rebuild

    else:

        force_rebuild = True

    readid2fast5 = readFast5(path_to_fast5, sequencing_summary)
    nano2readid = readNanoSum(nanopolish_summary_csv)
    ref_seqs = SeqIO.parse(open(path_to_reference), 'fasta')
    red_dict, omvs = createREDDict(ref_seqs)
    
    LOGGER.printLog(f'Found sequences: {list(red_dict.keys())}. Start getting event signal distributions per reference position and write red file for {sample_label}')

    if TIMEIT:
        start = perf_counter_ns()

    buildModels(red_dict, omvs, nano2readid, readid2fast5, nanopolish_result_csv, calculate_data_density, max_lines)
    del omvs # actively release memory

    if TIMEIT:
        end = perf_counter_ns()
        LOGGER.printLog(f'{ANSI.YELLOW}TIMED: Building distribution models took {pd.to_timedelta(end-start)}, {end - start} nanoseconds{ANSI.END}')

    LOGGER.printLog('Writing output files')
    writeOutput(red_file, red_dict)
        
    return red_file, force_rebuild

def createREDDict(fasta) -> tuple:
    '''
    Parameters
    ----------
    fasta : SeqIO.parse Generator

    Returns
    -------
    red_dict : dict
        stores the data that will be written to the RED files
    omvs : dict
        stores OnlineMeanVar objects to calculate mean and stdev using bayes
    '''

    # get all sequences from the reference file
    # sequences are stored as {reference: {pos: {strand: {base: (A|C|G|T), signal: []}}}}
    red_dict = {}
    omvs = {}
    omv_size = 30
    LOGGER.printLog(f'Using {omv_size} datapoints for OMV initialization.')
    for seq in fasta:
        
        # default values for mean, std etc is 0, default for expected model density is nan
        seq_size = len(seq.seq)
        # ([position], [+,-], [data])
        red_dict[seq.id] = np.zeros((seq_size, len(STRANDDECODER), len(REDENCODER)-1), dtype=float)
        omvs[seq.id] = np.zeros((seq_size, len(STRANDDECODER)), dtype=object)

        for pos in range(seq_size):
            omvs[seq.id][pos, 0] = omv.LocShift(omv_size)
            omvs[seq.id][pos, 1] = omv.LocShift(omv_size)

    return red_dict, omvs

def buildModels(red_dict : dict, omvs : dict, nano2readid : dict, readid2fast5 : dict, nanopolish_result_csv : str, calculate_data_density : bool, max_lines : int = None):

    for loop in ['building models', 'checking data']:
        max_mem = 0

        with open(nanopolish_result_csv, 'r') as nano_result:
            # skip header
            nano_result.readline()
            opened_readid = ''
            opened_fast5 = ''
            last_position = None
            last_contig = None
            event = {}
            strand = 0

            for lidx, line in enumerate(nano_result):

                if (lidx + 1) % 100000 == 0:
                    max_mem = max(max_mem, sizeof_fmt(PROCESS.memory_info().rss))
                    print(f'Line {lidx + 1}{f"/{max_lines}" if max_lines is not None else ""}, {loop}, max memory usage: {max_mem}\t\t', end = '\r')
                    if max_lines is not None and lidx >= max_lines:
                        break

                read_line(line, event) # read data from nanopolish and store into event dictionary
                if event['model_kmer'] == 'NNNNN': # maybe haplotypes end up here as NNNNN? -> actually mutations in the reads, nanopolish has no clue what to do?
                    continue
                
                # prepare signal data for new read
                readid = nano2readid[event['read_index']]
                # found new read, store last information
                if readid != opened_readid:

                    # new read -> count last position of previous read on previous contig (if not first read)
                    if last_position is not None:
                        red_dict[last_contig][last_position, strand, REDENCODER['n_reads']] += 1
                        last_position = None

                    if event['ref_kmer'] == event['model_kmer']:
                        strand = 0
                    else:
                        strand = 1

                    # only open new fast5 file, if necessary
                    if readid2fast5[readid] != opened_fast5:
                        if opened_fast5: # true if opened_fast5 != ''
                            fast5.close()

                        fast5 = h5py.File(readid2fast5[readid], 'r')
                        opened_fast5 = readid2fast5[readid]
                        
                    fast5_read = fast5['read_' + readid]
                    # get normalised read signal
                    norm_signal = getZNormSignal(fast5_read)
                    opened_readid = readid

                # extract segment signal from read
                segment = norm_signal[event['start_idx']:event['end_idx']]
                red = red_dict[event['contig']][event['position'], strand]
                omv = omvs[event['contig']][event['position'], strand]

                if loop == 'building models':
                
                    # online bayesian updating
                    omv.append(segment)
                    red[REDENCODER['n_datapoints']] += len(segment)
                    red[REDENCODER['n_segments']] += 1

                    # see a new position within a previously opened read: add 1 to n_reads counter
                    if event['position'] != last_position:
                        if last_position is not None:
                            red_dict[last_contig][last_position, strand, REDENCODER['n_reads']] += 1
                        last_position = event['position']
                        last_contig = event['contig']

                elif loop =='checking data':

                    red[REDENCODER['mean']], red[REDENCODER['std']] = omv.meanStdev()

                    r = 3 * red[REDENCODER['std']] # ~99% density of normal distribution
                    l = red[REDENCODER['mean']] - r # lower bound
                    u = red[REDENCODER['mean']] + r # upper bound

                    contained = ((l <= segment) & (segment <= u)).sum()
                    if contained == len(segment):
                        red[REDENCODER['contained_segments']] += 1
                    red[REDENCODER['contained_datapoints']] += contained
                    
                    if calculate_data_density:
                        red[REDENCODER['data_density']] += np.mean(stats.norm.pdf(segment, loc = red[REDENCODER['mean']], scale = red[REDENCODER['std']]))
        
            LOGGER.printLog(f'Line {lidx + 1}, {loop}, max memory usage: {max_mem}\t\t', newline_after=True)

    if calculate_data_density: 
        # normalize log density
        for contig in red_dict:
            for positions in red_dict[contig]:
                for strand in red_dict[contig][positions]:
                    if red_dict[contig][positions, strand, REDENCODER['n_segments']] != 0:
                        # compare data density with the expected model density -> how good describes my model the data
                        red_dict[contig][positions, strand, REDENCODER['data_density']] = red_dict[contig][positions, strand, REDENCODER['data_density']] / red_dict[contig][positions, strand, REDENCODER['n_segments']]

def writeOutput(red_file : str, red_dict : dict):

    nans = 0

    with open(red_file, 'w') as w:

        w.write('reference\tposition\tstrand\tsignal_mean\tsignal_std\tdata_density\texpected_model_density\tn_datapoints\tcontained_datapoints\tn_segments\tcontained_segments\tn_reads\n')

        for seq_id in red_dict:
            for pos in range(len(red_dict[seq_id])):
                for strand in [0, 1]:

                    data = red_dict[seq_id][pos, strand]

                    if data[REDENCODER['std']]: # True if std != 0
                        expected_model_density = 1 / (2 * np.sqrt(np.pi) * data[REDENCODER['std']])
                    else:
                        nans += 1
                        expected_model_density = np.nan

                    w.write(f'{seq_id}\t{pos}\t{STRANDDECODER[strand]}\t{data[REDENCODER["mean"]]:.8f}\t{data[REDENCODER["std"]]:.8f}\t{data[REDENCODER["data_density"]]:.8f}\t{expected_model_density:.8f}\t{data[REDENCODER["n_datapoints"]]:.0f}\t{data[REDENCODER["contained_datapoints"]]:.0f}\t{data[REDENCODER["n_segments"]]:.0f}\t{data[REDENCODER["contained_segments"]]:.0f}\t{data[REDENCODER["n_reads"]]:.0f}\n')

def main() -> None:
    
    args = parse()
    
    working_dir = args.working_dir
    path_to_fast5 = args.path_to_fast5
    path_to_reference = args.path_to_reference
    guppy_bin = args.guppy_bin
    guppy_model = args.guppy_model
    guppy_device = args.guppy_device
    sample_label = args.sample_label
    fast5_out = args.fast5_out
    threads = args.threads
    force_rebuild = args.force_rebuild
    path_to_basecalls = args.path_to_basecalls
    mx = args.minimap2x
    mk = args.minimap2k
    max_lines = args.max_lines
    calculate_data_density = args.calculate_data_density
    if args.error_prefix is not None:
        ERROR_PREFIX = args.error_prefix

    assert (guppy_bin is not None and guppy_model is not None) or path_to_basecalls is not None, 'Need at least the guppy binary path and model or path to basecalls'

    global TIMEIT
    TIMEIT = args.timeit
        
    log_file = os.path.join(working_dir, 'log', f'{sample_label}_magnipore_{datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}.log')
    if not os.path.exists(os.path.join(working_dir, 'log')):
        os.makedirs(os.path.join(working_dir, 'log'))
    global LOGGER
    LOGGER = Logger(open(log_file, 'w'))

    if not os.path.exists(working_dir):
        os.mkdir(working_dir)
            
    LOGGER.printLog(f'Starting magnipore pipeline. Writing log to {log_file}')

    sequencing_summary = None

    if path_to_basecalls is not None:
        sequencing_summary = os.path.join(path_to_basecalls, 'sequencing_summary.txt')
        path_to_basecalls = os.path.join(path_to_basecalls, sample_label)

        if not os.path.exists(path_to_basecalls + '.fastq'):
            if not os.path.exists(path_to_basecalls + '.fq'):
                LOGGER.error(f'{path_to_basecalls}.fastq or .fq NOT FOUND!', error_type=ERROR_PREFIX+'26')
            else:
                path_to_basecalls+='.fq'
        else:
            path_to_basecalls+='.fastq'

        if not os.path.exists(sequencing_summary):
            LOGGER.error(f'{sequencing_summary} NOT FOUND!', error_type=ERROR_PREFIX+'27')

    else:

        path_to_basecalls, sequencing_summary, force_rebuild = guppy_basecalling(guppy_bin, guppy_model, guppy_device, path_to_fast5, working_dir, sample_label, fast5_out, force_rebuild, sequencing_summary)

    # new alignment/mapping for nanopolish with the corrected reference
    alignment_bam, force_rebuild = minimap(path_to_reference, path_to_basecalls, working_dir, sample_label, threads, force_rebuild, mx, mk)
    
    nanopolish_summary_csv, nanopolish_result_csv, force_rebuild = nanopolish(path_to_fast5, sequencing_summary, path_to_basecalls, path_to_reference, alignment_bam, working_dir, sample_label, threads, force_rebuild)
    
    red_file, force_rebuild = aggregate_events(nanopolish_result_csv, nanopolish_summary_csv, path_to_fast5, path_to_reference, working_dir, sample_label, force_rebuild, sequencing_summary, calculate_data_density, max_lines)

    LOGGER.printLog(f'Aggregated reference event distributions for {sample_label} are stored in {red_file}')

if __name__ == '__main__':
    main()
