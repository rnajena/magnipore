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
import seaborn as sns
from matplotlib.pyplot import figure
from matplotlib import pyplot as plt
from matplotlib import rcParams
from time import perf_counter_ns
import os
import psutil

from magnipore.__init__ import __version__
from magnipore.Logger import Logger
from magnipore.Helper import ANSI, DATAENCODER, STRANDDECODER, complement
import magnipore.OnlineMeanVar as omv

PROCESS = psutil.Process(os.getpid())
LOGGER : Logger = None
TIMEIT = False

def sizeof_fmt(num, suffix="B"):
    for unit in ["", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi"]:
        if abs(num) < 1024.0:
            return f"{num:3.1f}{unit}{suffix}"
        num /= 1024.0
    return f"{num:.1f}Yi{suffix}"

def readFast5(fast5_path : str, sequencing_summary : str):
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
                    LOGGER.printLog(f'Indexing readid {ridx + 1}\r', newline_after=False)
            
                filename, read_id = line.strip().split('\t')[:2]
                readid2fast5[read_id] = os.path.join(fast5_path, filename)
    
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
        # description='Required tools can be installed using conda and : guppy, minimap2, nanopolish, h5py, samtools, scipy, mafft, numpy\nsee github https://github.com/JannesSP/magnipore'
        prog='Nanosherlock',
        ) 
    
    parser.add_argument('path_to_fast5', type = str, help='FAST5 file')
    parser.add_argument('path_to_reference', type = str, help='reference FASTA file')
    parser.add_argument('working_dir', type = str, help='Path to write all output files')
    parser.add_argument('sample_label', type = str, help='Name of the sample or pipeline run')

    parser.add_argument('--guppy_bin', type = str, default = None, help='Guppy binary')
    parser.add_argument('--guppy_model', type = str, default = None, help='Guppy model for basecalling')
    parser.add_argument('--guppy_device', type=str, default=None, help='Use the gpu to basecall with cuda:0')

    parser.add_argument('--path_to_basecalls', default = None, type = str, help = 'Path to existing basecalls and sequencing summary file. Basecalls must be in one single file with the name <sample_label>.fastq')
    parser.add_argument('--calculate_data_density', action = 'store_true', default = False, help = 'Will calculate data density after building the models. Will increase runtime!')
    parser.add_argument('-t', '--threads', type=int, help='Number of threads to use')
    parser.add_argument('-f5', '--fast5_out', action = 'store_true', help='Guppy generates FAST5 output (workspace folder)')
    parser.add_argument('-fr', '--force_rebuild', action = 'store_true', help='Run commands regardless if files are already present')
    parser.add_argument('-mx', '--minimap2x', default = 'splice', choices = ['map-ont', 'splice', 'ava-ont'], help = '-x parameter for minimap2')
    parser.add_argument('-mk', '--minimap2k', default = 14, type = int, help = '-k parameter for minimap2')
    parser.add_argument('--timeit', default = False, action = 'store_true', help = 'Measure and print time used by submodules')
    parser.add_argument('--max_lines', default=None, type=int, help='Only process first given number of lines from nanopolish eventalign')

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
    
    if TIMEIT:
        end = perf_counter_ns()

    if ret != 0:
        LOGGER.error(f'Error in guppy basecalling with error code {ret}')
        
    if TIMEIT:
        LOGGER.printLog(f'TIMED: Guppy basecalling took {pd.to_timedelta(end-start)}, {end - start} nanoseconds')

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
    
    if TIMEIT:
        end = perf_counter_ns()

    if ret != 0:
        LOGGER.error(f'Error in minimap2 run with error code {ret}')
    
    if TIMEIT:
        LOGGER.printLog(f'TIMED: minimap2 took {pd.to_timedelta(end-start)}, {end - start} nanoseconds')

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
        
        if TIMEIT:
            end = perf_counter_ns()

        if ret != 0:
            LOGGER.error(f'Error in samtools indexing for nanopolish with error code {ret}')
        
        if TIMEIT:
            LOGGER.printLog(f'TIMED: samtools indexing took {pd.to_timedelta(end-start)}, {end - start} nanoseconds')

    else:
        LOGGER.printLog(f"Using existing {os.path.exists(alignment_bam + '.bai')}")

    if not os.path.exists(path_to_basecalls + '.index') or force_rebuild:

        # LOGGER.printLog(f'Run nanopolish indexing')
        command = f'nanopolish index -d {path_to_fast5} -s {path_to_sequencing_summary} {path_to_basecalls}'
        LOGGER.printLog(f'nanopolish indexing command: {ANSI.GREEN}{command}{ANSI.END}')

        if TIMEIT:
            start = perf_counter_ns()

        ret = os.system(command)

        if TIMEIT:
            end = perf_counter_ns()
        
        if ret != 0:
            LOGGER.error(f'Error in nanopolish indexing with error code {ret}')

        if TIMEIT:
            LOGGER.printLog(f'TIMED: nanopolish indexing took {pd.to_timedelta(end-start)}, {end - start} nanoseconds')

    else:
        LOGGER.printLog(f"Using existing {os.path.exists(path_to_basecalls + '.index')}")

    
    # LOGGER.printLog(f'Run nanopolish eventalign')
    command = f'nanopolish eventalign --reads {path_to_basecalls} --bam {alignment_bam} --genome {path_to_reference} --summary={summary_csv} --scale-events --signal-index -t {threads} > {result_csv}'
    LOGGER.printLog(f'Nanopolish eventalign command: {ANSI.GREEN}{command}{ANSI.END}')
    
    if TIMEIT:
        start = perf_counter_ns()

    ret = os.system(command)
    
    if TIMEIT:
        end = perf_counter_ns()

    if ret != 0:
        LOGGER.error(f'Error in nanopolish eventalign with error code {ret}')
    
    if TIMEIT:
        LOGGER.printLog(f'TIMED: nanopolish eventalign took {pd.to_timedelta(end-start)}, {end - start} nanoseconds')

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
    sequences = createSeqDict(path_to_reference)
    
    LOGGER.printLog(f'Found sequences: {list(sequences.keys())}. Start getting event signal distributions per reference position and write red file for {sample_label}')

    if TIMEIT:
        start = perf_counter_ns()

    buildModels(sequences, nano2readid, readid2fast5, nanopolish_result_csv, calculate_data_density, max_lines)

    if TIMEIT:
        end = perf_counter_ns()

    if TIMEIT:
        LOGGER.printLog(f'TIMED: Building distribution models took {pd.to_timedelta(end-start)}, {end - start} nanoseconds')

    LOGGER.printLog('Writing output files')
    writeOutput(red_file, sequences, working_dir, sample_label, calculate_data_density)
        
    return red_file, force_rebuild

def createSeqDict(path_to_reference : str) -> dict:

    sequences = {}
    fasta = SeqIO.parse(open(path_to_reference), 'fasta')

    # get all sequences from the reference file
    # sequences are stored as {reference: {pos: {strand: {base: (A|C|G|T), signal: []}}}}
    for seq in fasta:
        
        # ([position], [+,-], [data])
        sequences[seq.id] = np.zeros((len(seq.seq), 2, 11), dtype=object)
        
        for pos, base in enumerate(seq.seq):
            
            sequences[seq.id][pos, 0, DATAENCODER['base']] = base
            sequences[seq.id][pos, 0, DATAENCODER['omv']] = omv.LocShift(30)

            sequences[seq.id][pos, 1, DATAENCODER['base']] = complement(base)
            sequences[seq.id][pos, 1, DATAENCODER['omv']] = omv.LocShift(30)

            if 3 <= pos <= sequences[seq.id].shape[0] - 4:
                sequences[seq.id][pos, 0, DATAENCODER['motif']] = seq.seq[pos-3:pos+4]
                sequences[seq.id][pos, 1, DATAENCODER['motif']] = complement(seq.seq[pos-3:pos+4])[::-1]
            elif 2 <= pos <= sequences[seq.id].shape[0] - 3:
                sequences[seq.id][pos, 0, DATAENCODER['motif']] = seq.seq[pos-2:pos+3]
                sequences[seq.id][pos, 1, DATAENCODER['motif']] = complement(seq.seq[pos-2:pos+3])[::-1]


    return sequences

def buildModels(sequences : dict, nano2readid : dict, readid2fast5 : dict, nanopolish_result_csv : str, calculate_data_density : bool, max_lines : int):

    for loop in ['building models', 'checking data']:
        # LOGGER.printLog(f'Starting {loop}')

        with open(nanopolish_result_csv, 'r') as nano_result:
            # skip header
            nano_result.readline()
            opened_readid = ''
            opened_fast5 = ''
            last_position = None
            event = {}
            strand = 0

            for lidx, line in enumerate(nano_result):

                if (lidx + 1) % 100000 == 0:
                    print(f'Line {lidx + 1}{f"/{max_lines}" if max_lines is not None else ""}, {loop}, memory usage: {sizeof_fmt(PROCESS.memory_info().rss)}\t\t', end = '\r')
                    if max_lines is not None and lidx >= max_lines:
                        break

                read_line(line, event) # read data from nanopolish and store into event dictionary
                if event['model_kmer'] == 'NNNNN': # maybe haplotypes end up here as NNNNN? -> actually mutations in the reads, nanopolish has no clue what to do?
                        continue
                
                # prepare signal data for new read
                readid = nano2readid[event['read_index']]
                if readid != opened_readid:

                    # new read -> count last position of last read if not first read
                    if last_position is not None:
                        sequences[event['contig']][last_position, strand, DATAENCODER['n_reads']] += 1
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
                red = sequences[event['contig']][event['position'], strand] 

                if loop == 'building models':
                
                    # online bayesian updating
                    red[DATAENCODER['omv']].append(segment)
                    red[DATAENCODER['n_datapoints']] += len(segment)
                    red[DATAENCODER['n_segments']] += 1

                    # see a new position within a previously opened read: add 1 to n_reads counter
                    if event['position'] != last_position:
                        if last_position is not None:
                            sequences[event['contig']][last_position, strand, DATAENCODER['n_reads']] += 1
                        last_position = event['position']

                elif loop =='checking data':

                    red[DATAENCODER['mean']], red[DATAENCODER['std']] = red[DATAENCODER['omv']].meanStdev()

                    r = 3 * red[DATAENCODER['std']] # ~99% density of normal distribution
                    l = red[DATAENCODER['mean']] - r # lower bound
                    u = red[DATAENCODER['mean']] + r # upper bound

                    contained = ((l <= segment) & (segment <= u)).sum()
                    if contained == len(segment):
                        red[DATAENCODER['contained_segments']] += 1
                    red[DATAENCODER['contained_datapoints']] += contained
                    
                    if calculate_data_density:
                        red[DATAENCODER['data_density']] += np.mean(stats.norm.pdf(
                            segment,
                            loc = red[DATAENCODER['mean']],
                            scale = red[DATAENCODER['std']]))
                    
        LOGGER.printLog(f'Line {lidx + 1}, {loop}, memory usage: {sizeof_fmt(PROCESS.memory_info().rss)}\t\t\r', newline_after=False)
        print()

    if calculate_data_density: 
        # normalize log density
        for contig in sequences:
            for positions in sequences[contig]:
                for strand in sequences[contig][positions]:
                    if sequences[contig][positions, strand, DATAENCODER['n_segments']] != 0:
                        # compare data density with the expected model density -> how good describes my model the data
                        sequences[contig][positions, strand, DATAENCODER['data_density']] = sequences[contig][positions, strand, DATAENCODER['data_density']] / sequences[contig][positions, strand, DATAENCODER['n_segments']]

def writeOutput(red_file : str, sequences : dict, working_dir : str, sample_label : str, calculate_data_density : bool):

    nans = 0
    plotting_data = pd.DataFrame(columns=['contig', 'position', 'strand', 'density difference', 'n_datapoints', 'n_segments', 'n_reads'])
    plotting_data = plotting_data.astype(
        {
            'contig' : 'str',
            'position' : 'int',
            'strand' : 'str',
            'density difference': 'float',
            'n_datapoints' : 'int',
            'n_segments' : 'int',
            'n_reads' : 'int'
        })
    
    with open(red_file, 'w') as w:

        w.write('reference\tposition\tstrand\tbase\tsignal_mean\tsignal_std\tmotif\tdata_density\texpected_model_density\tn_datapoints\tcontained_datapoints\tn_segments\tcontained_segments\tn_reads\n')

        for sequence in sequences:
            for position in range(len(sequences[sequence])):
                for strand in range(len(sequences[sequence][position])):

                    data = sequences[sequence][position, strand]

                    if data[DATAENCODER['std']]: # True if std != 0
                        expected_model_density = 1 / (2 * np.sqrt(np.pi) * data[DATAENCODER['std']])
                        new_entry = pd.DataFrame({
                            'contig': [sequence],
                            'position': [position],
                            'strand': [STRANDDECODER[strand]],
                            'density difference': [expected_model_density - data[DATAENCODER["data_density"]] if calculate_data_density else expected_model_density],
                            'n_datapoints' : [data[DATAENCODER["n_datapoints"]]],
                            'n_segments' : [data[DATAENCODER["n_segments"]]],
                            'n_reads' : [data[DATAENCODER["n_reads"]]]

                        })
                        plotting_data = pd.concat([plotting_data, new_entry], ignore_index=True)
                    else:
                        nans += 1
                        expected_model_density = np.nan

                    w.write(f'{sequence}\t{position}\t{STRANDDECODER[strand]}\t{data[DATAENCODER["base"]]}\t{data[DATAENCODER["mean"]]}\t{data[DATAENCODER["std"]]}\t{data[DATAENCODER["motif"]]}\t{data[DATAENCODER["data_density"]]}\t{expected_model_density}\t{data[DATAENCODER["n_datapoints"]]}\t{data[DATAENCODER["contained_datapoints"]]}\t{data[DATAENCODER["n_segments"]]}\t{data[DATAENCODER["contained_segments"]]}\t{data[DATAENCODER["n_reads"]]}\n')             

    # LOGGER.printLog(f'Positions without information: {nans}')

    plotStatistics(plotting_data, working_dir, sample_label, calculate_data_density)

def plotStatistics(dataFrame : pd.DataFrame, working_dir : str, sample_label : str, calculate_data_density : bool) -> None:
    LOGGER.printLog('Plotting data ...')
    working_dir = os.path.join(working_dir, 'magnipore', sample_label, 'plots')
    rcParams['agg.path.chunksize'] = 10000

    if not os.path.exists(working_dir):
        os.mkdir(working_dir)

    dataFrame.replace({'Strand':STRANDDECODER})
    dataFrame['Contig, Strand'] =  pd.Series(dataFrame.reindex(['contig', 'strand'], axis='columns').astype('str').values.tolist()).str.join(', ')

    try:
        figure(figsize = (12,8), dpi=2000)
        g = sns.lineplot(data = dataFrame, x = 'position', y = 'n_reads', hue = 'Contig, Strand')
        g.set_xlim((0, max(dataFrame['position'])))
        g.set_xticks(range(0, max(dataFrame['position']), max(dataFrame['position'])//10))
        plt.title(f'Read coverage of segmented signals for sample {sample_label}')
        plt.setp(g.get_legend().get_texts(), fontsize='6') # for legend text
        plt.setp(g.get_legend().get_title(), fontsize='6') # for legend title
        g.set_yscale("log")
        plt.grid(True, 'both', 'both', alpha=0.6, color='grey')
        plt.tight_layout()
        plt.savefig(os.path.join(working_dir, f'{sample_label}_readCoverage.png'))
        plt.savefig(os.path.join(working_dir, f'{sample_label}_readCoverage.pdf'))
        plt.close()
    except:
        plt.close()
        LOGGER.warning('Plotting read coverage failed')

    try:
        figure(figsize = (12,8), dpi=2000)
        g = sns.lineplot(data = dataFrame, x = 'position', y = 'n_segments', hue = 'Contig, Strand')
        g.set_xlim((0, max(dataFrame['position'])))
        g.set_xticks(range(0, max(dataFrame['position']), max(dataFrame['position'])//10))
        plt.title(f'Segment coverage of segmented signals for sample {sample_label}')
        plt.setp(g.get_legend().get_texts(), fontsize='6') # for legend text
        plt.setp(g.get_legend().get_title(), fontsize='6') # for legend title
        g.set_yscale("log")
        plt.grid(True, 'both', 'both', alpha=0.6, color='grey')
        plt.tight_layout()
        plt.savefig(os.path.join(working_dir, f'{sample_label}_segmentCoverage.png'))
        plt.savefig(os.path.join(working_dir, f'{sample_label}_segmentCoverage.pdf'))
        plt.close()
    except:
        plt.close()
        LOGGER.warning('Plotting segment coverage failed')

    try:
        figure(figsize = (12,8), dpi=2000)
        g = sns.lineplot(data = dataFrame, x = 'position', y = 'n_datapoints', hue = 'Contig, Strand')
        g.set_xlim((0, max(dataFrame['position'])))
        g.set_xticks(range(0, max(dataFrame['position']), max(dataFrame['position'])//10))
        plt.title(f'Signal coverage of segmented signals for sample {sample_label}')
        plt.setp(g.get_legend().get_texts(), fontsize='6') # for legend text
        plt.setp(g.get_legend().get_title(), fontsize='6') # for legend title
        g.set_yscale("log")
        plt.grid(True, 'both', 'both', alpha=0.6, color='grey')
        plt.tight_layout()
        plt.savefig(os.path.join(working_dir, f'{sample_label}_signalCoverage.png'))
        plt.savefig(os.path.join(working_dir, f'{sample_label}_signalCoverage.pdf'))
        plt.close()
    except:
        plt.close()
        LOGGER.warning('Plotting signal coverage failed')

    try:
        figure(figsize = (12,8), dpi=2000)
        g = sns.histplot(data = dataFrame, x = 'density difference', kde = True, hue = 'Contig, Strand', stat = 'density')
        plt.title('Distribution of model density vs data density difference of all position.' if calculate_data_density else 'Model density per position.')
        plt.setp(g.get_legend().get_texts(), fontsize='6') # for legend text
        plt.setp(g.get_legend().get_title(), fontsize='6') # for legend title
        plt.tight_layout()
        plt.savefig(os.path.join(working_dir, f'{sample_label}_densityDiff.png' if calculate_data_density else f'{sample_label}_modeldensity.png'))
        plt.savefig(os.path.join(working_dir, f'{sample_label}_densityDiff.pdf' if calculate_data_density else f'{sample_label}_modeldensity.pdf'))
        plt.close()
    except:
        plt.close()
        LOGGER.warning('Plotting density failed')

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
    # medaka_model = args.medaka_model
    calculate_data_density = args.calculate_data_density

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
        # LOGGER.printLog(f'Creating working directory {working_dir}')
            
    LOGGER.printLog(f'Starting magnipore pipeline. Writing log to {log_file}')

    sequencing_summary = None

    if path_to_basecalls is not None:
        sequencing_summary = os.path.join(path_to_basecalls, 'sequencing_summary.txt')
        path_to_basecalls = os.path.join(path_to_basecalls, sample_label)

        if not os.path.exists(path_to_basecalls + '.fastq'):
            if not os.path.exists(path_to_basecalls + '.fq'):
                LOGGER.error(f'{path_to_basecalls}.fastq or .fq NOT FOUND!')
            else:
                path_to_basecalls+='.fq'
        else:
            path_to_basecalls+='.fastq'

        if not os.path.exists(sequencing_summary):
            LOGGER.error(f'{sequencing_summary} NOT FOUND!')

    else:

        path_to_basecalls, sequencing_summary, force_rebuild = guppy_basecalling(guppy_bin, guppy_model, guppy_device, path_to_fast5, working_dir, sample_label, fast5_out, force_rebuild, sequencing_summary)

    # new alignment/mapping for nanopolish with the corrected reference
    alignment_bam, force_rebuild = minimap(path_to_reference, path_to_basecalls, working_dir, sample_label, threads, force_rebuild, mx, mk)
    
    nanopolish_summary_csv, nanopolish_result_csv, force_rebuild = nanopolish(path_to_fast5, sequencing_summary, path_to_basecalls, path_to_reference, alignment_bam, working_dir, sample_label, threads, force_rebuild)
    
    red_file, force_rebuild = aggregate_events(nanopolish_result_csv, nanopolish_summary_csv, path_to_fast5, path_to_reference, working_dir, sample_label, force_rebuild, sequencing_summary, calculate_data_density, max_lines)

    LOGGER.printLog(f'Aggregated reference event distributions for {sample_label} are stored in {red_file}')

    # return red_file

if __name__ == '__main__':
    main()
