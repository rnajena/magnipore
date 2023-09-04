#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import datetime
import os
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from pathlib import Path
from time import perf_counter_ns

import numpy as np
import pandas as pd
import psutil
from Bio import SeqIO
from read5 import read
from scipy import stats

import magnipore.OnlineMeanVar as omv
from magnipore.__init__ import __version__
from magnipore.Helper import ANSI, REDENCODER, STRANDDECODER, sizeof_fmt
from magnipore.Logger import Logger

PROCESS = psutil.Process(os.getpid())
LOGGER : Logger = None
TIMEIT = False
ERROR_PREFIX : str = '0'

def initLogger(file = None) -> None:
    global LOGGER
    LOGGER = Logger(file)

def mapFast5Files(raw_data_path : str, seq_sum : str = None) -> dict:
    '''
    Returns a dictionary mapping the readid to the corresponding fast5 file containing the read.
    
    Parameters
    ----------
    raw_data_path : str
        Path or directory containing the fast5 files
    seq_sum : str = None:
        Path to the sequencing summary file. Reduced time to index files.
    '''
    LOGGER.printLog(f'Loading readids from fast5 files in {raw_data_path} ... ')
    readid2file = {}

    if seq_sum is None:
        # recursively loop over the fast5 files in the given path
        fast5list = Path(raw_data_path).rglob('*.fast5')
        for fidx, fast5 in enumerate(fast5list):
            LOGGER.printLog(f'Indexing fast5 file {fidx + 1}\r', newline_after=False)
            # because path is object not string
            f5 = read(str(fast5))
            for readid in f5.getReads():
                readid2file[readid] = str(fast5)
            f5.close()

    else:
        with open(seq_sum, 'r') as seqsum:
            seqsum.readline()
            for ridx, line in enumerate(seqsum):
                if (ridx + 1) % 10000 == 0:
                    LOGGER.printLog(f'Indexing read {ridx + 1}\r', newline_after=False)
                filename, read_id = line.strip().split('\t')[:2]
                readid2file[read_id] = os.path.join(raw_data_path, filename)
    
    LOGGER.printLog(f'Indexed {len(readid2file) + 1} reads')
    return readid2file

def readSegSum(segSum : str):
    '''
    Read summary file with read_indices and read_ids.\n

    Parameters
    ----------
    segSum : str
        path to segmentation summary file

    Returns
    -------
    read_index2ID : dict
        mapping between segmentation index and read ID
    '''
    read_index2ID = {}
    with open(segSum, 'r') as summary:
        next(summary) # skip header
        for lidx, line in enumerate(summary):
            if (lidx + 1) % 10000 == 0:
                print(f'Reading line {lidx + 1}', end='\r')
            r_index, r_ID = line.strip().split()[:2]
            read_index2ID[int(r_index)] = r_ID
    LOGGER.printLog(f'Mapped {len(read_index2ID)} readids')
    return read_index2ID

def parse() -> Namespace:

    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description='Required tools: see github https://github.com/JannesSP/magnipore',
        prog='Nanosherlock',
        ) 
    
    parser.add_argument('raw_data', type = str, help='Parent directory of FAST5 files, can also be a direct path to a single SLOW5 or BLOW5 file, that contains all reads, if FASTQs are provided')
    parser.add_argument('reference', type = str, help='reference FASTA file, POSITIVE (+) or FORWARD strand')
    parser.add_argument('working_dir', type = str, help='Path to write all output files')
    parser.add_argument('sample_label', type = str, help='Name of the sample or pipeline run')
    
    parser.add_argument('--guppy_bin', type = str, default = None, help='Guppy binary')
    parser.add_argument('--guppy_model', type = str, default = None, help='Guppy model for basecalling')
    parser.add_argument('--guppy_device', type=str, default=None, help='Use the gpu to basecall with cuda:0')

    parser.add_argument('--basecalls', metavar='FASTQ', default = None, type = str, help = 'Path to existing basecalls. Basecalls must be in one single file.')
    parser.add_argument('--sequencing_summary', metavar='TXT', type = str, default = None, help = 'Use, when sequencing summary is not next to your FASTQ file. Path to existing sequencing summary file of sample.')
    parser.add_argument('--calculate_data_density', action = 'store_true', default = False, help = 'Will calculate data density after building the models. Will increase runtime!')
    parser.add_argument('-t', '--threads', type=int, help='Number of threads to use')
    parser.add_argument('-fr', '--force_rebuild', action = 'store_true', help='Run commands regardless if files are already present')
    parser.add_argument('-mx', '--minimap2x', default = 'map-ont', choices = ['map-ont', 'splice', 'ava-ont'], help = '-x parameter for minimap2')
    parser.add_argument('-mk', '--minimap2k', default = 14, type = int, help = '-k parameter for minimap2')
    parser.add_argument('--timeit', default = False, action = 'store_true', help = 'Measure and print time used by submodules')
    parser.add_argument('--max_lines', default=None, type=int, help='Only process first given number of lines from segmentation eventalign')
    parser.add_argument('-rna', default=False, action='store_true', help='Use when data is rna')
    parser.add_argument('-r10', default=False, action='store_true', help='Use when data is from R10.4.1 flowcell')
    parser.add_argument('-km', '--kmer_model', default = None, type=str, help='custom kmer model file for f5c eventalign')
    parser.add_argument('-e', '--error_prefix', type=str, default=None)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s' + f' {__version__}')

    return parser.parse_args()

def guppy_basecalling(guppy_bin : str, guppy_model : str,  guppy_device : str, raw_data : str, working_dir : str, sample_label : str, force_rebuild : bool, seq_sum : str) -> str:

    basecalls_path = os.path.join(working_dir, 'basecalls', sample_label)
    basecalls = os.path.join(basecalls_path, sample_label + ".fastq")

    if seq_sum is None:
        seq_sum = os.path.join(basecalls_path, 'sequencing_summary.txt')

    if os.path.exists(seq_sum) and os.path.exists(basecalls) and not force_rebuild:
        LOGGER.printLog(f'basecalls and sequencing summary already exist:\n-\t{basecalls}\n-\t{seq_sum}')
        return basecalls, seq_sum, force_rebuild
    else:
        force_rebuild = True

    if not os.path.exists(basecalls_path):
        os.makedirs(basecalls_path)
    
    command = f'{guppy_bin} -i {raw_data} -s {basecalls_path} --disable_pings --disable_qscore_filtering --calib_detect -c {guppy_model}'
    
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
    # compress fastq.gz with level 3
    os.system(f'gzip -f -3 {basecalls}')
        
    return basecalls+'.gz', seq_sum, force_rebuild

def mapping(reference : str, basecalls : str, working_dir : str, sample_label : str, threads : int, force_rebuild : bool, mx : str, mk : int) -> str:
    
    mapping_path = os.path.join(working_dir, 'mapping', sample_label)
    bam_path = os.path.join(mapping_path, sample_label) + '.bam'

    if os.path.exists(bam_path) and not force_rebuild:
        LOGGER.printLog(f'Bam alignment/mapping file already exists\n-\t{bam_path}')
        return bam_path, force_rebuild
    else:
        force_rebuild = True
    
    if not os.path.exists(mapping_path):
        os.makedirs(mapping_path)
        
    command = f'minimap2 -a -x {mx} -k{mk} -t {threads} {reference} {basecalls} | samtools view -hbF4 | samtools sort > {bam_path}'
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

def signalSegmentation(raw_data : str, file_format : str, basecalls : str, reference : str, alignment_bam : str, working_dir : str, sample_label : str, threads : int, force_rebuild : bool, rna : bool, r10 : bool, kmer_model : str) -> tuple:
    
    segmentation_path = os.path.join(working_dir, 'segmentation', sample_label)
    summary_csv = os.path.join(segmentation_path, "eventalign_summary.csv")
    result_csv = os.path.join(segmentation_path, "eventalign_result.csv")

    if os.path.exists(summary_csv) and os.path.exists(result_csv) and not force_rebuild:
        LOGGER.printLog(f'segmentation files already exist:\n-\t{summary_csv}\n-\t{result_csv}')
        return summary_csv, result_csv, force_rebuild
    else:
        force_rebuild = True

    if not os.path.exists(segmentation_path):
        os.makedirs(segmentation_path)

    # samtools indexing
    if not os.path.exists(alignment_bam + '.bai') or force_rebuild:
        command = f'samtools index {alignment_bam}'
        LOGGER.printLog(f'samtools index command: {ANSI.GREEN}{command}{ANSI.END}')
        if TIMEIT:
            start = perf_counter_ns()
        ret = os.system(command)
        if ret != 0:
            LOGGER.error(f'Error in samtools indexing with error code {ret}', error_type=ERROR_PREFIX+'23')
        if TIMEIT:
            end = perf_counter_ns()
            LOGGER.printLog(f'{ANSI.YELLOW}TIMED: samtools indexing took {pd.to_timedelta(end-start)}, {end - start} nanoseconds{ANSI.END}')
    else:
        LOGGER.printLog(f"Using existing {os.path.exists(alignment_bam + '.bai')}")

    # segmentation indexing
    if not os.path.exists(basecalls + '.index') or force_rebuild:
        command = f'f5c index {"--slow5" if file_format == ".slow5" else "-d"} {raw_data} {basecalls} -t {threads} --iop {max(1, threads//2)} '
        LOGGER.printLog(f'segmentation indexing command: {ANSI.GREEN}{command}{ANSI.END}')
        if TIMEIT:
            start = perf_counter_ns()
        ret = os.system(command)
        if ret != 0:
            LOGGER.error(f'Error in segmentation indexing with error code {ret}', error_type=ERROR_PREFIX+'24')
        if TIMEIT:
            end = perf_counter_ns()
            LOGGER.printLog(f'{ANSI.YELLOW}TIMED: segmentation indexing took {pd.to_timedelta(end-start)}, {end - start} nanoseconds{ANSI.END}')
    else:
        LOGGER.printLog(f"Using existing {os.path.exists(basecalls + '.index')}")

    # segmentation
    log_file = os.path.join(segmentation_path, "log.txt")
    command = f'f5c eventalign -r {basecalls} -b {alignment_bam} -g {reference} --summary={summary_csv} --scale-events --signal-index --secondary=no{f" --slow5 {raw_data} " if file_format == ".slow5" else " "}{"--rna " if rna else ""}--collapse-events --iop {max(1, threads//2)} -t {threads} -o {result_csv}'
    if kmer_model is not None:
        command += f' --kmer-model {kmer_model}'
    command += ' --pore r10' if r10 else ' --pore r9'
    command += f' 2> {log_file}'
    LOGGER.printLog(f'Segmentation command: {ANSI.GREEN}{command}{ANSI.END}')
    if TIMEIT:
        start = perf_counter_ns()
    ret = os.system(command)
    if ret != 0:
        LOGGER.error(f'Error in segmentation with error code {ret}. Look into {log_file} for a detailed error log.', error_type=ERROR_PREFIX+'25')
    if TIMEIT:
        end = perf_counter_ns()
        LOGGER.printLog(f'{ANSI.YELLOW}TIMED: segmentation took {pd.to_timedelta(end-start)}, {end - start} nanoseconds{ANSI.END}')

    return summary_csv, result_csv, force_rebuild

def read_line(line : str, r10 : bool) -> dict:
    '''
    Fills the event dict with new data from given line.
    Offset for r9 is +2, offset for r10 is +4.

    Parameters
    ----------
    line : str
        eventalign line
    event : dict
        event dictionary to store relevant data
    r10 : bool
        flag to switch between kmer-position offset
    '''
    line = line.strip().split()
    event = {}
    event['contig'] = line[0]
    # add offset to position to collect the signals for the middle base of the 5mer event
    # position describes the first base of the event kmer in the reference
    # offset = 2 for 5mer models in r9
    # offset = 4 for 9mer models in r10
    # position is 0-based
    event['position'] = int(line[1]) + 4 if r10 else int(line[1]) + 2
    event['ref_kmer'] = line[2]
    event['read_index'] = int(line[3])
    event['model_kmer'] = line[9]
    event['start_idx'] = int(line[13])
    event['end_idx'] = int(line[14])
    return event

def aggregate_events(seg_result : str, seg_sum: str, raw_data : str, file_format : str, reference : str, working_dir : str, sample_label : str, force_rebuild : bool, seq_sum : str, calculate_data_density : bool, r10 : bool, max_lines : int = None) -> str:
    '''
    Returns
    -------
    red_file : str
        path to read event distribution file
    '''
    # write distributions into reference event distribution file
    red_file = os.path.join(working_dir, 'magnipore', sample_label, f'{sample_label}.red')

    if not os.path.exists(os.path.join(working_dir, 'magnipore', sample_label)):
        os.makedirs(os.path.join(working_dir, 'magnipore', sample_label))
    elif os.path.exists(red_file) and not force_rebuild:
        return red_file

    readID2File = mapFast5Files(raw_data, seq_sum) if file_format == '.fast5' else raw_data
    idx2readid = readSegSum(seg_sum)
    ref_seqs = SeqIO.parse(open(reference), 'fasta')
    red_dict, omvs = createREDDict(ref_seqs)
    
    LOGGER.printLog(f'Found sequences: {list(red_dict.keys())}. Start getting event signal distributions per reference position and write red file for {sample_label}')

    if TIMEIT:
        start = perf_counter_ns()
    buildModels(red_dict, omvs, idx2readid, readID2File, seg_result, calculate_data_density, r10, max_lines)
    del omvs # actively release memory
    if TIMEIT:
        end = perf_counter_ns()
        LOGGER.printLog(f'{ANSI.YELLOW}TIMED: Building distribution models took {pd.to_timedelta(end-start)}, {end - start} nanoseconds{ANSI.END}')

    LOGGER.printLog('Writing output files')
    writeOutput(red_file, red_dict)
    return red_file

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

def getFile(read2FileMap, readid : str) -> str:
    '''
    Parameters
    ----------
    read2FileMap : dict or str
        dict in case of .fast5 files, str in case of .slow5 or blow5 files
    readid : str
        current readid

    Returns
    -------
    Path to file containing read with given readid
    '''
    if type(read2FileMap) is dict:
        # .fast5
        return read2FileMap[readid]
    # .slow5 | .blow5
    return read2FileMap

def getRead5Reader(read2FileMap, readid : str, r5):
    '''
    Parameters
    ----------
    read2FileMap : dict or str
        dict in case of .fast5 files, str in case of .slow5 or blow5 files
    readid : str
        current readid
    r5 : AbstractFileReader
        previous FileReader from read5: Fast5Reader, Slow5Reader or Pod5Reader

    Returns
    -------
    r5 : AbstractFileReader
        new FileReader from read5: Fast5Reader, Slow5Reader or Pod5Reader
    '''
    file = getFile(read2FileMap, readid)
    # open first file
    if r5 is None:
        return read(file)    
    # read is still in the same file, should be always true for slow5 and blow5 as a single file is provided
    if file == r5._filepath:
        return r5
    # different file to open -> close old one and open new file    
    if r5 is not None:
        r5.close()
    return read(file)

def buildModels(red_dict : dict, omvs : dict, nano2readid : dict, readID2File : dict, segmentation_result_csv : str, calculate_data_density : bool, r10 : bool, max_lines : int = None):

    for loop in ['building models', 'checking data']:
        max_mem = 0
        current_read = ''
        r5 = None
        start = perf_counter_ns()

        with open(segmentation_result_csv, 'r') as nano_result:
            # skip header
            nano_result.readline()

            for lidx, line in enumerate(nano_result):

                if (lidx + 1) % 100000 == 0:
                    end = perf_counter_ns()
                    max_mem = max(max_mem, PROCESS.memory_info().rss)
                    print(f'Line {lidx + 1}{f"/{max_lines}" if max_lines is not None else ""}, {loop}, max memory usage: {sizeof_fmt(max_mem)}, {pd.to_timedelta(end-start)}', end = '\r')
                    start = end
                    if max_lines is not None and (lidx + 1) >= max_lines:
                        break

                event = read_line(line, r10) # read data from segmentation and store into event dictionary
                if event['model_kmer'].count('N') == len(event['model_kmer']):
                    # maybe haplotypes end up here as NNNNN? -> actually mutations in the reads, segmentation has no clue what to do?
                    continue
                
                # prepare signal data for new read
                readid = nano2readid[event['read_index']]
                # found new read, store last information
                if readid != current_read:
                    strand = int(event['ref_kmer'] != event['model_kmer'])
                    r5 = getRead5Reader(readID2File, readid, r5)
                    norm_signal = r5.getZNormSignal(readid)
                    current_read = readid

                red_dict[event['contig']][event['position'], strand, REDENCODER['n_reads']] += 1

                # extract segment signal from read
                segment = norm_signal[event['start_idx']:event['end_idx']]
                red = red_dict[event['contig']][event['position'], strand]
                omv = omvs[event['contig']][event['position'], strand]

                if loop == 'building models':
                    # online bayesian updating
                    omv.append(segment)
                    red[REDENCODER['n_datapoints']] += len(segment)
                    red[REDENCODER['n_segments']] += 1

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
        
            LOGGER.printLog(f'Line {lidx + 1}, {loop}, max memory usage: {sizeof_fmt(max_mem)}\t\t', newline_after=True)

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
    with open(red_file, 'w') as red:
        red.write('reference\tposition\tstrand\tsignal_mean\tsignal_std\tdata_density\texpected_model_density\tn_datapoints\tcontained_datapoints\tn_segments\tcontained_segments\tn_reads\n')
        for seq_id in red_dict:
            for pos in range(len(red_dict[seq_id])):
                for strand in [0, 1]:
                    data = red_dict[seq_id][pos, strand]
                    if data[REDENCODER['std']]: # True if std != 0
                        expected_model_density = 1 / (2 * np.sqrt(np.pi) * data[REDENCODER['std']])
                    else:
                        nans += 1
                        expected_model_density = np.nan
                    red.write(f'{seq_id}\t{pos}\t{STRANDDECODER[strand]}\t{data[REDENCODER["mean"]]:.8f}\t{data[REDENCODER["std"]]:.8f}\t{data[REDENCODER["data_density"]]:.8f}\t{expected_model_density:.8f}\t{data[REDENCODER["n_datapoints"]]:.0f}\t{data[REDENCODER["contained_datapoints"]]:.0f}\t{data[REDENCODER["n_segments"]]:.0f}\t{data[REDENCODER["contained_segments"]]:.0f}\t{data[REDENCODER["n_reads"]]:.0f}\n')

def getFileFormat(raw_data : str) -> str:
    '''
    Checks file format for .fast5, .slow5, .blow5 or .pod5 extension.
    Can exit with error code x20 if no raw data found or the file format is unknown.

    Returns
    -------
    file_format : str
        one of ['.fast5', '.slow5', '.pod5'],
        slow5 == blow5 in this case
    '''
    if not os.path.exists(raw_data):
        LOGGER.error(f'{raw_data} not found', error_type=ERROR_PREFIX+'20')

    file_format = None
    if not (raw_data.endswith('.slow5') or raw_data.endswith('.blow5')):
        files = [f for f in os.listdir(raw_data) if os.path.isfile(os.path.join(raw_data, f))]
        contains_fast5 = False
        for file in files:
            if file.endswith('.fast5'):
                contains_fast5 = True
                file_format = '.fast5'
        if not contains_fast5:
            LOGGER.error(f'Unknown file format or no .fast5 files in {raw_data}', error_type=ERROR_PREFIX+'20')
    else:
        if raw_data.endswith('.slow5') or raw_data.endswith('.blow5'):
            file_format = '.slow5'
        elif raw_data.endswith('.pod5'):
            file_format = '.pod5'

    LOGGER.writeLog(f'Found raw data in {file_format} format.')
    return file_format

def setupLogger(working_dir : str, sample_label : str) -> None:
    '''
    Creates log file for the current magnipore run

    Returns
    -------
    LOGGER : Logger
        Logger object to write to log file, format console output and write error messages
    '''
    global LOGGER
    log_file = os.path.join(working_dir, 'log', f'{sample_label}_nanosherlock_{datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}.log')
    if not os.path.exists(os.path.join(working_dir, 'log')):
        os.makedirs(os.path.join(working_dir, 'log'))
    LOGGER = Logger(open(log_file, 'w'))
    LOGGER.printLog(f'Starting magnipore pipeline. Writing log to {log_file}')

def getBasecalls(file_format : str, basecalls : str, guppy_bin : str, guppy_model : str, guppy_device : str, raw_data : str, working_dir : str, sample_label : str, force_rebuild : bool, seq_sum : str = None):
    '''
    Checks input for basecalls, file formats and start guppy basecalling if needed.

    Returns
    -------
    basecalls : str
        Path to fastq file
    seq_sum : str
        Path to sequencing_summary.txt file. Can be None.
    force_rebuild : bool
        Flag if magnipore rebuilds already existing files
    '''
    if basecalls is not None:
        if not os.path.exists(basecalls):
            LOGGER.error(f'{basecalls} NOT FOUND!', error_type=ERROR_PREFIX+'26')
        if seq_sum is None and os.path.exists(os.path.join(os.path.dirname(basecalls), 'sequencing_summary.txt')):
            seq_sum = os.path.join(os.path.dirname(basecalls), 'sequencing_summary.txt')
        return basecalls, seq_sum, force_rebuild
    else:
        if file_format == '.slow5':
            LOGGER.error('Cannot basecall .slow5 or .blow5 with guppy - only .fast5 or .pod5!', error_type=ERROR_PREFIX+'19')
        return guppy_basecalling(guppy_bin, guppy_model, guppy_device, raw_data, working_dir, sample_label, force_rebuild, None)

def main() -> None:
    global TIMEIT, ERROR_PREFIX
    args = parse()
    working_dir = args.working_dir
    raw_data = args.raw_data
    reference = args.reference
    guppy_bin = args.guppy_bin
    guppy_model = args.guppy_model
    guppy_device = args.guppy_device
    sample_label = args.sample_label
    threads = args.threads
    force_rebuild = args.force_rebuild
    basecalls = args.basecalls
    seq_sum = args.sequencing_summary
    mx = args.minimap2x
    mk = args.minimap2k
    max_lines = args.max_lines
    calc_data_density = args.calculate_data_density
    rna = args.rna
    r10 = args.r10
    kmer_model = args.kmer_model
    TIMEIT = args.timeit
    if args.error_prefix is not None:
        ERROR_PREFIX = args.error_prefix

    assert (guppy_bin is not None and guppy_model is not None) or basecalls is not None, 'Need at least the guppy binary path and model or path to basecalls'

    if not os.path.exists(working_dir):
        os.mkdir(working_dir)
    setupLogger(working_dir, sample_label)

    file_format = getFileFormat(raw_data)
    basecalls, seq_sum, force_rebuild = getBasecalls(file_format, basecalls, guppy_bin, guppy_model, guppy_device, raw_data, working_dir, sample_label, force_rebuild, seq_sum)
    # new alignment/mapping for segmentation with the corrected reference
    alignment_bam, force_rebuild = mapping(reference, basecalls, working_dir, sample_label, threads, force_rebuild, mx, mk)
    seg_summary_csv, seg_result_csv, force_rebuild = signalSegmentation(raw_data, file_format, basecalls, reference, alignment_bam, working_dir, sample_label, threads, force_rebuild, rna, r10, kmer_model)
    red_file = aggregate_events(seg_result_csv, seg_summary_csv, raw_data, file_format, reference, working_dir, sample_label, force_rebuild, seq_sum, calc_data_density, r10, max_lines)
    LOGGER.printLog(f'Aggregated reference event distributions for {sample_label} are stored in {red_file}')

if __name__ == '__main__':
    main()
