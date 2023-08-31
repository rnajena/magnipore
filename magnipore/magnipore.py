#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import multiprocessing as mp
import os
import re
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from statistics import NormalDist
from time import perf_counter_ns

import numpy as np
import pandas as pd
from Bio import Seq, SeqIO
from scipy.stats import ks_2samp

from magnipore.__init__ import __version__
from magnipore.Helper import (ANSI, IUPAC, MAGNIPORE_COLUMNS, MUTDECODER,
                              REDENCODER, STRANDDECODER, STRANDENCODER,
                              complement, rev_complement)
from magnipore.Logger import Logger

LOGGER : Logger = None
TIMEIT : bool = False
# multiprocessing setup
mp.set_start_method('fork')
allQueue = mp.Manager().Queue()
magnQueue = mp.Manager().Queue()
alipQueue = mp.Manager().Queue()
num_muts = mp.Value('I', 0)
sign_pos = mp.Value('I', 0)
no_data = mp.Value('I', 0)
low_cov_count = mp.Value('I', 0)
num_pos = mp.Value('I', 0)

def initLogger(file = None) -> None:
    global LOGGER
    LOGGER = Logger(file)

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description='Required tools: see github https://github.com/JannesSP/magnipore',
        prog='Magnipore',
        )
    parser.add_argument("raw_data_first_sample", type = str, help='Parent directory of FAST5 files of first sample, can also be a single SLOW5 or BLOW5 file of first sample, that contains all reads, if FASTQs are provided')
    parser.add_argument("reference_first_sample", type = str, help='reference FASTA file of first sample, POSITIVE (+) or FORWARD strand, ATTENTION: can only contain a single sequence')
    parser.add_argument("label_first_sample", type = str, help='Name of the sample or pipeline run')
    parser.add_argument("raw_data_sec_sample", type = str, help='Parent directory of FAST5 files of second sample, can also be SLOW5 or BLOW5 file of second sample, that contains all reads, if FASTQs are provided')
    parser.add_argument("reference_sec_sample", type = str, help='reference FASTA file of second sample, POSITIVE (+) or FORWARD strand, ATTENTION: can only contain a single sequence')
    parser.add_argument("label_sec_sample", type = str, help='Name of the sample or pipeline run')
    parser.add_argument("working_dir", type = str, help='Path to write all output files')
    parser.add_argument("--guppy_bin", type = str, default = None, help='Guppy binary')
    parser.add_argument("--guppy_model", type = str, default = None, help='Guppy model used for basecalling')
    parser.add_argument('--guppy_device', type=str, default='cuda:0', help='Use the GPU to basecall "cuda:0" to use the GPU with ID 0')
    parser.add_argument('-b1', '--basecalls_first_sample', metavar='FASTQ', type = str, default = None, help = 'Path to existing basecalls of first sample. Basecalls must be in one single file.')
    parser.add_argument('-b2', '--basecalls_sec_sample', metavar='FASTQ', type = str, default = None, help = 'Path to existing basecalls of second sample. Basecalls must be in one single file.')
    parser.add_argument('-s1', '--sequencing_summary_first_sample', metavar='TXT', type = str, default = None, help = 'Use, when sequencing summary is not next to your FASTQ file. Path to existing sequencing summary file of second sample.')
    parser.add_argument('-s2', '--sequencing_summary_sec_sample', metavar='TXT', type = str, default = None, help = 'Use, when sequencing summary is not next to your FASTQ file. Path to existing sequencing summary file of first sample.')
    parser.add_argument('-d', '--calculate_data_density', action = 'store_true', default = False, help = 'Will calculate data density after building the models. Will increase runtime!')
    parser.add_argument('-t', "--threads", type=int, default=1, help='Number of threads to use')
    parser.add_argument('-fr', '--force_rebuild', action = 'store_true', help='Run commands regardless if files are already present')
    parser.add_argument('-mx', '--minimap2x', default = 'map-ont', choices = ['map-ont', 'splice', 'ava-ont'], help = '-x parameter for minimap2')
    parser.add_argument('-mk', '--minimap2k', default = 14, help = '-k parameter for minimap2')
    parser.add_argument('--timeit', default = False, action = 'store_true', help = 'Measure and print time used by submodules')
    parser.add_argument('-rna', '--rna', default=False, action='store_true', help='Use when data is rna')
    parser.add_argument('-r10', '--r10', default=False, action='store_true', help='Use when data is from R10.4.1 flowcell')
    parser.add_argument('-km', '--kmer_model', default = None, type=str, help='custom kmer model file for f5c eventalign')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s' + f' {__version__}')
    return parser.parse_args()

def align(ref_first_sample : str, ref_sec_sample : str, first_sample_label : str, sec_sample_label : str, working_dir : str, threads : int, force_rebuild : bool):
    
    # write both references into one file
    alignment_path = os.path.join(working_dir, 'alignment')
    ref_both_samples = os.path.join(alignment_path, f'{first_sample_label}_{sec_sample_label}.fa')
    ref_alignment = os.path.join(alignment_path, f'{first_sample_label}_{sec_sample_label}.aln')
    
    if os.path.exists(ref_alignment) and not force_rebuild:
        return ref_alignment
    
    if not os.path.exists(alignment_path):
        os.makedirs(alignment_path)
    
    # ensure newline between reference files
    command = f'(cat {ref_first_sample}; echo ""; cat {ref_sec_sample}) > {ref_both_samples}'
    
    # LOGGER.printLog(f'Writing both references into one file {ref_both_samples}')
    LOGGER.printLog(f'cat command: {ANSI.GREEN}{command}{ANSI.END}')
    ret = os.system(command)
    
    if ret != 0:
        LOGGER.error(f'Error in concatenating both reference files with error code {ret}', error_type=11)
    
    # provided the same reference for both samples
    if ref_first_sample == ref_sec_sample:
        command = 'awk \'BEGIN{FS=" "}{if(!/>/){print tolower($0)}else{print $1}}\' ' + f'{ref_both_samples} > {ref_alignment}'
        LOGGER.printLog(f'same reference for both samples: {ANSI.GREEN}{command}{ANSI.END}')

    # difference references for each sample
    else:
        command = f'mafft --auto --quiet --thread {threads} {ref_both_samples} > {ref_alignment}'
        LOGGER.printLog(f'mafft command: {ANSI.GREEN}{command}{ANSI.END}')

    if TIMEIT:
        start = perf_counter_ns()
    ret = os.system(command)
    if TIMEIT:
        end = perf_counter_ns()
    if ret != 0:
        LOGGER.error(f'Error in building alignment with mafft with error code {ret}', error_type=12)
    if TIMEIT:
        LOGGER.printLog(f'{ANSI.YELLOW}TIMED: mafft took {pd.to_timedelta(end-start)}, {end - start} nanoseconds{ANSI.END}')
    
    return ref_alignment

def getMapping(alignment_path : str, outpath : str, first_label : str, second_label : str):
    '''
    Get base mapping from reference to reference, split indels and substitutions/matches
    Only works on pairwise alignmend (max 2 sequences).

    Returns
    -------
    aligned_positions : dict
        {pos_1 : (pos_2, alignment_pos)}, stores a mapping between the bases of sequence 1 (pos_1) to the bases of sequence 2 (pos_2). Additionally stores the alignment position.
    unaligned_positions : dict
        {seq_id : [pos]}, stores a list of unaligned positions for each reference
    sequences : dict
        {seq_id : seq}, store the reference string for each reference
    '''
    fasta = SeqIO.parse(open(alignment_path), 'fasta')
    sequences = {}
    outfile = os.path.join(outpath, 'alignment', f'{first_label}_{second_label}_refdiffs.csv')
    w = open(outfile, 'w')
    w.write(f'type,{first_label}_pos,{second_label}_pos,base_{first_label},base_{second_label},alignment_motif_{first_label},alignment_motif_{second_label}\n')

    for seq in fasta:
        sequences[seq.id] = str(seq.seq)
    
    LOGGER.printLog(f'Found an alignment for sequences {list(sequences.keys())}')
    
    # {(pos_of_first_sample, base) : (pos_of_second_sample, base)}
    # different reference files for both samples
    if not 0 < len(sequences) < 3:
        LOGGER.error(f'Number of provided reference sequences is not equal 1 or 2: {len(sequences)}', 15)
    
    elif len(sequences) == 2:
        seq1_id, seq2_id = sequences.keys()
        seq1_seq, seq2_seq = sequences.values()
        i1 = i2 = 0
        first2sec_refpos_mapping = {}
        unaligned_positions = {seq : [] for seq in list(sequences.keys())}
        
        for alip, (base1, base2) in enumerate(zip(seq1_seq, seq2_seq)): # (seq1, base1), (seq2, base2) in zip(sequences.items()):
            motif1 = seq1_seq[max(0,alip-3):min(len(seq1_seq),alip+4)].upper()
            motif2 = seq2_seq[max(0,alip-3):min(len(seq2_seq),alip+4)].upper()
            base1 = base1.upper()
            base2 = base2.upper()

            if base1 == '-':
                unaligned_positions[seq2_id].append((i2, base2))
                w.write(f'insert_{second_label},{i1},{i2},{base1},{base2},{motif1},{motif2}\n')
                i2 += 1
                
            elif base2 == '-':
                unaligned_positions[seq1_id].append((i1, base1))
                w.write(f'insert_{first_label},{i1},{i2},{base1},{base2},{motif1},{motif2}\n')
                i1 += 1
            
            else:
                first2sec_refpos_mapping[i1] = (i2, alip)
                if base1 == 'N' or base2 == 'N':
                    w.write(f'N,{i1},{i2},{base1},{base2},{motif1},{motif2}\n')
                elif base1 != base2:
                    if base2 in IUPAC[base1] or base1 in IUPAC[base2]: # could be a mutation, could not be
                        w.write(f'possible_subsitution,{i1},{i2},{base1},{base2},{motif1},{motif2}\n')
                    else: # mutation
                        w.write(f'substitution,{i1},{i2},{base1},{base2},{motif1},{motif2}\n')
                i1 += 1
                i2 += 1

        assert i1 == len(seq1_seq.replace('-', '')), f'Iterator mismatch: {i1} != {len(seq1_seq.replace("-", ""))}\t{seq1_seq}'
        assert i2 == len(seq2_seq.replace('-', '')), f'Iterator mismatch: {i2} != {len(seq2_seq.replace("-", ""))}\t{seq2_seq}'

        LOGGER.printLog(f'Found {len(first2sec_refpos_mapping)} aligned positions, {len(unaligned_positions[seq1_id])} unaligned positions for sequence {seq1_id}, {len(unaligned_positions[seq2_id])} unaligned positions for sequence {seq2_id}')

        return first2sec_refpos_mapping, unaligned_positions, sequences

    # case: same reference fasta for both samples
    else:
        seq_id = list(sequences.keys())[0]
        seq = list(sequences.values())[0]

        first2sec_refpos_mapping = {i:(i,i) for i in range(len(seq))}
        unaligned_positions = {seq_id : []}

        return first2sec_refpos_mapping, unaligned_positions, sequences

def readRedFile(red_file : str, seq_dict : dict):
    '''
    Reads the RED file, stores data in a dictionary

    Returns
    -------
    red_sequences : dict
        read data from RED files
        {seq_id:[pos, strand, features]}
    '''
    LOGGER.printLog(f'Reading RED file {red_file}')
    # sequences are stored as {reference: [pos, strand, REDENCODER]}
    red_sequences = {}

    with open(red_file, 'r') as red:
        for line in red:
            if line.startswith('reference\t'):
                continue

            # w.write('reference\tposition\tstrand\tsignal_mean\tsignal_std\tdata_density\texpected_model_density\tn_datapoints\tcontained_datapoints\tn_segments\tcontained_segments\tn_reads\n')
            ref_id, pos, strand, mean, std, data_density, expected_model_density, n_datapoints, contained_datapoints, n_segments, contained_segments, n_reads = line.strip().split('\t')
            pos = int(pos)

            if ref_id not in red_sequences:
                red_sequences[ref_id] = np.zeros((len(seq_dict[ref_id]), len(STRANDENCODER), len(REDENCODER)), dtype=float)

            # Order of REDENCODER + expected_model_density at the end
            red_sequences[ref_id][pos, STRANDENCODER[strand]] = [mean, std, data_density, n_datapoints, contained_datapoints, n_segments, contained_segments, n_reads, expected_model_density]
    return red_sequences

def asyncCompareSignals(strand : int, base1 : str, base2 : str, motif1 : str, motif2 : str, alip : int, mut_context : bool, data_pos1 : np.ndarray, data_pos2 : np.ndarray, seqs_ids : tuple, pos1 : int, pos2 : int, sidx : int, max_sidx : int) -> None:
    global num_muts, sign_pos, no_data, low_cov_count, allQueue, magnQueue, alipQueue

    m1 = data_pos1[strand, REDENCODER['mean']]
    s1 = data_pos1[strand, REDENCODER['std']]
    m2 = data_pos2[strand, REDENCODER['mean']]
    s2 = data_pos2[strand, REDENCODER['std']]
    # check if positions have a distribution -> stdev for both are non-0
    hasData = s1 and s2
    low_cov = data_pos1[strand, REDENCODER['n_reads']] < 10 or data_pos2[strand, REDENCODER['n_reads']] < 10

    if strand == 1:
        base1 = complement(base1)
        motif1 = rev_complement(motif1)
        base2 = complement(base2)
        motif2 = rev_complement(motif2)

    bayesian_p = np.nan
    kl_divergence = np.nan
    td = np.nan
    significant = False

    if hasData:
        dist1 = NormalDist(m1, s1)
        dist2 = NormalDist(m2, s2)
        mDiff = abs(m1 - m2)
        sAvg = (s1 + s2)/2
        bayesian_p = dist1.overlap(dist2)
        kl_divergence = kullback_leibler_normal(m1, s1, m2, s2)
        td = td_score(mDiff, sAvg)
        significant = td>=1

    outline = f'{STRANDDECODER[strand]}\t{td:.8f}\t{kl_divergence:.8f}\t{bayesian_p:.8f}\t{MUTDECODER[mut_context]}\t{seqs_ids[0]}\t{pos1}\t{base1}\t{motif1}\t{m1:.8f}\t{s1:.8f}\t{data_pos1[strand, REDENCODER["n_datapoints"]]:.0f}\t{data_pos1[strand, REDENCODER["contained_datapoints"]]:.0f}\t{data_pos1[strand, REDENCODER["n_segments"]]:.0f}\t{data_pos1[strand, REDENCODER["contained_segments"]]:.0f}\t{data_pos1[strand, REDENCODER["n_reads"]]:.0f}\t{seqs_ids[1]}\t{pos2}\t{base2}\t{motif2}\t{m2:.8f}\t{s2:.8f}\t{data_pos2[strand, REDENCODER["n_datapoints"]]:.0f}\t{data_pos2[strand, REDENCODER["contained_datapoints"]]:.0f}\t{data_pos2[strand, REDENCODER["n_segments"]]:.0f}\t{data_pos2[strand, REDENCODER["contained_segments"]]:.0f}\t{data_pos2[strand, REDENCODER["n_reads"]]:.0f}\n'

    allQueue.put((sidx, max_sidx, outline))

    if significant:
        magnQueue.put(outline)
        alipQueue.put(alip)
        with num_muts.get_lock(): num_muts.value += mut_context
        with sign_pos.get_lock(): sign_pos.value += 1

    with no_data.get_lock(): no_data.value += 1 - int(hasData)
    with low_cov_count.get_lock(): low_cov_count.value += low_cov
    with num_pos.get_lock(): num_pos.value += 1

def asyncWriter_withPrint(file : str, q : mp.Queue) -> None:
    with open(file, 'a') as f:
        while 1:
            entry = q.get()
            if entry is None:
                break
            lidx, max_lidx, line = entry
            if (lidx + 1) % 1000 == 0:
                print(f'\t{lidx + 1}/{max_lidx}', end='\r')
            f.write(line)
            f.flush()

def asyncWriter(file : str, q : mp.Queue) -> None:
    with open(file, 'a') as f:
        while 1:
            line = q.get()
            if line is None:
                break
            f.write(line)
            f.flush()

def asyncStockholmWriter(stkStrings : list, q : mp.Queue, alignment_path : str, first_sample_label : str, sec_sample_label : str, outdir : str) -> None:
    while 1:
        pos = q.get()
        if pos is None:
            break
        for stk in stkStrings:
            stk[pos] = 'X'
    
    LOGGER.printLog(f'Writing stockholm file containing magnipore results')
    records = []
    ids = []
    # remove duplicate sequences (happens when both samples have the same reference fasta)
    for record in SeqIO.parse(open(alignment_path), 'fasta'):
        if record.id not in ids:
            records.append(record)
            ids.append(record.id)
    records.extend([
        SeqIO.SeqRecord(
        Seq.Seq((''.join(seq)).upper()),
        id ='magnipore_marked_' + str((first_sample_label, sec_sample_label)[i]),
        name ='magnipore_' + str((first_sample_label, sec_sample_label)[i]),
        description='X=significant signal change, .=not significant') for i, seq in enumerate(stkStrings)
    ])
    SeqIO.write(records, os.path.join(outdir, first_sample_label + '_' + sec_sample_label + '_marked.stk'), 'stockholm')

def magnipore(mapping : dict, unaligned : dict, seq_dict : dict, aln_dict: dict, alignment_path : str, red1 : dict, red2 : dict, first_sample_label : str, sec_sample_label : str, working_dir : str, pore_type : str, processes : int) -> None:
    
    working_dir = os.path.join(working_dir, 'magnipore', f'{first_sample_label}_{sec_sample_label}')
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    magnipore_sign_file = os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}.magnipore')
    magnipore_all_file = os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}.all')
    with open(magnipore_sign_file, 'w') as w:
        w.write('\t'.join(MAGNIPORE_COLUMNS) + '\n')
    with open(magnipore_all_file, 'w') as w:
        w.write('\t'.join(MAGNIPORE_COLUMNS) + '\n')

    # setup multiprocessing
    global num_muts, sign_pos, no_data, low_cov_count, allQueue, magnQueue, alipQueue
    processes = max(4, processes)
    pool = mp.Pool(processes-3)
    allWriter = mp.Process(target=asyncWriter_withPrint, args=(magnipore_all_file, allQueue))
    magnWriter = mp.Process(target=asyncWriter, args=(magnipore_sign_file, magnQueue))
    reformat = lambda seq : list(re.sub(r"[^-]", ".", seq))
    magnipore_strings = list(map(reformat, aln_dict.values()))
    stockholmWriter = mp.Process(target=asyncStockholmWriter, args=(magnipore_strings, alipQueue, alignment_path, first_sample_label, sec_sample_label, working_dir))
    allWriter.start()
    magnWriter.start()
    stockholmWriter.start()

    seqs_ids = list(seq_dict.keys())
    # when providing the same reference for both samples seqs_ids is only length 1, add the same sequence id again
    if len(seqs_ids) == 1:
        seqs_ids.append(seqs_ids[0])
    seq1 = seq_dict[seqs_ids[0]].upper()
    seq2 = seq_dict[seqs_ids[1]].upper()

    if pore_type == 'r9':
        pore_range = 5
    elif pore_type == 'r10':
        pore_range = 9
    else:
        LOGGER.error('Unknown Pore Type', 16)
    num_indels = 0
    
    LOGGER.printLog(f'Start Signal comparison processes ({processes}) ...')
    # compare distributions of aligned positions
    for sidx, (pos1, (pos2, alip)) in enumerate(mapping.items()):
        
        # red: sequences are stored as {reference: {pos: {base: ('A'|'C'|'G'|'T'), mean: float, std: float}}}
        data_pos1 = red1[seqs_ids[0]][pos1]
        data_pos2 = red2[seqs_ids[1]][pos2]
        base1 = seq1[pos1]
        base2 = seq2[pos2]
        # no data for first and last two positions in red file, these positions should never be significant
        if not min(pos1, pos2) >= pore_range//2 and pos1+pore_range//2+1 <= len(seq1) and pos2+pore_range//2+1 <= len(seq2):
            continue
        motif1 = seq1[pos1-pore_range//2:pos1+pore_range//2+1]
        motif2 = seq2[pos2-pore_range//2:pos2+pore_range//2+1]
        mut_context = motif1 != motif2
        
        for strand in [0, 1]:
            pool.apply_async(asyncCompareSignals, args=(strand, base1, base2, motif1, motif2, alip, mut_context, data_pos1, data_pos2, seqs_ids, pos1, pos2, sidx, len(mapping)), error_callback=callbackError)

    LOGGER.printLog('Waiting for Signal comparison processes to finish ...')
    pool.close()
    pool.join()
    # Send termination signal to queues and collect writers
    allQueue.put(None)
    magnQueue.put(None)
    alipQueue.put(None)
    allWriter.join()
    magnWriter.join()
    stockholmWriter.join()
    # get counts
    with num_muts.get_lock(): num_muts_val = num_muts.value
    with sign_pos.get_lock(): sign_pos_val = sign_pos.value
    with no_data.get_lock(): no_data_val = no_data.value
    with low_cov_count.get_lock(): low_cov_count_val = low_cov_count.value

    LOGGER.printLog('Writing indels file')
    indelFile = os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}.indels')
    with open(indelFile, 'w') as f:
        f.write('type\tstrand\tref\tpos\tbase\n')
        for seq in unaligned:
            for position, base in unaligned[seq]:
                f.write(f'insert\t+\t{seq}\t{position}\t{base}\n')
                num_indels += 1

    print(f'\t{sidx + 1}/{len(mapping)}')
    LOGGER.printLog(f'Number of indels: {ANSI.YELLOW}{num_indels}{ANSI.END}\n'\
               f'Number of significant positions: {ANSI.YELLOW}{sign_pos_val}{ANSI.END} - Classified as mutations: {ANSI.YELLOW}{num_muts_val}{ANSI.END}\n'\
               f'Positions with no data {ANSI.YELLOW}{no_data_val}{ANSI.END}, at least one aligned position without information (no signals)\n'\
               f'Number of positions with low coverage in at least one sample: {ANSI.YELLOW}{low_cov_count_val}{ANSI.END} - I recommend filtering out these positions in the .magnipore file.\nPositions with now data or low coverage can be high if one strand has no aligned reads!\n'\
               f'Wrote {magnipore_sign_file}')

    with open(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}.txt'), 'w') as w:
        w.write(f'Number of indels: {num_indels}\n'\
                f'Number of significant positions: {sign_pos_val} - Classified as mutations: {num_muts_val}\n'\
                f'Positions with no data {no_data_val}, at least one aligned position without information (no signals)\n'\
                f'Number of positions with low coverage in at least one sample: {low_cov_count_val} - I recommend filtering out these positions in the .magnipore file.\nPositions with now data or low coverage can be high if one strand has no aligned reads!\n')

    return magnipore_all_file

def callbackError(error):
    LOGGER.error(f'Error in multiprocessing signal comparison: {error}', 17)

def ks_test(dist1 : tuple, dist2 : tuple) -> tuple:
    data1 = np.random.normal(dist1[0], dist1[1], 100)
    data2 = np.random.normal(dist2[0], dist2[1], 100)
    return ks_2samp(data1, data2)

def td_score(mDiff, sAvg) -> tuple:
    '''
    Calculates the td-score mDiff/sAvg

    Returns
    -------
    tdscore : float
    '''
    return mDiff/sAvg
    # return np.abs(mDiff - sAvg) / np.sqrt(2)

def kullback_leibler_normal(m0 : float, s0 : float, m1 : float, s1 : float) -> float:
    if not s1 or not s0: # if s0 or s1 are 0 cannot calculate kl divergence
        return np.nan
    return (np.square(s0/s1) + np.square(m1-m0)/np.square(s1) - 1 + np.log(np.square(s1)/np.square(s0))) / 2

def callNanosherlock(working_dir : str, sample_label : str, reference_path : str, fast5_path : str, basecalls_path : str, seq_sum_path : str, threads : int, mx : int, mk : int, guppy_bin : str, guppy_model : str, guppy_device : str, calculate_data_density : bool, rna : bool, r10 : bool, kmer_model : str, force_rebuild : bool) -> str:

    red_file_path = os.path.join(working_dir, 'magnipore', sample_label, f'{sample_label}.red')
    if not os.path.exists(red_file_path) or not os.path.exists(reference_path) or force_rebuild:
        command = f'python -m magnipore.nanosherlock {fast5_path} {reference_path} {working_dir} {sample_label} -t {threads} -mx {mx} -mk {mk} -e 1'

        if force_rebuild:
            command += ' --force_rebuild'
        if calculate_data_density:
            command += ' --calculate_data_density'
        if rna:
            command += ' -rna'
        if r10:
            command += ' -r10'
        if kmer_model is not None:
            command += f' -kmer {kmer_model}'
        if basecalls_path is not None:
            command += f' --basecalls {basecalls_path}'
            if seq_sum_path is not None:
                command += f' --sequencing_summary {seq_sum_path}'
        else:
            assert guppy_bin is not None and guppy_model is not None, 'Need at least the guppy binary path and model or path to basecalls'
            command += f' --guppy_bin {guppy_bin} --guppy_model {guppy_model} --guppy_device {guppy_device}'
        
        if TIMEIT:
            command = command + ' --timeit'
            start = perf_counter_ns()
        LOGGER.printLog(f'Pipeline command: {ANSI.GREEN}{command}{ANSI.END}')
        ret = os.system(command)
        if TIMEIT:
            end = perf_counter_ns()
        if ret != 0:
            LOGGER.error(f'Error in nanosherlock for sample {sample_label} with error code {ret}', error_type=13)
        if TIMEIT:
            LOGGER.printLog(f'{ANSI.YELLOW}TIMED: Calculating distributions of sample {sample_label} took {pd.to_timedelta(end-start)}, {end-start} nanoseconds{ANSI.END}')

    else:
        LOGGER.printLog(f'{sample_label} RED file already exists:\n-\t{red_file_path}')

    return red_file_path

def callMagniplot(magnipore_file : str, label_first_sample : str, label_sec_sample : str, threads : int) -> None:
    plot_dir = os.path.join(os.path.dirname(magnipore_file), 'plots')
    num_lines = num_pos.value + 1 # + header
    command = f'magniplot {magnipore_file} {plot_dir} {label_first_sample} {label_sec_sample} -t {threads} -nl {num_lines}'
    LOGGER.printLog(f'Creating plots from {magnipore_file}.')

    if TIMEIT:
        start = perf_counter_ns()
    LOGGER.printLog(f'Magniplot command: {ANSI.GREEN}{command}{ANSI.END}')
    ret = os.system(command)
    if TIMEIT:
        end = perf_counter_ns()
    if ret != 0:
        LOGGER.error(f'Error in magniplot with error code {ret}', 18)
    if TIMEIT:
        LOGGER.printLog(f'{ANSI.YELLOW}TIMED: Plotting took {pd.to_timedelta(end-start)}, {end-start} nanoseconds{ANSI.END}')

def main():
    
    args = parse()
    
    raw_data_first_sample = args.raw_data_first_sample
    ref_first_sample = args.reference_first_sample
    label_first_sample = args.label_first_sample
    
    raw_data_sec_sample = args.raw_data_sec_sample
    ref_sec_sample = args.reference_sec_sample
    label_sec_sample = args.label_sec_sample
    
    basecalls_first_sample = args.basecalls_first_sample
    basecalls_sec_sample = args.basecalls_sec_sample
    seqsum_first_sample = args.sequencing_summary_first_sample
    seqsum_sec_sample = args.sequencing_summary_sec_sample

    working_dir = args.working_dir
    guppy_bin = args.guppy_bin
    guppy_model = args.guppy_model
    guppy_device = args.guppy_device
    
    threads = args.threads
    force_rebuild = args.force_rebuild
    calculate_data_density = args.calculate_data_density
    rna = args.rna
    r10 = args.r10

    mx = args.minimap2x
    mk = args.minimap2k

    kmer_model = args.kmer_model

    global TIMEIT 
    TIMEIT = args.timeit

    global LOGGER
    log_file = os.path.join(working_dir, 'log', f'{label_first_sample}_{label_sec_sample}_magnipore.log')
    if not os.path.exists(os.path.join(working_dir, 'log')):
        os.makedirs(os.path.join(working_dir, 'log'))
    LOGGER = Logger(open(log_file, 'w'))
    
    # first sample
    red_first_sample = callNanosherlock(working_dir, label_first_sample, ref_first_sample, raw_data_first_sample, basecalls_first_sample, seqsum_first_sample, threads, mx, mk, guppy_bin, guppy_model, guppy_device, calculate_data_density, rna, r10, kmer_model, force_rebuild)

    # second sample
    red_sec_sample = callNanosherlock(working_dir, label_sec_sample, ref_sec_sample, raw_data_sec_sample, basecalls_sec_sample, seqsum_sec_sample, threads, mx, mk, guppy_bin, guppy_model, guppy_device, calculate_data_density, rna, r10, kmer_model, force_rebuild)

    # mafft alignment
    alignment_path = align(ref_first_sample, ref_sec_sample, label_first_sample, label_sec_sample, working_dir, threads, force_rebuild)
    mapping, unaligned, aln_dict = getMapping(alignment_path, working_dir, label_first_sample, label_sec_sample)
    seq_dict = {key:seq.replace('-', '') for key,seq in aln_dict.items()}

    # reading RED files
    red_first_sample = readRedFile(red_first_sample, seq_dict)
    red_sec_sample = readRedFile(red_sec_sample, seq_dict)
    
    if TIMEIT:
        start = perf_counter_ns()
    magnipore_all_file = magnipore(mapping, unaligned, seq_dict, aln_dict, alignment_path, red_first_sample, red_sec_sample, label_first_sample, label_sec_sample, working_dir, 'r10' if r10 else 'r9', threads)
    if TIMEIT:
        end = perf_counter_ns()
        LOGGER.printLog(f'{ANSI.YELLOW}TIMED: Evaluating distributions took {pd.to_timedelta(end-start)}, {end - start} nanoseconds{ANSI.END}')
    callMagniplot(magnipore_all_file, label_first_sample, label_sec_sample, threads)
    LOGGER.printLog('Done')

if __name__ == '__main__':
    main()
