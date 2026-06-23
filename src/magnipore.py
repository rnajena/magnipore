#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import multiprocessing as mp
from multiprocessing.sharedctypes import Synchronized
from multiprocessing.synchronize import Lock
from os import system
from os.path import join, dirname, basename, exists
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from statistics import NormalDist
import read5_ont
import numpy as np
from scipy.stats import median_abs_deviation as mad
from Bio import Seq, SeqIO
from scipy.stats import ks_2samp, norm
from pysam import AlignmentFile
from tqdm import tqdm
from src.Red import Red
from src.__init__ import __version_str__
from src.Helper import (ANSI, IUPAC, MAGNIPORE_COLUMNS, MUTDECODER,
                              STRANDDECODER, STRANDENCODER, PORE2K,
                              complement, rev_complement)
from src.Logger import Logger

# might wanna disable these for bug fixing
from warnings import simplefilter
np.seterr(divide='ignore', invalid='ignore')
simplefilter("ignore", category=RuntimeWarning)

LOGGER : Logger = None
COV_THRESHOLD = 10

def init_Logger(file : str = None):
    """
    Initialize the global logger.
    
    Args:
        file (str): Path to the file to write logs to. If None, no file is written.
    """
    global LOGGER
    LOGGER = Logger(file)

def callbackErrorRed(error):
    """
    Error callback for multiprocessing pool error in red building.
    
    Args:
        error: Exception raised by the worker process.
    """
    LOGGER.error(f'Error in multiprocessing red building: {error}', 3)

def callbackErrorComparison(error):
    """
    Error callback for multiprocessing pool error in signal comparison.
    
    Args:
        error: Exception raised by the worker process.
    """
    LOGGER.error(f'Error in multiprocessing magnipore signal comparison: {error}', 4)

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description='Required tools: see github https://github.com/JannesSP/magnipore',
        prog='magnipore run',
        )
    parser.add_argument('raw_data_first_sample', type = str, help='Path to POD5 file of first sample.')
    parser.add_argument('raw_data_sec_sample', type = str,  help='Path to POD5 file of second sample')
    parser.add_argument('basecalls_data_first_sample', type = str, help='Path to basecalled BAM file from Dorado of first sample.')
    parser.add_argument('basecalls_data_sec_sample', type = str,  help='Path to basecalled BAM file from Dorado of second sample')
    parser.add_argument('uncalled4_first_sample', type=str, help='Path to Uncalled4 output of first sample')
    parser.add_argument('uncalled4_sec_sample', type=str, help='Path to Uncalled4 output of second sample')
    parser.add_argument('alignment', type = str, help='Path to reference alignment file')
    parser.add_argument('outdir', type = str, help='Path to write all output files')
    parser.add_argument('pore',  type=str, choices=["rna_r9", "dna_r9", "rna_rp4", "dna_r10_260bps", "dna_r10_400bps"], help='Pore generation used to sequence the data')
    parser.add_argument('-l1', '--label_first_sample', type = str, default = 'sample_1', help='Name of the first sample')
    parser.add_argument('-l2', '--label_sec_sample', type = str, default='sample_2', help='Name of the second sample')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use')
    parser.add_argument('-d', '--calculate_data_density', action = 'store_true', default = False, help = 'Will calculate data density after building the models. Will increase runtime!')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s' + f' {__version_str__}')
    return parser.parse_args()

def getMapping(alignment : str, out : str, l1 : str, l2 : str) -> tuple[dict[int : tuple[int, int]], dict[str : list[tuple[int, str]]], dict[str : str], dict[str : str]]:
    """
    Reads a pairwise alignment file and returns a mapping of positions between the two sequences and a dictionary of unaligned positions for each sequence.
    
    Args:
        alignment (str): Path to the pairwise alignment file.
        out (str): Path to the output directory.
        l1 (str): Label for the first sequence.
        l2 (str): Label for the second sequence.
    
    Returns
    -------
    aligned_positions (dict[int \\: tuple[int, int]])
        A dictionary mapping positions of the first sequence to positions of the second sequence.
    unaligned_positions (dict[str \\: list[tuple[int, str]]])
        A dictionary of unaligned positions for each sequence.
    alignment (dict[str \\: str])
        A dictionary of the sequences with alignment gaps.
    sequences (dict[str \\: str])
        A dictionary of the sequences.
    """
    outfile = join(out, f'{l1}_{l2}_refdiffs.csv')

    with open(outfile, "w") as w:
        w.write(f'type,{l1}_pos,{l2}_pos,base_{l1},base_{l2},alignment_motif_{l1},alignment_motif_{l2}\n')

        alignment = {seq.id: str(seq.seq) for seq in SeqIO.parse(alignment, "fasta")}
        sequences = {id : seq.replace('-', '') for id, seq in alignment.items()}

        LOGGER.printLog(f'Found a pairwise alignment for {list(sequences.keys())}')

        # {(pos_of_first_sample, base) : (pos_of_second_sample, base)}
        num_seqs = len(alignment)
        if not 0 < num_seqs < 3:
            LOGGER.error(f'Number of provided reference sequences is not equal 1 or 2: {num_seqs}', 2)
        
        if num_seqs == 1:
            seq_id, seq = next(iter(sequences.items()))
            return {i: (i, i) for i in range(len(seq))}, {seq_id: []}, alignment, sequences

        seq1_id, seq2_id = alignment.keys()
        seq1_seq, seq2_seq = map(str.upper, alignment.values())  # Convert both sequences to uppercase
        i1 = i2 = 0
        mapping = {}
        unaligned = {seq1_id: [], seq2_id: []}
        
        for alip, (base1, base2) in enumerate(zip(seq1_seq, seq2_seq)):
            motif1 = seq1_seq[max(0, alip-3) : alip+4]
            motif2 = seq2_seq[max(0, alip-3) : alip+4]

            if base1 == '-':
                unaligned[seq2_id].append((i2, base2))
                w.write(f'insert_{l2},{i1},{i2},{base1},{base2},{motif1},{motif2}\n')
                i2 += 1
                
            elif base2 == '-':
                unaligned[seq1_id].append((i1, base1))
                w.write(f'insert_{l1},{i1},{i2},{base1},{base2},{motif1},{motif2}\n')
                i1 += 1
            
            else:
                mapping[i1] = (i2, alip)
                if 'N' in (base1, base2):
                    w.write(f'N,{i1},{i2},{base1},{base2},{motif1},{motif2}\n')
                elif base1 != base2:
                    mutation_type = (
                        "possible_substitution"
                        if base2 in IUPAC.get(base1, "") or base1 in IUPAC.get(base2, "")
                        else "substitution"
                    )
                    w.write(f'{mutation_type},{i1},{i2},{base1},{base2},{motif1},{motif2}\n')

                i1 += 1
                i2 += 1

        assert i1 == len(seq1_seq.replace("-", "")), f'Iterator mismatch: {i1} != {len(seq1_seq.replace("-", ""))}'
        assert i2 == len(seq2_seq.replace("-", "")), f'Iterator mismatch: {i2} != {len(seq2_seq.replace("-", ""))}'

        LOGGER.printLog(
            f'Found {len(mapping)} aligned positions, '
            f'{len(unaligned[seq1_id])} unaligned positions for sequence {seq1_id}, '
            f'{len(unaligned[seq2_id])} unaligned positions for sequence {seq2_id}'
        )

    return mapping, unaligned, alignment, sequences

def get_readid_map(basecalls : str) -> dict:
    """
    Creates a mapping of read IDs to their processed identifiers from a basecalled BAM or SAM file.

    Args:
        basecalls (str): Path to the basecalled BAM or SAM file.

    Returns:
        dict: A dictionary mapping each read's query name to its processed identifier (if available)
              or the query name itself if the processed identifier is not present.
    """
    read_id_map = {}
    with AlignmentFile(basecalls, "r" if basecalls.endswith('.sam') else "rb", check_sq=False) as samfile:
        for basecalled_read in samfile.fetch(until_eof=True):
            read_id_map[basecalled_read.query_name] = basecalled_read.get_tag("pi") if basecalled_read.has_tag("pi") else basecalled_read.query_name
    return read_id_map

def checker_task(queue : mp.Queue, reds : list[list[Red]], num_updaters : int, raw : str, read_id_map : dict, cal_data_density : bool, processed_counter, lock) -> list[list[Red]]:
    """
    Updater process: Reads from queue and updates shared `reds` in-place.
    """
    r5 = read5_ont.read(raw)
    oldid = None
    signal = None
    local_counter = 0  # Local counter to batch updates

    while 1:
        line = queue.get()
        if line is None:
            r5.close()
            with lock:
                processed_counter.value += local_counter
            return reds

        try:
            _, pos, strand, _, _, readid, start, length, _, _ = line.split('\t')
            pos, strand = int(pos), STRANDENCODER[strand]
            start, length = int(start), int(length)
        except ValueError: # Handle unaligned kmers
            continue

        signalid = read_id_map[readid]
        if oldid != signalid:
            signal = r5.getZNormSignal(signalid)
            oldid = signalid
        segment = signal[start : start + length]

        red = reds[pos // num_updaters][strand]

        mean, stdev = red.get_signal_mean_stdev()
        r = 3 * stdev  # ~99% density of normal distribution
        lower_bound, upper_bound = mean - r, mean + r

        contained = np.count_nonzero((segment >= lower_bound) & (segment <= upper_bound))
        if contained == len(segment):
            red.add_contained_segments(1)
        red.add_contained_datapoints(contained)
        
        if cal_data_density:
            red.add_data_density(np.mean(norm.pdf(segment, loc = mean, scale = stdev)))

        # Increment the local counter
        local_counter += 1

        # Batch update the processed counter
        if local_counter >= 1024:  # Adjust batch size as needed
            with lock:
                processed_counter.value += local_counter
            local_counter = 0

    r5.close()
    return reds

def updater_task(queue : mp.Queue, reds : list[list[Red]], num_updaters : int, raw : str, read_id_map : dict, processed_counter, lock) -> list[list[Red]]:
    """
    Reads line from the queue and updates RED objects.
    """
    r5 = read5_ont.read(raw)
    oldid = None
    signal = None
    local_counter = 0  # Local counter to batch updates

    while 1:
        line = queue.get()
        if line is None:
            r5.close()
            with lock:
                processed_counter.value += local_counter
            return reds

        try:
            #! uncalled4
            _, pos, strand, _, _, readid, start, length, _, _ = line.split('\t')
            pos = int(pos)
            strand = STRANDENCODER[strand]
            start, length = int(start), int(length)

            #! dynamont
            # readid, signalid, start, end, basepos, base, motif, state, posterior_probability, polish = line.split(',')
            # pos = int(basepos)
            # strand = STRANDENCODER[strand]
            # start, length = int(start), int(end) - int(start)
        except ValueError: # Handle unaligned kmers
            continue

        signalid = read_id_map[readid]
        if oldid != signalid:
            signal = r5.getZNormSignal(signalid)
            oldid = signalid
        #? extract summary statistics instead of raw signal? mean or median, var or mad, for each segment
        segment = signal[start : start + length]
        entry = [[np.mean(segment), np.std(segment), np.median(segment), mad(segment)]]  # Calculate mean and std of the segment

        #! extract segment from signal
        # segment = signal[start : start + length]
        
        # update RED object
        red = reds[pos // num_updaters][strand]
        red.add_reads(1)
        red.append(entry)
        red.add_datapoints(len(segment))
        red.add_segments(1)

        # Increment the local counter
        local_counter += 1

        # Batch update the processed counter
        if local_counter >= 1024:  # Adjust batch size as needed
            with lock:
                processed_counter.value += local_counter
            local_counter = 0

    r5.close()
    return reds

def progress_updater(progress : tqdm, processed_counter : Synchronized, lock : Lock):
    """
    Continuously updates the progress bar based on the processed count.

    Parameters
    ----------
    progress : tqdm
        A tqdm progress bar instance to update.
    processed_counter : Synchronized[int]
        A synchronized integer representing the number of processed items.
    lock : Lock
        A lock to ensure thread-safe access to the processed counter.
    """
    last_count = 0
    while True:
        with lock:
            new_count = processed_counter.value
        if new_count < 0:
            break
        progress.update(new_count - last_count)
        last_count = new_count

def reconstruct_reds(red_chunks : list[list[Red]], red_len : int, t : int) -> list[list[Red]]:
    """
    Reconstructs a list of Red objects from divided chunks.

    Parameters
    ----------
    red_chunks : list[list[Red]]
        A nested list where each sublist contains a portion of Red objects.
    red_len : int
        The total length of the reconstructed list.
    t : int
        The number of chunks or segments the list is divided into.

    Returns
    -------
    list[list[Red]]
        A reconstructed nested list of Red objects, reassembled from the chunks.
    """
    return [red_chunks[pos % t][pos // t] for pos in range(red_len)]

def build_red_objects(reds : list[list[Red]], raw : str, read_id_map : dict, segmentation : str, calculate_data_density : bool, t : int, max_lines=None) -> list[list[Red]]:
    """
    Spawns worker processes for data extraction and a single updater process for `reds` updates.
    """
    LOGGER.printLog(f"Using {t} worker processes...")
    manager = mp.Manager()
    lock = manager.Lock()

    # Shared progress counter
    processed_counter = manager.Value('i', 0) 
    # Start a background process to update the progress bar
    progress = tqdm(total=max_lines if max_lines else None, desc="Building Models", unit=" lines", initial=0, leave=False)
    progress_process = mp.Process(target=progress_updater, args=(progress, processed_counter, lock))
    progress_process.start()

    queues = [manager.Queue() for _ in range(t)]
    reds_chunks = np.array_split(reds, t) # Split `reds` evenly across updaters
    worker_pool = mp.Pool(t)
    worker = [ # Start updater processes, each handling a subset of `reds`
        worker_pool.apply_async(updater_task, (queues[i], reds_chunks[i], t, raw, read_id_map, processed_counter, lock), error_callback=callbackErrorRed)
        for i in range(t)
    ]

    # LOGGER.printLog(f"Start reading {segmentation}")
    with open(segmentation, 'r') as file:
        next(file)  # Skip header

        # Read file line-by-line without loading everything in memory
        for line_count, line in enumerate(file, 1):

            #! uncalled4
            pos = int(line.split('\t')[1])  # Extract position from line
            #! dynamont
            # pos = int(line.split(',')[4])  # Extract position from line

            queues[pos % t].put(line)
            
            if max_lines and line_count >= max_lines:
                break

    # LOGGER.printLog("Done reading file, waiting for updaters to finish...", newline_before=True)
    # Signal updaters to stop processing
    for queue in queues:
        queue.put(None)

    # Collect updated `reds` from all updaters
    reds_chunks = [updater.get() for updater in worker]

    worker_pool.close()
    worker_pool.join()

    with lock:
        total = processed_counter.value
        processed_counter.value = -1
    progress.update(total)
    # Ensure progress is updated fully before closing
    progress_process.join()
    progress.close()
    LOGGER.printLog("Finished data collection")

    ###############! check data density and contained segments

    LOGGER.printLog("Calculating additional statistics")

    # reset progress for new progress bar
    with lock:
        processed_counter.value = 0
    # Start a background process to update the progress bar
    progress = tqdm(total=total, desc="Calculate Statistics", unit=" lines", initial=0, leave=False)
    progress_process = mp.Process(target=progress_updater, args=(progress, processed_counter, lock))
    progress_process.start()

    # Start updater processes, each handling a subset of `reds`
    worker_pool = mp.Pool(t)
    worker = [
        worker_pool.apply_async(checker_task, (queues[i], reds_chunks[i], t, raw, read_id_map, calculate_data_density, processed_counter, lock), error_callback=callbackErrorRed)
        for i in range(t)
    ]

    # Start worker pool
    with open(segmentation, 'r') as file:
        next(file)  # Skip header

        # Read file line-by-line without loading everything in memory
        for line_count, line in enumerate(file, 1):
            pos = int(line.split("\t")[1])  # Extract position from line
            queues[pos % t].put(line)

            if max_lines and line_count >= max_lines:
                break

    # Signal updaters to stop processing
    for queue in queues:
        queue.put(None)

    # Collect updated `reds` from all updaters
    reds_chunks = [updater.get() for updater in worker]
    # Merge updated results
    reds = reconstruct_reds(reds_chunks, len(reds), t)
    worker_pool.close()
    worker_pool.join()

    with lock:
        # total = processed_counter.value
        processed_counter.value = -1
    progress.update(total)
    # Ensure progress is updated fully before closing
    progress_process.join()
    progress.close()

    return reds

def write_red_file(red_file : str, reds : list[list[Red]]):
    """
    Writes RED data to a file.

    Parameters
    ----------
    red_file : str
        path to the output RED file
    reds : list[list[Red]]
        stores the data that will be written to the RED files
    """
    total_entries = sum(len(row) for row in reds)  # Total elements for tqdm
    with open(red_file, 'w') as file:
        # dwell_time_mean\tdwell_time_std\t
        file.write('strand\tposition\tsignal_mean\tsignal_std\tdata_density\texpected_model_density\tn_datapoints\tcontained_datapoints\tn_segments\tcontained_segments\tn_reads\n')
        
        with tqdm(total=total_entries, desc="Writing RED file", unit=" entries", initial=1, leave=False) as progress:
            for pos, row in enumerate(reds):
                for strand, red in enumerate(row):
                    file.write(f'{STRANDDECODER[strand]}\t{pos}\t{red}\n')
                    progress.update(1)  # Update progress
            progress.close()

def pickle_reds(red_file: str, reds: list[list[Red]]):
    """
    Serializes and writes RED data to a binary file.

    Parameters
    ----------
    red_file : str
        Path to the output RED file.
    reds : list[list[Red]]
        Nested list of Red objects to be serialized and stored.
    """
    import pickle
    with open(red_file, 'wb') as file:
        pickle.dump(reds, file)

def read_red_file(red_file: str, seq: str) -> list[list[Red]]:
    """
    Reads a RED file and constructs a list of Red objects for each position and strand.

    Parameters
    ----------
    red_file : str
        Path to the RED file.
    seq : str
        Reference sequence.

    Returns
    -------
    reds : list[list[Red]]
        Nested list of Red objects with dimensions (len(seq), len(STRANDENCODER)).
    """
    LOGGER.printLog(f"Reading RED file {red_file}")

    # Initialize the list of lists with empty Red objects
    reds = [[Red() for _ in range(len(STRANDENCODER))] for _ in range(len(seq))]

    with open(red_file, "r") as file:
        next(file)  # Skip header

        for line in tqdm(file, desc="Reading RED file", leave=False, unit=" lines", initial=1):
            # strand  position        signal_mean     signal_std      data_density    expected_model_density  n_datapoints    contained_datapoints    n_segments      contained_segments      n_reads
            values = line.strip().split("\t")
            pos = int(values[1])
            strand_idx = STRANDENCODER[values[0]]

            # Extract feature values
            #  dwell_time_mean, dwell_time_std
            signal_mean, signal_std, data_density, expected_model_density, n_datapoints, contained_datapoints, n_segments, contained_segments, n_reads = map(float, values[2:])

            # Populate the Red object
            red = reds[pos][strand_idx]
            red.signal_stats.mean = signal_mean
            red.signal_stats.std = signal_std
            # red.dwell_time_stats.mean = dwell_time_mean
            # red.dwell_time_stats.std = dwell_time_std
            red.n = int(n_datapoints)
            red.data_density = data_density
            red.n_datapoints = int(n_datapoints)
            red.contained_datapoints = int(contained_datapoints)
            red.n_segments = int(n_segments)
            red.contained_segments = int(contained_segments)
            red.n_reads = int(n_reads)

    return reds

def nanosherlock(outdir : str, label : str, pod5 : str, bam : str, uncalled4 : str, calculate_data_density : bool, seq : str, t : int, max_lines : int = None) -> list[list[Red]]:
    red_file = join(outdir, f'{label}.red')
    if exists(red_file):
        return read_red_file(red_file, seq)
    
    readidMap = get_readid_map(bam)
    LOGGER.printLog("Initiliazing RED models...")
    reds = [[Red(initlen=(1000, 4)) for _ in range(len(STRANDENCODER))] for _ in range(len(seq))]
    LOGGER.printLog(f"Updating RED models with {basename(uncalled4)}...")
    reds = build_red_objects(reds, pod5, readidMap, uncalled4, calculate_data_density, t, max_lines)
    LOGGER.printLog(f"Writing RED models to {basename(red_file)}...")
    # write_red_file(red_file, reds)

    red_file = join(outdir, f'{label}.pickle')
    pickle_reds(red_file, reds)

    return reds

def ks_test(dist1 : tuple, dist2 : tuple) -> tuple:
    """
    Perform a Kolmogorov-Smirnov test on two normal distributions.

    Parameters
    ----------
    dist1 : tuple
        (mean, std) of the first normal distribution
    dist2 : tuple
        (mean, std) of the second normal distribution

    Returns
    -------
    D : float
        The test statistic
    p : float
        The p-value of the test
    """
    return ks_2samp(
        np.random.normal(*dist1, 100),
        np.random.normal(*dist2, 100)
    )

# def cohens_d(mu1 : float, mu2 : float, s1 : float, s2 : float) -> float:
#     """
#     Calculate Cohen's d effect size.

#     Parameters
#     ----------
#     mu1 : float
#         Mean of the first group
#     mu2 : float
#         Mean of the second group
#     s1 : float
#         Standard deviation of the first group
#     s2 : float
#         Standard deviation of the second group

#     Returns
#     -------
#     d : float
#         Cohen's d effect size
#     """
#     return (mu1 - mu2) / np.sqrt((s1 ** 2 + s2 ** 2) / 2)

def cohens_d(mu1 : float, mu2 : float, s1 : float, s2 : float, n1 : int, n2 : int):
    """
    Calculate Cohen's d effect size.

    Parameters
    ----------
    mu1 : float
        Mean of the first group
    mu2 : float
        Mean of the second group
    s1 : float
        Standard deviation of the first group
    s2 : float
        Standard deviation of the second group

    Returns
    -------
    d : float
        Cohen's d effect size
    """
    # https://novustat.com/statistik-blog/cohens-d-effektstaerke-berechnen.html
    if (n1 + n2 <= 2) or (n1 <= 1) or (n2 <= 1) or (s1 == 0 and s2 == 0):
        return 0
    denominator = np.sqrt(((n1 - 1) * s1 ** 2 + (n2 - 1) * s2 ** 2) / (n1 + n2 - 2))
    if denominator == 0:
        return 0
    return (mu1 - mu2) / denominator

def is_significant(cd : float) -> bool:
    """
    Check if Cohen's d is significant.

    Parameters
    ----------
    cd : float
        Cohen's d effect size

    Returns
    -------
    bool
        True if significant, False otherwise
    """
    return abs(cd) >= 0.8  # Common threshold for large effect size

# def td_score(mDiff, sAvg) -> tuple:
#     '''
#     Calculates the td-score mDiff/sAvg

#     Returns
#     -------
#     tdscore : float
#     '''
#     return mDiff/sAvg
#     # return np.abs(mDiff - sAvg) / np.sqrt(2)

def kullback_leibler_normal(m0 : float, s0 : float, m1 : float, s1 : float) -> float:
    if not s1 or not s0: # if s0 or s1 are 0 cannot calculate kl divergence
        return np.nan
    ratio = (s0 / s1) ** 2
    return (ratio + (m1-m0) ** 2/s1 ** 2 - 1 + np.log(ratio)) / 2

def compare_signals(args):
    strand, base1, base2, motif1, motif2, alip, data_pos1, data_pos2, seqs_ids, pos1, pos2, num_muts, sign_pos, no_data, low_cov_count, num_pos, lock, all_queue, sign_queue, stk_queue = args
    m1, s1 = data_pos1.get_signal_mean_stdev()
    m2, s2 = data_pos2.get_signal_mean_stdev()
    # check if positions have a distribution -> stdev for both are non-0
    hasData = s1 and s2

    # skip positions where at least one sample has no data
    if not hasData:
        with lock:
            no_data.value += 1
        return

    mut_context = motif1 != motif2

    if strand:
        base1, motif1, base2, motif2 = complement(base1), rev_complement(motif1), complement(base2), rev_complement(motif2)

    bayesian_p = NormalDist(m1, s1).overlap(NormalDist(m2, s2))
    kl_divergence = kullback_leibler_normal(m1, s1, m2, s2)
    # score = td_score(abs(m1 - m2), (s1 + s2) / 2)
    score = cohens_d(m1, m2, s1, s2, data_pos1.n_reads, data_pos2.n_reads)
    significant = is_significant(score)

    outline = (
        f"{STRANDDECODER[strand]}\t{score:.4f}\t{kl_divergence:.4f}\t{bayesian_p:.4f}\t{MUTDECODER[mut_context]}\t"
        f"{seqs_ids[0]}\t{pos1}\t{base1}\t{motif1}\t" # data_pos1
        f"{m1:.7f}\t{s1:.7f}\t{data_pos1.magnipore_string()}\t"
        f"{seqs_ids[-1]}\t{pos2}\t{base2}\t{motif2}\t{m2:.7f}\t{s2:.7f}\t" # data_pos2
        f"{data_pos2.magnipore_string()}\n"
    )

    all_queue.put(outline)

    if significant:
        sign_queue.put(outline)
        stk_queue.put(alip)
        with lock:
            num_muts.value += mut_context
            sign_pos.value += 1

    with lock:
        low_cov_count.value += min(data_pos1.n_reads, data_pos2.n_reads) < COV_THRESHOLD
        num_pos.value += 1

def async_writer(file, queue):
    with open(file, 'a') as f:
        for line in iter(queue.get, None):
            f.write(line)

def async_stk_writer(stks : list, stk_queue : mp.Queue, alignment : dict, first_sample_label : str, sec_sample_label : str, outdir : str):
    while 1:
        pos = stk_queue.get()
        if pos is None:
            break
        for stk in stks:
            stk[pos] = 'X'
    
    LOGGER.printLog('Writing stockholm file containing magnipore results')
    # loading alignment sequences
    records = [
        SeqIO.SeqRecord(
            Seq.Seq(seq),
            id = id,
            name = id,
            description = f"alignment {id}"
        ) for id, seq in alignment.items()
    ]
    # Loading magnipore strings
    records.extend([
        SeqIO.SeqRecord(
            Seq.Seq((''.join(seq)).upper()),
            id ='magnipore_marked_' + label,
            name ='magnipore_' + label,
            description='X=significant signal change, .=not significant'
        ) for label, seq in zip((first_sample_label, sec_sample_label), stks)
    ])
    SeqIO.write(records, join(outdir, first_sample_label + '_' + sec_sample_label + '_marked.stk'), 'stockholm')

def reformat(seq):
    """
    Reformat a sequence to stockholm format by replacing all non-gap characters with dots.
    
    Parameters
    ----------
    seq : str
        The sequence to reformat
    
    Returns
    -------
    list
        A list of characters where all non-gaps are replaced with dots
    """
    from re import sub
    return list(sub(r"[^-]", ".", seq))  # replaces all non-gaps with dots for stockhol format

def magnipore(mapping : dict, unaligned : dict, seqs : dict[str : str], alignment : dict, red1 : list[list[Red]], red2 : list[list[Red]], l1 : str, l2 : str, outdir : str, pore_range : int, threads : int) -> tuple[str, int]:
    """
    Compares the distributions of aligned positions between two samples.

    Parameters
    ----------
    mapping : dict
        mapping of aligned positions between the two samples
    unaligned : dict
        dictionary of unaligned positions for each sample
    seqs : dict[str : str]
        dictionary of sequence names to sequences
    alignment : dict
        dictionary of alignment sequences
    red1 : list[list[Red]]
        list of Red objects for sample 1
    red2 : list[list[Red]]
        list of Red objects for sample 2
    l1 : str
        sample label 1
    l2 : str
        sample label 2
    outdir : str
        output directory for magnipore results
    pore_range : int
        number of positions to consider for each aligned position
    threads : int
        number of threads to use for parallel processing

    Returns
    -------
    tuple
        tuple containing the name of the output file and the number of positions written to the file
    """
    half_range = pore_range // 2
    sign_file = join(outdir, f'{l1}_{l2}.magnipore')
    all_file = join(outdir, f'{l1}_{l2}.all')
    with open(sign_file, 'w') as w:
        w.write('\t'.join(MAGNIPORE_COLUMNS) + '\n')
    with open(all_file, 'w') as w:
        w.write('\t'.join(MAGNIPORE_COLUMNS) + '\n')

    # setup multiprocessing
    manager = mp.Manager()
    pool = mp.Pool(threads)

    # initialize shared variables
    lock = manager.Lock()
    all_queue, sign_queue, stk_queue = manager.Queue(), manager.Queue(), manager.Queue()
    num_muts, sign_pos, no_data, low_cov_count, num_pos = [manager.Value('I', 0) for _ in range(5)]

    sign_writer = mp.Process(target=async_writer, args=(sign_file, sign_queue))
    all_writer = mp.Process(target=async_writer, args=(all_file, all_queue))
    magnipore_strings = list(map(reformat, alignment.values()))
    stk_writer = mp.Process(target=async_stk_writer, args=(magnipore_strings, stk_queue, alignment, l1, l2, outdir))
    all_writer.start()
    sign_writer.start()
    stk_writer.start()

    # when providing the same reference for both samples seqs_ids is only length 1, add the same sequence id again
    seq1, seq2 = seqs[list(seqs.keys())[0]], seqs[list(seqs.keys())[-1]]
    
    num_indels = 0
    
    LOGGER.printLog(f'Start signal comparison with {threads} processes')
    # compare distributions of aligned positions
    task_iter = (
        (
            strand,
            seq1[pos1],
            seq2[pos2],
            seq1[max(0, pos1 - half_range) : pos1 + half_range + 1],
            seq2[max(0, pos2 - half_range) : pos2 + half_range + 1],
            alip,
            red1[pos1][strand],
            red2[pos2][strand],
            list(seqs.keys()),
            pos1,
            pos2,
            num_muts,
            sign_pos,
            no_data,
            low_cov_count,
            num_pos,
            lock,
            all_queue,
            sign_queue,
            stk_queue
        )
        for pos1, (pos2, alip) in mapping.items()
        for strand in [0, 1]
    )
    
    # chunksize, speed (by eye and only correct for one system)
    # 32       , ~ 4000 pos/s
    # 64       , ~ 6000 pos/s
    # 128      , ~ 6000 pos/s
    # 256      , ~ 6000 pos/s
    # 512      , ~ 5000 pos/s
    # 1024     , ~ pos/s
    # Parallel processing with progress tracking
    list(tqdm(
        pool.imap_unordered(func = compare_signals, iterable = task_iter, chunksize = 128),
        total=len(mapping) * 2,
        desc="Comparing signals",
        leave=True,
        unit=" positions"
        ))

    pool.close()
    pool.join()
    # Send termination signal to queues and collect writers
    for q in [all_queue, sign_queue, stk_queue]:
        q.put(None)
    
    all_writer.join()
    sign_writer.join()
    stk_writer.join()
    
    # get counts
    with lock:
        num_muts_val = num_muts.value
        sign_pos_val = sign_pos.value
        no_data_val = no_data.value
        low_cov_count_val = low_cov_count.value

    indel_file = join(outdir, f'{l1}_{l2}.indels')
    LOGGER.printLog(f'Writing indel file to {dirname(indel_file)}')
    with open(indel_file, 'w') as f:
        f.write('type\tstrand\tref\tpos\tbase\n')
        for seq in unaligned:
            for position, base in unaligned[seq]:
                f.write(f'insert\t+\t{seq}\t{position}\t{base}\n')
                num_indels += 1
    
    LOGGER.printLog(f'Total alignment positions: {len(mapping)*2}\n'\
                    f'Indels: {ANSI.YELLOW}{num_indels}{ANSI.END}\n'\
                    f'Significant positions: {ANSI.YELLOW}{sign_pos_val}{ANSI.END}\n'\
                    f'Classified as mutations: {ANSI.YELLOW}{num_muts_val}{ANSI.END}\n'\
                    f'Positions with no data {ANSI.YELLOW}{no_data_val}{ANSI.END}, at least one sample at aligned position has no data\n'\
                    f'Positions with coverage < {COV_THRESHOLD} in at least one sample: {ANSI.YELLOW}{low_cov_count_val}{ANSI.END} - filtering recommended.\n'\
                     'Positions with no data or low coverage can be high if one strand has no aligned reads!\n'\
                    f'Wrote {sign_file}')

    with open(join(outdir, f'{l1}_{l2}.txt'), 'w') as w:
        w.write(f'Total alignment positions: {len(mapping)*2}\n'\
                f'Indels: {ANSI.YELLOW}{num_indels}{ANSI.END}\n'\
                f'Significant positions: {ANSI.YELLOW}{sign_pos_val}{ANSI.END}\n'\
                f'Classified as mutations: {ANSI.YELLOW}{num_muts_val}{ANSI.END}\n'\
                f'Positions with no data {ANSI.YELLOW}{no_data_val}{ANSI.END}, at least one sample at aligned position has no data\n'\
                f'Positions with coverage < {COV_THRESHOLD} in at least one sample: {ANSI.YELLOW}{low_cov_count_val}{ANSI.END} - filtering recommended.\n'\
                 'Positions with no data or low coverage can be high if one strand has no aligned reads!\n')

    with lock:
        return all_file, num_pos.value + 1

def call_magnipore_plot(magnipore_file : str, label_first_sample : str, label_sec_sample : str, threads : int, num_lines : int):
    plot_dir = join(dirname(magnipore_file), 'plots')
    command = f'magniplot {magnipore_file} {plot_dir} {label_first_sample} {label_sec_sample} -t {threads} -nl {num_lines}'
    LOGGER.printLog(f'Creating plots from {magnipore_file}.')

    LOGGER.printLog(f'Magniplot command: {ANSI.GREEN}{command}{ANSI.END}')
    ret = system(command)
    if ret != 0:
        LOGGER.error(f'Error in magniplot with error code {ret}', 5)

def main():
    
    init_Logger()

    args = parse()
    
    r1 = args.raw_data_first_sample
    r2 = args.raw_data_sec_sample
    b1 = args.basecalls_data_first_sample
    b2 = args.basecalls_data_sec_sample
    u1 = args.uncalled4_first_sample
    u2 = args.uncalled4_sec_sample
    al = args.alignment
    outdir = args.outdir
    l1 = args.label_first_sample
    l2 = args.label_sec_sample
    t = args.threads
    d = args.calculate_data_density

    k = PORE2K.get(args.pore, None)
    if k is None:
        LOGGER.error('Unknown Pore Type', 1)

    mapping, unaligned, alignment, sequences = getMapping(al, outdir, l1, l2)
    #! sequences must be in same order as input
    red1 = nanosherlock(outdir, l1, r1, b1, u1, d, list(sequences.values())[0], t)
    red2 = nanosherlock(outdir, l2, r2, b2, u2, d, list(sequences.values())[-1], t)
    
    # TODO remove later, just for testing reservoir right now
    exit(1000)
    magnipore_all_file, num_lines = magnipore(mapping, unaligned, sequences, alignment, red1, red2, l1, l2, outdir, k, t)

    call_magnipore_plot(magnipore_all_file, l1, l2, t, num_lines)
    LOGGER.printLog('Done')

if __name__ == '__main__':
    main()
