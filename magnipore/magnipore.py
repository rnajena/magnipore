#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import seaborn as sns
import os
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace

import numpy as np
from Bio import SeqIO, Seq
from matplotlib import pyplot as plt
from scipy.stats import ks_2samp
import pandas as pd
from statistics import NormalDist
import re
import matplotlib.lines as mlines
from time import perf_counter_ns
from magnipore.__init__ import __version__
from magnipore.Helper import ANSI, MAGNIPORE_COLUMNS, REDENCODER, STRANDENCODER, STRANDDECODER, IUPAC, complement, rev_complement
from magnipore.Logger import Logger

LOGGER : Logger = None
FONTSIZE = 18
TIMEIT = False
SUBSCRIPT = 'nanosherlock'

def initLogger(file = None) -> None:
    global LOGGER
    LOGGER = Logger(file)

def parse() -> Namespace:

    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description='Required tools: see github https://github.com/JannesSP/magnipore',
        prog='Magnipore',
        )
    
    parser.add_argument("path_to_fast5_first_sample", type = str, help='FAST5 file of first sample')
    parser.add_argument("path_to_reference_first_sample", type = str, help='reference FASTA file of first sample, POSITIVE (+) or FORWARD strand, ATTENTION: can only contain a single sequence')
    parser.add_argument("first_sample_label", type = str, help='Name of the sample or pipeline run')
    
    parser.add_argument("path_to_fast5_sec_sample", type = str, help='FAST5 file of second sample')
    parser.add_argument("path_to_reference_sec_sample", type = str, help='reference FASTA file of second sample, POSITIVE (+) or FORWARD strand, ATTENTION: can only contain a single sequence')
    parser.add_argument("sec_sample_label", type = str, help='Name of the sample or pipeline run')
    
    parser.add_argument("working_dir", type = str, help='Path to write all output files')
    
    parser.add_argument("--guppy_bin", type = str, default = None, help='Guppy binary')
    parser.add_argument("--guppy_model", type = str, default = None, help='Guppy model used for basecalling')
    parser.add_argument('--guppy_device', type=str, default='cuda:0', help='Use the GPU to basecall "cuda:0" to use the GPU with ID 0')

    parser.add_argument('--path_to_first_basecalls', metavar='FASTQ_DIR', type = str, default = None, help = 'Path to existing basecalls and sequencing summary file for first sample. Basecalls must be in one single file with the name <first_sample_label>.fastq')
    parser.add_argument('--path_to_sec_basecalls', metavar='FASTQ_DIR', type = str, default = None, help = 'Path to existing basecalls and sequencing summary file for second sample. Basecalls must be in one single file with the name <sec_sample_label>.fastq')
    parser.add_argument('--calculate_data_density', action = 'store_true', default = False, help = 'Will calculate data density after building the models. Will increase runtime!')

    parser.add_argument('-t', "--threads", type=int, default=1, help='Number of threads to use')
    parser.add_argument('-f5', '--fast5_out', action = 'store_true', help='Guppy generates FAST5 output (workspace folder) of Guppy')
    parser.add_argument('-fr', '--force_rebuild', action = 'store_true', help='Run commands regardless if files are already present')
    parser.add_argument('-mx', '--minimap2x', default = 'splice', choices = ['map-ont', 'splice', 'ava-ont'], help = '-x parameter for minimap2')
    parser.add_argument('-mk', '--minimap2k', default = 14, help = '-k parameter for minimap2')
    parser.add_argument('--timeit', default = False, action = 'store_true', help = 'Measure and print time used by submodules')

    parser.add_argument('-v', '--version', action='version', version='%(prog)s' + f' {__version__}')

    return parser.parse_args()

def mafft(ref_first_sample : str, ref_sec_sample : str, first_sample_label : str, sec_sample_label : str, working_dir : str, threads : int, force_rebuild : bool):
    
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
    
    # LOGGER.printLog(f'Building alignment in {ref_alignment}')
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
    '''
    
    fasta = SeqIO.parse(open(alignment_path), 'fasta')
    sequences = {}

    outfile = os.path.join(outpath, 'alignment', f'{first_label}_{second_label}_refdiffs.csv')
    w = open(outfile, 'w')
    w.write(f'type,{first_label}_pos,{second_label}_pos,base_{first_label},base_{second_label},alignment_motif_{first_label},alignment_motif_{second_label}\n')

    for seq in fasta:

        sequences[seq.id] = str(seq.seq)
        
    LOGGER.printLog(f'Found an alignment for sequences {list(sequences.keys())}')
    # LOGGER.printLog(f'Lengths are {list(map(len, sequences.values()))}')
    
    # {(pos_of_first_sample, base) : (pos_of_second_sample, base)}
    # different reference files for both samples
    try:
        seq1_id, seq2_id = sequences.keys()
        seq1_seq, seq2_seq = sequences.values()
        i1 = i2 = 0
        first2sec_refpos_mapping = {}
        unaligned_positions = {seq : [] for seq in list(sequences.keys())}
        
        for alip, (base1, base2) in enumerate(zip(seq1_seq, seq2_seq)): # (seq1, base1), (seq2, base2) in zip(sequences.items()):
            motif1 = seq1_seq[max(0,alip-3):min(len(seq1_seq),alip+4)].upper()
            motif2 = seq2_seq[max(0,alip-3):min(len(seq2_seq),alip+4)].upper()
            base1 = base2.upper()
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
                elif base2 not in IUPAC[base1]: # mutation
                    w.write(f'substitution,{i1},{i2},{base1},{base2},{motif1},{motif2}\n')
                i1 += 1
                i2 += 1
                
        assert i1 == len(seq1_seq.replace('-', '')), 'Mapping iterator for sample 1 does not match sequence length'
        assert i2 == len(seq2_seq.replace('-', '')), 'Mapping iterator for sample 2 does not match sequence length'

        LOGGER.printLog(f'Found {len(first2sec_refpos_mapping)} aligned positions, {len(unaligned_positions[seq1_id])} unaligned positions for sequence {seq1_id}, {len(unaligned_positions[seq2_id])} unaligned positions for sequence {seq2_id}')

        return first2sec_refpos_mapping, unaligned_positions, sequences

    # case: same reference fasta for both samples
    except:
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

# TODO adjust seq_ids for multiple references? -> how to align multiple segments/chromosomes between two samples?
def magnipore(mapping : dict, unaligned : dict, seq_dict : dict, aln_dict: dict, red_first_sample : dict, red_sec_sample : dict, first_sample_label : str, sec_sample_label : str, working_dir : str) -> tuple:

    seqs_ids = list(seq_dict.keys())
    working_dir = os.path.join(working_dir, 'magnipore', f'{first_sample_label}_{sec_sample_label}')

    if not os.path.exists(working_dir):
        os.makedirs(working_dir)

    magnipore_file = os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}.magnipore')
    magnipore = open(magnipore_file, 'w')
    magnipore.write('\t'.join(MAGNIPORE_COLUMNS) + '\n')
    all = open(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}.all'), 'w')
    all.write('\t'.join(MAGNIPORE_COLUMNS) + '\n')

    # red: sequences are stored as {reference: {pos: {base: ('A'|'C'|'G'|'T'), mean: float, std: float}}}
    num_indels, sign_pos, nans = 0, 0, 0
    # replace every nucleotide character with a dot
    reformat = lambda seq : list(re.sub(r"[^-]", ".", seq))
    # magnipore_strings = [list(re.sub(r"[^-]", ".", alignment_sequences[0])), list(re.sub(r"[^-]", ".", alignment_sequences[1]))]
    magnipore_strings = list(map(reformat, aln_dict.values()))

    # TODO add some quality value
    plotting_data = pd.DataFrame(columns=['mean_diff', 'first_std', 'sec_std', 'avg_std', 'mut_context', 'td_score', 'kl_divergence'])
    plotting_data = plotting_data.astype(
        {
            'mean_diff': 'float64',
            'first_std': 'float64',
            'sec_std': 'float64',
            'avg_std': 'float64',
            'mut_context': 'str',
            'td_score':'float64',
            'kl_divergence':'float64'
        })
    num_muts = 0

    LOGGER.printLog(f'Start comparing sequences position wise ...')
    first_seq = seq_dict[seqs_ids[0]].upper()
    sec_seq = seq_dict[seqs_ids[1]].upper()

    # compare distributions of aligned positions
    for sample_idx, (pos_first_sample, (pos_sec_sample, alip)) in enumerate(mapping.items()):
        
        if (sample_idx + 1) % 1000 == 0:
            print(f'\t{sample_idx + 1}/{len(mapping)}', end='\r')
        
        dist_first_sample = red_first_sample[seqs_ids[0]][pos_first_sample]
        dist_sec_sample = red_sec_sample[seqs_ids[1]][pos_sec_sample]
        first_base = first_seq[pos_first_sample]
        sec_base = sec_seq[pos_sec_sample]
        # no data for first and last two positions in red file, these positions should never be significant
        # motifs can have different lengths, rare case in start and end of reference
        if min(pos_first_sample, pos_sec_sample) >= 3 and pos_first_sample+r <= len(first_seq) and pos_sec_sample+r <= len(sec_seq):
            r = 3
        else:
            r = 2
        first_motif = first_seq[pos_first_sample-r:pos_first_sample+r+1]
        sec_motif = sec_seq[pos_sec_sample-r:pos_sec_sample+r+1]
        mut_context = first_motif != sec_motif
        
        for strand in [0, 1]:

            m0 = dist_first_sample[strand, REDENCODER['mean']]
            s0 = dist_first_sample[strand, REDENCODER['std']]
            m1 = dist_sec_sample[strand, REDENCODER['mean']]
            s1 = dist_sec_sample[strand, REDENCODER['std']]
            firstDist = NormalDist(m0, s0)
            secDist = NormalDist(m1, s1)
            mDiff = abs(m0 - m1)
            sAvg = (s0 + s1)/2
            # functions check if positions have a distribution -> stdev for both are non-0
            kl_divergence = kullback_leibler_normal(m0, s0, m1, s1)
            td, isnan = td_score(mDiff, sAvg)
            nans += isnan
            
            if strand == 1:
                first_base = complement(first_base)
                first_motif = rev_complement(first_motif)
                sec_base = complement(sec_base)
                sec_motif = rev_complement(sec_motif)

            new_entry = pd.DataFrame({
                        'first_mean': [m0],
                        'sec_mean': [m1],
                        'first_std' : [s0],
                        'sec_std' : [s1],
                        'mean_diff' : [mDiff],
                        'avg_std' : [sAvg],
                        'mut_context' : ['mutation' if mut_context else 'matching reference'],
                        'significant' : ['significant' if mDiff > sAvg else 'insignificant'],
                        'td_score' : [td],
                        'kl_divergence' : [kl_divergence]
                })
    
            plotting_data = pd.concat([plotting_data, new_entry], ignore_index=True)
            outline = f'{STRANDDECODER[strand]}\t{td:.8f}\t{kl_divergence:.8f}\t{firstDist.overlap(secDist) if s0 and s1 else np.nan:.8f}\t{"mut" if mut_context else "mod"}\t{seqs_ids[0]}\t{pos_first_sample}\t{first_base}\t{first_motif}\t{m0:.8f}\t{s0:.8f}\t{dist_first_sample[strand, REDENCODER["n_datapoints"]]:.0f}\t{dist_first_sample[strand, REDENCODER["contained_datapoints"]]:.0f}\t{dist_first_sample[strand, REDENCODER["n_segments"]]:.0f}\t{dist_first_sample[strand, REDENCODER["contained_segments"]]:.0f}\t{dist_first_sample[strand, REDENCODER["n_reads"]]:.0f}\t{seqs_ids[1]}\t{pos_sec_sample}\t{sec_base}\t{sec_motif}\t{m1:.8f}\t{s1:.8f}\t{dist_sec_sample[strand, REDENCODER["n_datapoints"]]:.0f}\t{dist_sec_sample[strand, REDENCODER["contained_datapoints"]]:.0f}\t{dist_sec_sample[strand, REDENCODER["n_segments"]]:.0f}\t{dist_sec_sample[strand, REDENCODER["contained_segments"]]:.0f}\t{dist_sec_sample[strand, REDENCODER["n_reads"]]:.0f}\n'
            all.write(outline)

            # distance between both means is greater than the average std of both distributions
            if mDiff > sAvg and s0 and s1:
                num_muts += mut_context
                # X == interesting position
                magnipore_strings[0][alip] = 'X'
                magnipore_strings[1][alip] = 'X'
                sign_pos += 1
                magnipore.write(outline)

    all.close()
    magnipore.close()

    with open(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}.indels'), 'w') as indels:
        indels.write(f'type\tstrand\tref\tpos\tbase\n')
        for seq in unaligned:
            for position, base in unaligned[seq]:                
                indels.write(f'insert\t+\t{seq}\t{position}\t{base}\n')
                num_indels += 1

    print(f'\t{sample_idx + 1}/{len(mapping)}')
    LOGGER.printLog(f'Number of indels: {num_indels}\n'\
               f'Number of significant positions: {sign_pos}\n'\
               f'Number of significant positions with reference differences: {num_muts}\n'\
               f'Number of nans {nans}, at least one aligned position without information (no signals)\nCan be high if one strand has no information!\n'\
               f'Wrote {magnipore_file}')

    with open(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}.txt'), 'w') as w:
        w.write(f'Number of indels: {num_indels}\n'\
                f'Number of significant positions: {sign_pos}\n'\
                f'Number of significant positions with reference differences: {num_muts}\n'\
                f'Number of nans {nans}, at least one aligned position without information (no signals)\nCan be high if one strand has no information!\n')

    LOGGER.printLog('Writing indels file')
    return plotting_data, magnipore_strings

def writeStockholm(magnipore_strings : list, alignment_path : str, first_sample_label : str, sec_sample_label : str, working_dir : str) -> None:
    LOGGER.printLog(f'Writing stockholm with magnipore markers')
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
        description='X=significant signal change, .=not significant') for i, seq in enumerate(magnipore_strings)
    ])
    SeqIO.write(records, os.path.join(working_dir, 'magnipore', f'{first_sample_label}_{sec_sample_label}', first_sample_label + '_' + sec_sample_label + '_marked.stk'), 'stockholm')
    LOGGER.printLog('Done with magnipore')

def plotStatistics(plotting_data : pd.DataFrame, working_dir : str, first_sample_label : str, sec_sample_label : str) -> None:
    plot_dir = os.path.join(working_dir, 'magnipore', f'{first_sample_label}_{sec_sample_label}', 'plots')
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    ### Mean Dist vs Std Avg plot
    LOGGER.printLog('Plotting Mean vs Stdev')
    plotMeanDiffStdAvg(plotting_data, plot_dir, first_sample_label, sec_sample_label)

    ### plot scores
    LOGGER.printLog(f'Plotting TD score and KL divergence')
    plotScores(plotting_data, plot_dir, first_sample_label, sec_sample_label)

def plotScores(dataframe : pd.DataFrame, working_dir : str, first_sample_label : str, sec_sample_label : str) -> None:
    
    colors = {
        'matching reference, insignificant':'wheat',
        'matching reference, significant':'darkorange',
        'mutation, insignificant':'skyblue',
        'mutation, significant':'darkblue'}
    dataframe['Context, Significance'] =  pd.Series(dataframe.reindex(['mut_context', 'significant'], axis='columns').astype('str').values.tolist()).str.join(', ')

    plt.figure(figsize = (12,8), dpi=300)
    plt.title(f'TD score for all positions\n{first_sample_label} vs. {sec_sample_label}')
    sns.histplot(data=dataframe, x='td_score', hue=dataframe['Context, Significance'], log_scale=(True, True), multiple="stack", palette=colors)
    plt.grid(True,  'both', 'both', alpha=0.6, linestyle='--')
    plt.tight_layout()
    plt.savefig(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}_td_score.png'))
    plt.savefig(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}_td_score.pdf'))
    plt.close()

    plt.figure(figsize = (12,8), dpi=300)
    plt.title(f'Kullback-Leibler divergence for all positions\n{first_sample_label} vs. {sec_sample_label}')
    sns.histplot(data=dataframe, x='kl_divergence', hue=dataframe['Context, Significance'], log_scale=(True, True), multiple="stack", palette=colors)
    plt.grid(True,  'both', 'both', alpha=0.6, linestyle='--')
    plt.tight_layout()
    plt.savefig(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}_kl_div.png'))
    plt.savefig(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}_kl_div.pdf'))
    plt.close()
    
def plotMeanDiffStdAvg(dataframe : pd.DataFrame, working_dir : str, first_sample_label : str, sec_sample_label : str) -> None:
    
    marker = lambda mut_context: 'D' if mut_context == 'mutation' else 'o'
    color = lambda mut_context: 'blue' if mut_context == 'mutation' else '#d95f02' 

    ### Mean Dist vs Std Avg plot
    plt.figure(figsize = (12,12), dpi=300)
    plt.rcParams.update({
        'font.size': FONTSIZE,
        })
    label1 = first_sample_label.replace("_", " ")
    label2 = sec_sample_label.replace("_", " ")

    g = sns.JointGrid(x = 'mean_diff', y = 'avg_std', data = dataframe, hue = 'mut_context', marginal_ticks=True, palette=['blue', '#d95f02'], hue_order=['mutation', 'matching reference'], height = 10)
    g.plot_joint(func=sns.scatterplot, s = 8)
    g.ax_joint.cla()
    for _, row in dataframe.iterrows():
        g.ax_joint.plot(row['mean_diff'], row['avg_std'], color = color(row['mut_context']), marker = marker(row['mut_context']), markersize=3, alpha = 0.6)
    
    g.fig.suptitle(f'{len(dataframe.index)} compared bases mean difference against\naverage standard deviation\n{label1} and {label2}', y=0.98)
    g.ax_joint.grid(True, 'both', 'both', alpha = 0.4, linestyle = '--', linewidth = 0.5)

    lims = np.array([
        [-.02, max(dataframe['mean_diff']) + 0.1],
        [-.02, max(dataframe['avg_std']) + 0.1]
    ])

    y1 = np.arange(min(lims[:, 0]), max(lims[:, 1]) + 0.01, 0.01)
    y2 = np.repeat(max(lims[:, 1]), len(y1))

    g.ax_joint.fill_between(y1, lims[1,0], y1, color = '#1b9e77', alpha = 0.15, label = 'significant, TD>=1')
    g.ax_joint.fill_between(y1, y1, y2, color = '#7570b3', alpha = 0.15, label = 'insignificant, TD<1')

    g.ax_joint.set_xlim(tuple(lims[0]))
    g.ax_joint.set_ylim(tuple(lims[1]))

    xlabel = 'Mean Difference'
    ylabel = 'Average Standard Deviation'

    g.ax_joint.set_xlabel(xlabel)
    g.ax_joint.set_ylabel(ylabel)
    legend_mut = mlines.Line2D([], [], color='blue', marker='D', linestyle='None', markersize=10, label='mutation')
    legend_mod = mlines.Line2D([], [], color='#d95f02', marker='o', linestyle='None', markersize=10, label='matching reference')
    sign = mlines.Line2D([], [], color='#1b9e77', marker='s', linestyle='None', markersize=10, label='significant, TD>=1')
    insign = mlines.Line2D([], [], color='#7570b3', marker='s', linestyle='None', markersize=10, label='insignificant, TD<1')
    handles = [legend_mut, legend_mod, sign, insign]
    g.ax_joint.legend(handles = handles, fontsize = 'small', framealpha = 0.3)

    g.plot_marginals(sns.histplot, binwidth = 0.005, kde = True, linewidth = 0)
    g.ax_marg_x.grid(True, 'both', 'both', alpha = 0.4, linestyle = '-', linewidth = 0.5)
    g.ax_marg_y.grid(True, 'both', 'both', alpha = 0.4, linestyle = '-', linewidth = 0.5)
    g.fig.tight_layout()
    g.fig.subplots_adjust(top=0.95)

    plt.savefig(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}_meanDiffStdAvgDist.png'))
    plt.savefig(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}_meanDiffStdAvgDist.pdf'))

    plt.close()

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
    isnan : bool
    '''
    if not sAvg:
        return np.nan, True
    return mDiff/sAvg, 0
    # return np.abs(mDiff - sAvg) / np.sqrt(2)

def kullback_leibler_normal(m0 : float, s0 : float, m1 : float, s1 : float) -> float:
    if not s1 or not s0: # if s0 or s1 are 0 cannot calculate kl divergence
        return np.nan
    return (np.square(s0/s1) + np.square(m1-m0)/np.square(s1) - 1 + np.log(np.square(s1)/np.square(s0))) / 2

def callNanosherlock(working_dir : str, sample_label : str, reference_path : str, fast5_path : str, basecalls_path : str, threads : int, mx : int, mk : int, guppy_bin : str, guppy_model : str, guppy_device : str, fast5_out : bool, calculate_data_density : bool, force_rebuild : bool) -> str:

    red_file_path = os.path.join(working_dir, 'magnipore', sample_label, f'{sample_label}.red')
    if not os.path.exists(red_file_path) or not os.path.exists(reference_path) or force_rebuild:
        command = f'{SUBSCRIPT} {fast5_path} {reference_path} {working_dir} {sample_label} -t {threads} -mx {mx} -mk {mk} -e 1'

        if fast5_out:
            command += ' --fast5_out'
        if force_rebuild:
            command += ' --force_rebuild'
        if calculate_data_density:
            command += ' --calculate_data_density'
        if basecalls_path is not None:
            command += f' --path_to_basecalls {basecalls_path}'
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
            LOGGER.error(f'Error in {SUBSCRIPT} for sample {sample_label} with error code {ret}', error_type=13)
        if TIMEIT:
            LOGGER.printLog(f'{ANSI.YELLOW}TIMED: Calculating distributions of sample {sample_label} took {pd.to_timedelta(end-start)}, {end-start} nanoseconds{ANSI.END}')

    else:
        LOGGER.printLog(f'{sample_label} RED file already exists:\n-\t{red_file_path}')

    return red_file_path

def main():
    
    args = parse()
    
    path_to_fast5_first_sample = args.path_to_fast5_first_sample
    path_to_reference_first_sample = args.path_to_reference_first_sample
    first_sample_label = args.first_sample_label
    
    path_to_fast5_sec_sample = args.path_to_fast5_sec_sample
    path_to_reference_sec_sample = args.path_to_reference_sec_sample
    sec_sample_label = args.sec_sample_label
    
    path_to_first_basecalls = args.path_to_first_basecalls
    path_to_sec_basecalls = args.path_to_sec_basecalls

    working_dir = args.working_dir
    guppy_bin = args.guppy_bin
    guppy_model = args.guppy_model
    guppy_device = args.guppy_device
    
    threads = args.threads
    fast5_out = args.fast5_out
    force_rebuild = args.force_rebuild
    calculate_data_density = args.calculate_data_density

    mx = args.minimap2x
    mk = args.minimap2k

    global TIMEIT 
    TIMEIT = args.timeit

    global LOGGER
    log_file = os.path.join(working_dir, 'log', f'{first_sample_label}_{sec_sample_label}_magnipore.log')
    if not os.path.exists(os.path.join(working_dir, 'log')):
        os.makedirs(os.path.join(working_dir, 'log'))
    LOGGER = Logger(open(log_file, 'w'))
    
    # first sample
    red_first_sample = callNanosherlock(working_dir, first_sample_label, path_to_reference_first_sample, path_to_fast5_first_sample, path_to_first_basecalls, threads, mx, mk, guppy_bin, guppy_model, guppy_device, fast5_out, calculate_data_density, force_rebuild)

    # second sample
    red_sec_sample = callNanosherlock(working_dir, sec_sample_label, path_to_reference_sec_sample, path_to_fast5_sec_sample, path_to_sec_basecalls, threads, mx, mk, guppy_bin, guppy_model, guppy_device, fast5_out, calculate_data_density, force_rebuild)

    # mafft alignment
    alignment_path = mafft(path_to_reference_first_sample, path_to_reference_sec_sample, first_sample_label, sec_sample_label, working_dir, threads, force_rebuild)
    mapping, unaligned, aln_dict = getMapping(alignment_path, working_dir, first_sample_label, sec_sample_label)
    seq_dict = {key:seq.replace('-', '') for key,seq in aln_dict.items()}

    # reading RED files
    red_first_sample = readRedFile(red_first_sample, seq_dict)
    red_sec_sample = readRedFile(red_sec_sample, seq_dict)
    
    if TIMEIT:
        start = perf_counter_ns()
    plotting_data, magnipore_strings = magnipore(mapping, unaligned, seq_dict, aln_dict, red_first_sample, red_sec_sample, first_sample_label, sec_sample_label, working_dir)
    plotStatistics(plotting_data, working_dir, first_sample_label, sec_sample_label)
    writeStockholm(magnipore_strings, alignment_path, first_sample_label, sec_sample_label, working_dir)
    if TIMEIT:
        end = perf_counter_ns()
        LOGGER.printLog(f'{ANSI.YELLOW}TIMED: Evaluating distributions took {pd.to_timedelta(end-start)}, {end - start} nanoseconds{ANSI.END}')

if __name__ == '__main__':
    main()
