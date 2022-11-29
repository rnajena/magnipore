#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import seaborn as sns
import os
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import matplotlib.colors as colors
import numpy as np
from Bio import SeqIO, Seq
from matplotlib import pyplot as plt
from scipy.stats import ks_2samp
import pandas as pd
from Helper import ANSI, MAGNIPORE_COLUMNS
from statistics import NormalDist
from Logger import Logger
import re
from time import perf_counter_ns

LOGGER : Logger = None
FONTSIZE = 18
TIMEIT = False
SUBSCRIPT = 'nanosherlock.py'

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def parse() -> Namespace:

    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        # description='Required tools in environment: guppy, minimap2, nanopolish, h5py, samtools, scipy and mafft\nsee github https://github.com/JannesSP/magnipore'
        )
    
    parser.add_argument("path_to_fast5_first_sample", type = str, help='FAST5 file of first sample')
    parser.add_argument("path_to_reference_first_sample", type = str, help='reference FASTA file of first sample')
    parser.add_argument("first_sample_label", type = str, help='Name of the sample or pipeline run')
    
    parser.add_argument("path_to_fast5_sec_sample", type = str, help='FAST5 file of second sample')
    parser.add_argument("path_to_reference_sec_sample", type = str, help='reference FASTA file of second sample')
    parser.add_argument("sec_sample_label", type = str, help='Name of the sample or pipeline run')
    
    parser.add_argument("working_dir", type = str, help='Path to write all output files')
    
    parser.add_argument("--guppy_bin", type = str, default = None, help='Guppy binary')
    parser.add_argument("--guppy_model", type = str, default = None, help='Guppy model used for basecalling')
    parser.add_argument('--guppy_device', type=str, default='cuda:0', help='Use the GPU to basecall "cuda:0" to use the GPU with ID 0')

    parser.add_argument('--path_to_first_basecalls', type = str, default = None, help = 'FASTQ file to use: <first_sample_label>.fastq\nWill skip basecalling for first sample')
    parser.add_argument('--path_to_sec_basecalls', type = str, default = None, help = 'FASTQ file to use: <sec_sample_label>.fastq\nWill skip basecalling for second sample')
    parser.add_argument('--calculate_data_density', action = 'store_true', default = False, help = 'Will calculate data density after building the models. Will increase runtime!')

    parser.add_argument('-t', "--threads", type=int, default=1, help='Number of threads to use')
    parser.add_argument('-f5', '--fast5_out', action = 'store_true', help='Guppy generates FAST5 output (workspace folder) of Guppy')
    parser.add_argument('-fr', '--force_rebuild', action = 'store_true', help='Run commands regardless if files are already present')
    parser.add_argument('--strict', action = 'store_true', default = False, help = 'Do not write positions with a mutational context into .magnipore files')
    parser.add_argument('-r2', '--range2', action = 'store_true', default = False, help = 'Use range 2 instead of range 3 for the mutational context check')
    parser.add_argument('-mx', '--minimap2x', default = 'splice', choices = ['map-ont', 'splice', 'ava-ont'], help = '-x parameter for minimap2')
    parser.add_argument('-mk', '--minimap2k', default = 14, help = '-k parameter for minimap2')
    parser.add_argument('--timeit', default = False, action = 'store_true', help = 'Measure and print time used by submodules')

    return parser.parse_args()

def mafft(ref_first_sample : str, ref_sec_sample : str, first_sample_label : str, sec_sample_label : str, working_dir : str, threads : int):
    
    # write both references into one file
    alignment_path = os.path.join(working_dir, 'alignment')
    ref_both_samples = os.path.join(alignment_path, f'{first_sample_label}_{sec_sample_label}.fa')
    ref_alignment = os.path.join(alignment_path, f'{first_sample_label}_{sec_sample_label}.aln')
    
    if os.path.exists(ref_alignment):
        return ref_alignment
    
    if not os.path.exists(alignment_path):
        os.makedirs(alignment_path)
    
    # ensure newline between reference files
    command = f'(cat {ref_first_sample}; echo ""; cat {ref_sec_sample}) > {ref_both_samples}'
    
    LOGGER.printLog(f'Writing both references into one file {ref_both_samples}')
    LOGGER.printLog(f'cat command: {ANSI.RED}{command}{ANSI.END}')
    ret = os.system(command)
    
    if ret != 0:
        LOGGER.error('Error in concatenating both reference files')
    
    LOGGER.printLog(f'Building alignment in {ref_alignment}')
    LOGGER.printLog(f'Mafft alignment command: {ANSI.RED}{command}{ANSI.END}')
    
    command = f'mafft --auto --thread {threads} {ref_both_samples} > {ref_alignment}'

    if TIMEIT:
        start = perf_counter_ns()

    ret = os.system(command)

    if TIMEIT:
        end = perf_counter_ns()

    if ret != 0:
        LOGGER.error('Error in building alignment with mafft')

    if TIMEIT:
        LOGGER.printLog(f'TIMED: mafft took {pd.to_timedelta(end-start)}, {end - start} nanoseconds')
    
    return ref_alignment

def getMapping(alignment_path : str, outpath : str, first_label : str, second_label : str):
    '''
    Get base mapping from reference to reference, split indels and substitutions/matches
    '''
    
    fasta = SeqIO.parse(open(alignment_path), 'fasta')
    sequences = {}

    outfile = os.path.join(outpath, 'alignment', f'{first_label}_{second_label}_refdiffs.csv')
    w = open(outfile, 'w')
    w.write(f'type,{first_label},{second_label},base1,base2\n')

    for seq in fasta:

        sequences[seq.id] = str(seq.seq)
        
    LOGGER.printLog(f'Found {len(sequences)} sequences in the alignment.')
    LOGGER.printLog(f'Found an alignment for sequences {list(sequences.keys())}')
    LOGGER.printLog(f'Lengths are {list(map(len, sequences.values()))}')
    
    # {(pos_of_first_sample, base) : (pos_of_second_sample, base)}
    
    try:
        seq1_id, seq2_id = sequences.keys()
        seq1_seq, seq2_seq = sequences.values()
        i1 = i2 = 0
        first2sec_refpos_mapping = {}
        unaligned_positions = {seq : [] for seq in list(sequences.keys())}
        
        for alip, (base1, base2) in enumerate(zip(seq1_seq, seq2_seq)): # (seq1, base1), (seq2, base2) in zip(sequences.items()):
            
            if base1 == '-':
                
                unaligned_positions[seq2_id].append((i2, base2.upper()))
                w.write(f'insert_{second_label},{i1},{i2},{base1},{base2.upper()}\n')
                i2 += 1
                
            elif base2 == '-':
                
                unaligned_positions[seq1_id].append((i1, base1.upper()))
                w.write(f'insert_{first_label},{i1},{i2},{base1.upper()},{base2}\n')
                i1 += 1
            
            else:
            
                first2sec_refpos_mapping[i1] = (i2, alip)
                if base1.lower() != base2.lower():
                    w.write(f'substitution,{i1},{i2},{base1.upper()},{base2.upper()}\n')
                elif base1.lower() == base2.lower() == 'n':
                    w.write(f'N,{i1},{i2},{base1.upper()},{base2.upper()}\n')
                i1 += 1
                i2 += 1
                
        assert i1 == len(seq1_seq.replace('-', '')), 'Mapping iterator for sample 1 does not match sequence length'
        assert i2 == len(seq2_seq.replace('-', '')), 'Mapping iterator for sample 2 does not match sequence length'

        LOGGER.printLog(f'Found {len(first2sec_refpos_mapping)} aligned positions')
        LOGGER.printLog(f'Found {len(unaligned_positions[seq1_id])} unaligned positions for sequence {seq1_id}')
        LOGGER.printLog(f'Found {len(unaligned_positions[seq2_id])} unaligned positions for sequence {seq2_id}')

        return first2sec_refpos_mapping, unaligned_positions, (seq1_id, seq2_id), (seq1_seq, seq2_seq)

    # case: same reference fasta for both samples
    except:
        seq_id = list(sequences.keys())[0]
        seq = list(sequences.values())[0]

        first2sec_refpos_mapping = {i:(i,i) for i in range(len(seq))}
        unaligned_positions = {seq_id : []}

        return first2sec_refpos_mapping, unaligned_positions, (seq_id, seq_id), (seq, seq)

def readRedFile(red_file : str):
    '''
    Reads the red file, stores data in a dictionary
    '''
    
    LOGGER.printLog(f'Reading red file {red_file}')
    
    # sequences are stored as {reference: {pos: {base: ('A'|'C'|'G'|'T'), mean: float, std: float}}}
    red_sequences = {}

    with open(red_file, 'r') as red:
        
        for line in red:
            
            if line.startswith('reference\t'):
                
                continue

            # w.write('reference\tposition\tstrand\tbase\tsignal_mean\tsignal_std\tmotif\tdata_density\texpected_model_density\tn_datapoints\tcontained_datapoints\tn_segments\tcontained_segments\tn_reads\n')
            reference, position, strand, base, mean, std, motif, data_density, expected_model_density, n_datapoints, contained_datapoints, n_segments, contained_segments, n_reads = line.strip().split('\t')
            position = int(position)

            if reference not in red_sequences:
                
                red_sequences[reference] = {}
                
            if position not in red_sequences[reference]:
                
                red_sequences[reference][position] = {
                    '+': {'base' : '', 'motif':'', 'mean' : np.nan, 'std' : np.nan},
                    '-': {'base' : '', 'motif':'', 'mean' : np.nan, 'std' : np.nan}
                    }
                                
            red_sequences[reference][position][strand]['base'] = base
            red_sequences[reference][position][strand]['motif'] = motif
            red_sequences[reference][position][strand]['mean'] = float(mean)
            red_sequences[reference][position][strand]['std'] = float(std)
            red_sequences[reference][position][strand]['data_density'] = float(data_density)
            red_sequences[reference][position][strand]['expected_model_density'] = float(expected_model_density)
            red_sequences[reference][position][strand]['contained_datapoints'] = int(contained_datapoints)
            red_sequences[reference][position][strand]['n_datapoints'] = int(n_datapoints)
            red_sequences[reference][position][strand]['contained_segments'] = int(contained_segments)
            red_sequences[reference][position][strand]['n_segments'] = int(n_segments)
            red_sequences[reference][position][strand]['n_reads'] = int(n_reads)

    return red_sequences

def magnipore(mapping : dict, unaligned : dict, seqs_ids : tuple, alignment_sequences : tuple, alignment_path : str, red_first_sample : dict, red_sec_sample : dict, first_sample_label : str, sec_sample_label : str, working_dir : str, strict : bool, r2 : bool):
    
    suffix = f'{"_strict" if strict else ""}{"_r2" if r2 else ""}'
    working_dir = os.path.join(working_dir, 'magnipore', f'{first_sample_label}_{sec_sample_label}{suffix}')

    if not os.path.exists(working_dir):
        os.makedirs(working_dir)

    all_file = os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}{suffix}.all')
    magnipore_file = os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}{suffix}.magnipore')
    indels_file = os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}{suffix}.indels')

    indels = open(indels_file, 'w')
    indels.write(f'type\tstrand\tref\tpos\tbase\n')
    
    magnipore = open(magnipore_file, 'w')
    magnipore.write('\t'.join(MAGNIPORE_COLUMNS) + '\n')
    all = open(all_file, 'w')
    all.write('\t'.join(MAGNIPORE_COLUMNS) + '\n')

    # red: sequences are stored as {reference: {pos: {base: ('A'|'C'|'G'|'T'), mean: float, std: float}}}
    num_indels, sign_pos, nans, alignmentGapCorrection = 0, 0, 0, 0
    # replace every nucleotide character with a dot
    alignment_sequences = [list(re.sub(r"[^-]", ".", alignment_sequences[0])), list(re.sub(r"[^-]", ".", alignment_sequences[1]))]

    # TODO add some quality value
    plotting_data = pd.DataFrame(columns=['mean_diff', 'first_std', 'sec_std', 'avg_std', 'mut_context', 'td_score', 'kl_divergence'])
    plotting_data = plotting_data.astype(
        {
            'mean_diff': 'float',
            'first_std': 'float',
            'sec_std': 'float',
            'avg_std': 'float',
            'mut_context': 'str',
            'td_score':'float',
            'kl_divergence':'float'
        })
    num_motif_diff = 0

    LOGGER.printLog(f'Start comparing sequences position wise ...')
    mappingLength = len(mapping)

    # compare distributions of aligned positions
    for sample_idx, (pos_first_sample, (pos_sec_sample, alip)) in enumerate(mapping.items()):
        
        if (sample_idx + 1) % 1000 == 0:
            print(f'\t{sample_idx + 1}/{mappingLength}', end='\r')

        # no information for comparison, rare case in start and end of reference/contig (nanopolish artefact, limit of 5-mer models)
        if pos_first_sample not in red_first_sample[seqs_ids[0]] or pos_sec_sample not in red_sec_sample[seqs_ids[1]]:
            
            nans += 1
            continue
        
        dist_first_sample = red_first_sample[seqs_ids[0]][pos_first_sample]
        dist_sec_sample = red_sec_sample[seqs_ids[1]][pos_sec_sample]
        
        for strand in ['+', '-']:
            
            # check if both aligned positions have a distributions -> stdev for both are non-0
            if dist_first_sample[strand]['std'] and dist_sec_sample[strand]['std']:

                m0 = dist_first_sample[strand]['mean']
                s0 = dist_first_sample[strand]['std']
                m1 = dist_sec_sample[strand]['mean']
                s1 = dist_sec_sample[strand]['std']

                firstDist = NormalDist(m0, s0)
                secDist = NormalDist(m1, s1)
                
                mDiff = abs(m0 - m1)
                sAvg = (s0 + s1)/2
                kl_divergence = kullback_leibler_normal(m0, s0, m1, s1)
                td = td_score(mDiff, sAvg)

                first_motif = dist_first_sample[strand]['motif']
                sec_motif = dist_sec_sample[strand]['motif']

                # only look at context of range 2
                if r2:
                    first_motif = first_motif[len(first_motif)//2 - 2 : len(first_motif)//2 + 3]
                    sec_motif = sec_motif[len(sec_motif)//2 - 2 : len(sec_motif)//2 + 3]
    
                    mut_context = first_motif != sec_motif

                # range of 3
                else:

                    if len(first_motif) == len(sec_motif):

                        mut_context = first_motif != sec_motif
                    
                    # motifs have different lengths, rare case in start and end of reference/contig
                    else:

                        mut_context = first_motif[len(first_motif)//2 - 2 : len(first_motif)//2 + 3] != sec_motif[len(sec_motif)//2 - 2 : len(sec_motif)//2 + 3]

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

                all.write(f'{strand}\t{td}\t{kl_divergence}\t{firstDist.overlap(secDist)}\t{"mut" if mut_context else "mod"}\t{seqs_ids[0]}\t{pos_first_sample}\t{dist_first_sample[strand]["base"]}\t{dist_first_sample[strand]["motif"]}\t{m0}\t{s0}\t{dist_first_sample[strand]["n_datapoints"]}\t{dist_first_sample[strand]["contained_datapoints"]}\t{dist_first_sample[strand]["n_segments"]}\t{dist_first_sample[strand]["contained_segments"]}\t{dist_first_sample[strand]["n_reads"]}\t{seqs_ids[1]}\t{pos_sec_sample}\t{dist_sec_sample[strand]["base"]}\t{dist_sec_sample[strand]["motif"]}\t{m1}\t{s1}\t{dist_sec_sample[strand]["n_datapoints"]}\t{dist_sec_sample[strand]["contained_datapoints"]}\t{dist_sec_sample[strand]["n_segments"]}\t{dist_sec_sample[strand]["contained_segments"]}\t{dist_sec_sample[strand]["n_reads"]}\n')

                # distance between both means is greater than the average std of both distributions
                if mDiff > sAvg:

                    # 7mer motif or if at border 5mer motif
                    # both 7mers or both 5mers or one 7mer one 5mer
                    if mut_context:
                        
                        num_motif_diff += 1

                        # skip substitutions
                        if strict:
                            continue
                    
                    # X == interesting position
                    alignment_sequences[0][alip] = 'X'
                    alignment_sequences[1][alip] = 'X'

                    sign_pos += 1
                    
                    magnipore.write(f'{strand}\t{td}\t{kl_divergence}\t{firstDist.overlap(secDist)}\t{"mut" if mut_context else "mod"}\t{seqs_ids[0]}\t{pos_first_sample}\t{dist_first_sample[strand]["base"]}\t{dist_first_sample[strand]["motif"]}\t{m0}\t{s0}\t{dist_first_sample[strand]["n_datapoints"]}\t{dist_first_sample[strand]["contained_datapoints"]}\t{dist_first_sample[strand]["n_segments"]}\t{dist_first_sample[strand]["contained_segments"]}\t{dist_first_sample[strand]["n_reads"]}\t{seqs_ids[1]}\t{pos_sec_sample}\t{dist_sec_sample[strand]["base"]}\t{dist_sec_sample[strand]["motif"]}\t{m1}\t{s1}\t{dist_sec_sample[strand]["n_datapoints"]}\t{dist_sec_sample[strand]["contained_datapoints"]}\t{dist_sec_sample[strand]["n_segments"]}\t{dist_sec_sample[strand]["contained_segments"]}\t{dist_sec_sample[strand]["n_reads"]}\n')

            # one of the aligned and compared position is not present in the .red file
            else:
                nans += 1


    all.close()
    magnipore.close()
    print(f'\t{sample_idx + 1}/{mappingLength}')
    LOGGER.printLog(f'Number of indels: {num_indels}\n'\
               f'Number of significant positions: {sign_pos}\n'\
               f'Number of significant positions with reference differences: {num_motif_diff}\n'\
               f'Number of nans {nans}, at least one aligned position without information (no signals)\nCan be high if one strand has no information!\n'\
               f'Wrote {magnipore_file}')

    with open(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}{suffix}.txt'), 'w') as w:
        w.write(f'Number of indels: {num_indels}\n'\
                f'Number of significant positions: {sign_pos}\n'\
                f'Number of significant positions with reference differences: {num_motif_diff}\n'\
                f'Number of nans {nans}, at least one aligned position without information (no signals)\nCan be high if one strand has no information!\n')

    LOGGER.printLog('Start writing files ...')

    for seq in unaligned:
        
        for position, base in unaligned[seq]:
            
            indels.write(f'insert\t+\t{seq}\t{position}\t{base}\n')
            num_indels += 1
            
    indels.close()
    LOGGER.printLog('Done with indels')
    
    plot_dir = os.path.join(working_dir, 'plots')
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    ### Mean Dist vs Std Avg plot
    LOGGER.printLog('Plotting Mean vs Stdev')
    plotMeanDiffStdAvg(plotting_data, plot_dir, first_sample_label, sec_sample_label, suffix)

    ### std vs std plot with mean_diff colored scatter
    # LOGGER.printLog(f'Plotting Stdev {first_sample_label} vs Stdev {sec_sample_label}')
    # plotStdStdMean(plotting_data, plot_dir, first_sample_label, sec_sample_label, suffix)

    ### mean vs mean plot with avg_std colored scatter
    # LOGGER.printLog(f'Plotting mean {first_sample_label} vs mean {sec_sample_label}')
    # plotMeanMeanStd(plotting_data, plot_dir, first_sample_label, sec_sample_label, suffix)

    ### plot scores
    LOGGER.printLog(f'Plotting TD score and KL divergence')
    plotScores(plotting_data, plot_dir, first_sample_label, sec_sample_label, suffix)

    ### plot mean stdev distributions
    # LOGGER.printLog(f'Plotting Mean and Stdev for {first_sample_label}')
    # plotMeanStdDist(plotting_data, 'first_mean', 'first_std', plot_dir, first_sample_label, suffix)

    # LOGGER.printLog(f'Plotting Mean and Stdev for {sec_sample_label}')
    # plotMeanStdDist(plotting_data, 'sec_mean', 'sec_std', plot_dir, sec_sample_label, suffix)

    LOGGER.printLog(f'Writing stockholm with magnipore markers')
    fasta = []
    ids = []
    # remove duplicate sequences (happens when both samples have the same reference fasta)
    for record in SeqIO.parse(open(alignment_path), 'fasta'):
        if record.id not in ids:
            fasta.append(record)
            ids.append(record.id)
    fasta += [SeqIO.SeqRecord(
                Seq.Seq(''.join(seq)),
                id ='magnipore_marked_' + str((first_sample_label, sec_sample_label)[i]),
                name ='magnipore_' + str((first_sample_label, sec_sample_label)[i]),
                description='X=significant signal change, N=not significant') for i, seq in enumerate(alignment_sequences)]   
    stk = open(os.path.join(working_dir, first_sample_label + '_' + sec_sample_label + suffix + '_marked.stk'), 'w')
    SeqIO.write(fasta, stk, 'stockholm')
    
    LOGGER.printLog('Done with magnipore')

def plotMeanStdDist(dataframe : pd.DataFrame, mean_column : str, stdev_column : str, working_dir : str, sample_label : str, suffix : str) -> None:
    
    plt.figure(figsize = (12,8), dpi=300)
    plt.title(f'{sample_label} mean and stdev distribution')
    sns.kdeplot(data = dataframe, x = mean_column, y = stdev_column, fill = True)
    plt.grid(True,  'both', 'both')
    plt.xlabel('\u03BC Mean')
    plt.ylabel('\u03C3 Standard deviation')
    plt.tight_layout()
    plt.savefig(os.path.join(working_dir, f'{sample_label}{suffix}_msDist.png'))
    plt.savefig(os.path.join(working_dir, f'{sample_label}{suffix}_msDist.pdf'))
    plt.close()

    plt.figure(figsize = (12,8), dpi=300)
    plt.title(f'{sample_label} mean and stdev distribution')
    sns.kdeplot(data = dataframe, x = mean_column, y = stdev_column, fill = True, log_scale=(False, True))
    plt.grid(True,  'both', 'both')
    plt.xlabel('\u03BC Mean')
    plt.ylabel('\u03C3 Standard deviation')
    plt.tight_layout()
    plt.savefig(os.path.join(working_dir, f'{sample_label}{suffix}_msDist_logscale.png'))
    plt.savefig(os.path.join(working_dir, f'{sample_label}{suffix}_msDist_logscale.pdf'))
    plt.close()

def plotScores(dataframe : pd.DataFrame, working_dir : str, first_sample_label : str, sec_sample_label : str, suffix : str) -> None:
    
    dataframe[', '.join(['Context', 'Significance'])] =  pd.Series(dataframe.reindex(['mut_context', 'significant'], axis='columns').astype('str').values.tolist()).str.join(', ')

    plt.figure(figsize = (12,8), dpi=300)
    plt.title(f'TD score for all positions\n{first_sample_label} vs. {sec_sample_label}')
    sns.histplot(data = dataframe, x = 'td_score', hue=dataframe['Context, Significance'], log_scale = (True, True), multiple="stack")
    plt.grid(True,  'both', 'both')
    plt.tight_layout()
    plt.savefig(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}{suffix}_td_score.png'))
    plt.savefig(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}{suffix}_td_score.pdf'))
    plt.close()

    plt.figure(figsize = (12,8), dpi=300)
    plt.title(f'Kullback-Leibler divergence for all positions\n{first_sample_label} vs. {sec_sample_label}')
    sns.histplot(data = dataframe, x = 'kl_divergence', hue=dataframe['Context, Significance'], log_scale = (True, True), multiple="stack")
    plt.grid(True,  'both', 'both')
    plt.tight_layout()
    plt.savefig(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}{suffix}_kl_div.png'))
    plt.savefig(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}{suffix}_kl_div.pdf'))
    plt.close()

def plotStdStdMean(dataframe : pd.DataFrame, working_dir : str, first_sample_label : str, sec_sample_label : str, suffix : str) -> None:
    
    plt.figure(figsize = (12,12), dpi=300)
    plt.rcParams.update({'font.size': FONTSIZE})
    g = sns.JointGrid(x = 'first_std', y = 'sec_std', data = dataframe, marginal_ticks=True, height = 10)
    hue_norm = plt.Normalize(min(dataframe['mean_diff']), max(dataframe['mean_diff']))
    cmap = truncate_colormap(plt.get_cmap('Reds'), minval = 0.4)
    g.plot_joint(sns.scatterplot, hue = dataframe['mean_diff'], hue_norm = hue_norm, palette = cmap, legend = False, s = 12, alpha = 0.6)
    g.fig.suptitle(f'{len(dataframe.index)} compared bases\nstdev of {first_sample_label} and {sec_sample_label}\ncolored with mean difference')
    g.ax_joint.grid(True, 'both', 'both', alpha = 0.4, linestyle = '-', linewidth = 0.5)

    plt.xlabel(f'{first_sample_label} stdev')
    plt.ylabel(f'{sec_sample_label} stdev')

    g.plot_marginals(sns.histplot, binwidth = 0.005, kde = True, linewidth = 0)
    g.ax_marg_x.grid(True, 'both', 'both', alpha = 0.4, linestyle = '-', linewidth = 0.5)
    g.ax_marg_y.grid(True, 'both', 'both', alpha = 0.4, linestyle = '-', linewidth = 0.5)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=hue_norm)
    plt.colorbar(sm, label = 'mean diff', aspect = 50)
    g.fig.subplots_adjust(top=0.95)
    g.fig.tight_layout()
    plt.savefig(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}{suffix}_stds.png'))
    plt.savefig(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}{suffix}_stds.pdf'))
    g.ax_joint.set_xscale('log')
    g.ax_joint.set_yscale('log')
    plt.savefig(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}{suffix}_stds_logscale.png'))
    plt.savefig(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}{suffix}_stds_logscale.pdf'))

    plt.close()
    
def plotMeanDiffStdAvg(dataframe : pd.DataFrame, working_dir : str, first_sample_label : str, sec_sample_label : str, suffix : str) -> None:
    
    ### Mean Dist vs Std Avg plot
    plt.figure(figsize = (12,12), dpi=300)
    plt.rcParams.update({
        'font.size': FONTSIZE,
        "text.usetex": True,
        "font.family": "sans-serif",
        "font.sans-serif": "Helvetica",
        })
    label1 = first_sample_label.replace("_", " ")
    label2 = sec_sample_label.replace("_", " ")

    g = sns.JointGrid(x = 'mean_diff', y = 'avg_std', data = dataframe, hue = 'mut_context', marginal_ticks=True, palette=['blue', 'orange'], hue_order=['mutation', 'matching reference'], height = 10)
    g.plot_joint(sns.scatterplot, s = 12) #, markers = ['.', '^']
    g.fig.suptitle(f'{len(dataframe.index)} compared bases\nmean difference against average standard deviation\n{label1} and {label2}')
    g.ax_joint.grid(True, 'both', 'both', alpha = 0.4, linestyle = '--', linewidth = 0.5)

    lims = np.array([
        [-.02, max(dataframe['mean_diff']) + 0.1],
        [-.02, max(dataframe['avg_std']) + 0.1]
    ])

    y1 = np.arange(min(lims[:, 0]), max(lims[:, 1]) + 0.01, 0.01)
    y2 = np.repeat(max(lims[:, 1]), len(y1))

    g.ax_joint.fill_between(y1, lims[1,0], y1, color = 'green', alpha = 0.08, label = 'significant signals, TD(.,.)$\ge1$')
    g.ax_joint.fill_between(y1, y1, y2, color = 'red', alpha = 0.08, label = 'insignificant signals, TD(.,.)$<1$')

    g.ax_joint.set_xlim(tuple(lims[0]))
    g.ax_joint.set_ylim(tuple(lims[1]))

    # \u03BC is mu
    # \u03C3 is sigma
    # alternativ \mbox{} \footnotesize

    xlabel = '$\mu_{\mbox{\small ' + label1 + '}} - \mu_{\mbox{\small ' + label2 + '}}$ mean difference'
    ylabel = '$\\frac{\sigma_{\mbox{\small ' + label1 + '}}\mbox{ }+\mbox{ }\sigma_{\mbox{\small ' + label2 + '}}}{2}$ average standard deviation'

    # print(xlabel)
    # print(ylabel)

    g.ax_joint.set_xlabel(xlabel)
    g.ax_joint.set_ylabel(ylabel)
    g.ax_joint.legend(fontsize = 'small', framealpha = 0.3)

    g.plot_marginals(sns.histplot, binwidth = 0.005, kde = True, linewidth = 0)
    g.ax_marg_x.grid(True, 'both', 'both', alpha = 0.4, linestyle = '-', linewidth = 0.5)
    g.ax_marg_y.grid(True, 'both', 'both', alpha = 0.4, linestyle = '-', linewidth = 0.5)
    g.fig.tight_layout()
    g.fig.subplots_adjust(top=0.95)

    plt.savefig(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}{suffix}_meanDiffStdAvgDist.png'))
    plt.savefig(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}{suffix}_meanDiffStdAvgDist.pdf'))

    g.ax_joint.set_ylim(bottom = min(dataframe['avg_std']))
    g.ax_joint.set_yscale('log')

    plt.savefig(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}{suffix}_meanDiffStdAvgDist_logscale.png'))
    plt.savefig(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}{suffix}_meanDiffStdAvgDist_logscale.pdf'))

    plt.close()

def plotMeanMeanStd(dataframe : pd.DataFrame, working_dir : str, first_sample_label : str, sec_sample_label : str, suffix : str) -> None:
    
    plt.figure(figsize = (12,12), dpi=300)
    plt.rcParams.update({'font.size': FONTSIZE})
    g = sns.JointGrid(x = 'first_mean', y = 'sec_mean', data = dataframe, marginal_ticks=True, height = 10)

    log_avg_std = np.log(dataframe['avg_std'])

    hue_norm = plt.Normalize(min(log_avg_std), max(log_avg_std))
    cmap = truncate_colormap(plt.get_cmap('Reds'), minval = 0.4)
    g.plot_joint(sns.scatterplot, hue = log_avg_std, hue_norm = hue_norm, palette = cmap, legend = False, s = 12, alpha = 0.6)
    g.fig.suptitle(f'{len(dataframe.index)} compared bases:\n{first_sample_label} mean against {sec_sample_label} mean')
    g.ax_joint.grid(True, 'both', 'both', alpha = 0.4, linestyle = '-', linewidth = 0.5)

    plt.xlim((min(dataframe['first_mean']), max(dataframe['first_mean']) + 0.05))
    plt.ylim((min(dataframe['sec_mean']), max(dataframe['sec_mean']) + 0.05))

    plt.xlabel(f'{first_sample_label} mean')
    plt.ylabel(f'{sec_sample_label} mean')

    g.plot_marginals(sns.histplot, binwidth = 0.005, kde = True, linewidth = 0)
    g.ax_marg_x.grid(True, 'both', 'both', alpha = 0.4, linestyle = '-', linewidth = 0.5)
    g.ax_marg_y.grid(True, 'both', 'both', alpha = 0.4, linestyle = '-', linewidth = 0.5)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=hue_norm)
    plt.colorbar(sm, label = 'log avg_std', aspect = 50)
    g.fig.subplots_adjust(top=0.95)
    g.fig.tight_layout()
    g.fig.subplots_adjust(top=0.95)

    plt.tight_layout()

    plt.savefig(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}{suffix}_mmDist.png'))
    plt.savefig(os.path.join(working_dir, f'{first_sample_label}_{sec_sample_label}{suffix}_mmDist.pdf'))

    plt.close()

def ks_test(dist1 : tuple, dist2 : tuple) -> tuple:
    
    data1 = np.random.normal(dist1[0], dist1[1], 100)
    data2 = np.random.normal(dist2[0], dist2[1], 100)
    
    return ks_2samp(data1, data2)

def td_score(mDiff, sAvg) -> float:
    '''
        Calculates the td-score mDiff/sAvg
    '''
    return mDiff/sAvg
    # return np.abs(mDiff - sAvg) / np.sqrt(2)

def kullback_leibler_normal(m0 : float, s0 : float, m1 : float, s1 : float) -> float:
    if not s1 or not s0: # if s0 or s1 are 0 cannot calculate kl divergence
        return np.nan
    return (np.square(s0/s1) + np.square(m1-m0)/np.square(s1) - 1 + np.log(np.square(s1)/np.square(s0))) / 2

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
    strict = args.strict
    r2 = args.range2
    calculate_data_density = args.calculate_data_density

    mx = args.minimap2x
    mk = args.minimap2k

    global TIMEIT 
    TIMEIT = args.timeit

    basedir = os.path.dirname(__file__)

    global LOGGER
    log_file = os.path.join(working_dir, 'log', f'{first_sample_label}_{sec_sample_label}_magnipore.log')
    if not os.path.exists(os.path.join(working_dir, 'log')):
        os.makedirs(os.path.join(working_dir, 'log'))
    LOGGER = Logger(open(log_file, 'w'))
    
    # first sample
    red_first_sample = os.path.join(working_dir, 'magnipore', first_sample_label, f'{first_sample_label}.red')

    if not os.path.exists(red_first_sample) or not os.path.exists(path_to_reference_first_sample) or force_rebuild:
    
        command_first_sample = f'python3 {os.path.join(basedir, SUBSCRIPT)} {path_to_fast5_first_sample} {path_to_reference_first_sample} {working_dir} {first_sample_label} -t {threads} -mx {mx} -mk {mk}'

        if fast5_out:
            command_first_sample += ' --fast5_out'

        if force_rebuild:
            command_first_sample += ' --force_rebuild'

        if calculate_data_density:
            command_first_sample += ' --calculate_data_density'

        if path_to_first_basecalls is not None:
            command_first_sample += f' --path_to_basecalls {path_to_first_basecalls}'
        else:
            assert guppy_bin is not None and guppy_model is not None, 'Need at least the guppy binary path and model or path to basecalls'
            command_first_sample += f' --guppy_bin {guppy_bin} --guppy_model {guppy_model} --guppy_device {guppy_device}'
        
        if TIMEIT:
            command_first_sample = command_first_sample + ' --timeit'
            start = perf_counter_ns()

        LOGGER.printLog(f'Pipeline command: {ANSI.RED}{command_first_sample}{ANSI.END}')
        ret = os.system(command_first_sample)

        if TIMEIT:
            end = perf_counter_ns()

        if ret != 0:
            LOGGER.error(f'Error in {SUBSCRIPT} for sample {first_sample_label}')

        if TIMEIT:
            LOGGER.printLog(f'TIMED: Calculating distributions of sample {first_sample_label} took {pd.to_timedelta(end-start)}, {end-start} nanoseconds')

    else:
        LOGGER.printLog(f'{first_sample_label} red file already exists:\n-\t{red_first_sample}')
        
    # second sample
    red_sec_sample = os.path.join(working_dir, 'magnipore', sec_sample_label, f'{sec_sample_label}.red')
    
    if not os.path.exists(red_sec_sample) or not os.path.exists(path_to_reference_sec_sample) or force_rebuild:
        command_sec_sample = f'python3 {os.path.join(basedir, SUBSCRIPT)} {path_to_fast5_sec_sample} {path_to_reference_sec_sample} {working_dir} {sec_sample_label} -t {threads} -mx {mx} -mk {mk}'

        if fast5_out:
            command_sec_sample += ' --fast5_out'

        if force_rebuild:
            command_sec_sample += ' --force_rebuild'

        if calculate_data_density:
            command_sec_sample += ' --calculate_data_density'

        if path_to_sec_basecalls is not None:
           command_sec_sample += f' --path_to_basecalls {path_to_sec_basecalls}'
        else:
            assert guppy_bin is not None and guppy_model is not None, 'Need at least the guppy binary path and model or path to basecalls'
            command_sec_sample += f' --guppy_bin {guppy_bin} --guppy_model {guppy_model} --guppy_device {guppy_device}'

        if TIMEIT:
            command_sec_sample = command_sec_sample + ' --timeit'
            start = perf_counter_ns()

        LOGGER.printLog(f'Pipeline command: {ANSI.RED}{command_sec_sample}{ANSI.END}')
        ret = os.system(command_sec_sample)

        if TIMEIT:
            end = perf_counter_ns()
        
        if ret != 0:
            LOGGER.error(f'Error in {SUBSCRIPT} for sample {sec_sample_label}')
    
        if TIMEIT:
            LOGGER.printLog(f'TIMED: Calculating distributions of sample {sec_sample_label} took {pd.to_timedelta(end-start)}, {end-start} nanoseconds')

    else:
        LOGGER.printLog(f'{sec_sample_label} red file already exists:\n-\t{red_sec_sample}')

    # mafft alignment
    alignment_path = mafft(path_to_reference_first_sample, path_to_reference_sec_sample, first_sample_label, sec_sample_label, working_dir, threads)
    mapping, unaligned, seq_ids, aligned_sequences = getMapping(alignment_path, working_dir, first_sample_label, sec_sample_label)
    
    red_first_sample = readRedFile(red_first_sample)
    red_sec_sample = readRedFile(red_sec_sample)
    
    if TIMEIT:
        start = perf_counter_ns()

    magnipore(mapping, unaligned, seq_ids, aligned_sequences, alignment_path, red_first_sample, red_sec_sample, first_sample_label, sec_sample_label, working_dir, strict, r2)

    if TIMEIT:
        end = perf_counter_ns()
        LOGGER.printLog(f'TIMED: Evaluating distributions took {pd.to_timedelta(end-start)}, {end - start} nanoseconds')

if __name__ == '__main__':
    main()
