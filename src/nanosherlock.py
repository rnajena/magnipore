# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import datetime
import OnlineMeanVar as omv
import os
import subprocess
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from pathlib import Path
import h5py
import numpy as np
from scipy import stats
from Bio import SeqIO
from constants import ANSI, IUPAC, COMPLEMENT
import pandas as pd
import seaborn as sns
from matplotlib.pyplot import figure
from matplotlib import pyplot as plt
from Logger import Logger

LOGGER : Logger = None
        
def readFast5(fast5_path : str, sequencing_summary : str):
    '''
    Returns a dictionary mapping the readid to the corresponding fast5 file containing the read.
    '''

    LOGGER.printLog(f'Loading readids from fast5 files in {fast5_path} ... ')

    readid2fast5 = {}

    if sequencing_summary is None:

        fast5list = Path(fast5_path).rglob('*.fast5')

        for fidx, fast5 in enumerate(fast5list):

            LOGGER.printLog(f'Indexing fast5 file {fidx + 1}\r', newline_after=False)
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
    LOGGER.printLog(f'Start loading ids from {nanoSum_path} ...')
    with open(nanoSum_path, 'r') as summary:
        next(summary)

        for line in summary:
            r_index, r_ID = line.strip().split()[:2]
            read_index2ID[int(r_index)] = r_ID

    LOGGER.printLog(f'Done: found {len(read_index2ID)} readids')
    return read_index2ID

def parse() -> Namespace:

    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description='Required tools in environment: see github https://github.com/JannesSP/magnipore'
        ) 
    
    parser.add_argument('path_to_fast5', type = str, help='FAST5 directory')
    parser.add_argument('path_to_reference', type = str, help='reference FASTA file')
    parser.add_argument('working_dir', type = str, help='Path to write all output files')
    parser.add_argument('guppy_bin', type = str, help='Guppy binary')
    parser.add_argument('guppy_model', type = str, help='Guppy model for basecalling')
    parser.add_argument('sample_label', type = str, help='Name of the sample or pipeline run')

    parser.add_argument('--path_to_basecalls', default = None, type = str, help = 'Path to existing basecalls and sequencing summary file. Basecalls must be in one single file with the name <sample_label>.fastq')
    parser.add_argument('--medaka_model', type = str, default = 'r941_min_sup_g507', help = 'Medaka model, see more in Medaka\'s documentation')
    parser.add_argument('--calculate_data_density', action = 'store_true', default = False, help = 'Will calculate data density after building the models. Will increase runtime!')

    parser.add_argument('-t', '--threads', type=int, help='Number of threads to use')
    parser.add_argument('-f5', '--fast5_out', action = 'store_true', help='Guppy generates FAST5 output (workspace folder)')
    parser.add_argument('-fr', '--force_rebuild', action = 'store_true', help='Run commands regardless if files are already present')
    parser.add_argument('-c', '--consensus_tool', default = 'skip', choices=['ococo', 'medaka', 'skip'], help = 'Choose which consensus caller to use to refine the reference.')
    parser.add_argument('-mx', '--minimap2x', default = 'splice', choices = ['map-ont', 'splice', 'ava-ont'], help = '-x parameter for minimap2')
    parser.add_argument('-mk', '--minimap2k', default = 14, type = int, help = '-k parameter for minimap2')

    return parser.parse_args()

def guppy_basecalling(guppy_bin : str, guppy_model : str,  path_to_fast5 : str, working_dir : str, sample_label : str, fast5_out : bool, force_rebuild : bool, sequencing_summary : str) -> str:

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
        
        LOGGER.printLog(f'Creating sub-directory {basecalls_path}')
        os.makedirs(basecalls_path)
    
    LOGGER.printLog(f'Run guppy basecalling on {path_to_fast5}')
    command = f'{guppy_bin} -i {path_to_fast5} -s {basecalls_path} --disable_pings --device cuda:0 --disable_qscore_filtering --calib_detect -c {guppy_model}'

    if fast5_out:
        
        command += ' --fast5_out'
        
    LOGGER.printLog(f'guppy command: {ANSI.RED}{command}{ANSI.END}')

    ret = os.system(command)
    
    if ret != 0:
        
        LOGGER.error('Error in guppy basecalling!')
        
    os.system(f'cat {os.path.join(basecalls_path, "*.fastq")} > {os.path.join(basecalls_path, "*.tmpfastq")}')
    os.system(f'rm {os.path.join(basecalls_path, "*.fastq")}')
    os.system(f'mv {os.path.join(basecalls_path, "*.tmpfastq")} {basecalls}')
        
    if fast5_out:

        basecalled_fast5_path = os.path.join(working_dir, "basecalled_fast5", sample_label)
        LOGGER.printLog(f'Move basecalled fast5 files to {basecalled_fast5_path}')
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
        LOGGER.printLog(f'Creating sub-directory {minimap_path}')
        os.makedirs(minimap_path)
        
    LOGGER.printLog(f'Run minimiap2 on basecalls: {path_to_basecalls} and reference: {path_to_reference}')
    command = f'minimap2 -a -x {mx} -k{mk} -t {threads} {path_to_reference} {path_to_basecalls} | samtools view -hbF4 | samtools sort > {bam_path}'
    LOGGER.printLog(f'minimap2 command: {ANSI.RED}{command}{ANSI.END}')
    
    ret = os.system(command)
    
    if ret != 0:
        
        LOGGER.error('Error in minimap2 run!')
        
    return bam_path, force_rebuild

def ococo(alignment_bam : str, path_to_reference : str, working_dir : str, sample_label : str, force_rebuild : bool) -> str:
    
    ococo_path = os.path.join(working_dir, 'ococo', sample_label)
    ococo_out_files = os.path.join(ococo_path, sample_label)

    if os.path.exists(ococo_out_files + '_consensus.fa') and os.path.exists(ococo_out_files + '.vcf') and not force_rebuild:

        LOGGER.printLog(f'ococo files already exists:\n-\t{ococo_out_files + "_consensus.fa"}\n-\t{ococo_out_files + ".vcf"}')

        return ococo_out_files, force_rebuild

    else:

        force_rebuild = True
    
    if not os.path.exists(ococo_path):
        LOGGER.printLog(f'Creating sub-directory {ococo_path}')
        os.makedirs(ococo_path)

    LOGGER.printLog(f'Run ococo on bam file: {alignment_bam} and reference: {path_to_reference}')
    command = f'ococo -i {alignment_bam} -f {path_to_reference} --vcf-cons {ococo_out_files}.vcf --fasta-cons {ococo_out_files}_consensus.fa --stats-out {ococo_out_files}.stats --pileup {ococo_out_files}.pileup'
    LOGGER.printLog(f'ococo command: {ANSI.RED}{command}{ANSI.END}')
    
    ret = os.system(command)
    
    if ret != 0:
        LOGGER.error('Error in ococo run!')

    return ococo_out_files, force_rebuild

def medaka(path_to_basecalls : str, path_to_reference : str, working_dir : str, sample_label : str, force_rebuild : bool, threads : int, model : str = 'r941_min_sup_g507') -> str:
    
    medaka_path = os.path.join(working_dir, 'medaka', sample_label)
    medaka_out_file = os.path.join(medaka_path, 'consensus.fasta')
    sample_out_file = os.path.join(medaka_path, sample_label + '.fasta')
    medaka_hdf_file = os.path.join(medaka_path, 'consensus_probs.hdf')

    if os.path.exists(sample_out_file) and not force_rebuild:

        LOGGER.printLog(f'medaka file already exists:\n-\t{sample_out_file}')

        return sample_out_file, force_rebuild

    else:

        force_rebuild = True

    if not os.path.exists(medaka_path):
        LOGGER.printLog(f'Creating sub-directory {medaka_path}')
        os.makedirs(medaka_path)

    LOGGER.printLog(f'Run medaka consensus on basecalls {path_to_basecalls} and reference {path_to_reference} with model {model}')
    command = f'medaka_consensus -t {threads} -i {path_to_basecalls} -d {path_to_reference} -o {medaka_path} -m {model}'
    LOGGER.printLog(f'medaka_consensus command: {ANSI.RED}{command}{ANSI.END}')
    subprocess.run(f'. ~/anaconda3/etc/profile.d/conda.sh && conda activate medaka && {command} && conda deactivate', shell=True) # , executable='/bin/bash')
    
    LOGGER.printLog(f'Run medaka variant on reference {path_to_reference} and the consensus {medaka_out_file}')
    command = f'medaka variant {path_to_reference} {medaka_hdf_file} {os.path.join(medaka_path, sample_label + ".vcf")}'
    LOGGER.printLog(f'medaka variant command: {ANSI.RED}{command}{ANSI.END}')
    subprocess.run(f'. ~/anaconda3/etc/profile.d/conda.sh && conda activate medaka && {command} && conda deactivate', shell=True)
    
    with open(medaka_out_file, 'r') as r:        
        with open(sample_out_file, 'w') as w:
            for line in r:
                if line.startswith('>'):
                    w.write(f'>{sample_label}\n')
                else:
                    w.write(line)
    
    return sample_out_file, force_rebuild

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

        LOGGER.printLog(f'Creating sub-directory {nanopolish_path}')
        os.makedirs(nanopolish_path)

    if not os.path.exists(alignment_bam + '.bai') or force_rebuild:

        LOGGER.printLog(f'Indexing {alignment_bam} with samtools')
        command = f'samtools index {alignment_bam}'
        LOGGER.printLog(f'samtools index command: {ANSI.RED}{command}{ANSI.END}')
    
        ret = os.system(command)
        
        if ret != 0:
            LOGGER.error('Error in samtools indexing for nanopolish!')
    
    else:

        LOGGER.printLog(f"Using existing {os.path.exists(alignment_bam + '.bai')}")

    if not os.path.exists(path_to_basecalls + '.index') or force_rebuild:

        LOGGER.printLog(f'Run nanopolish indexing')
        command = f'nanopolish index -d {path_to_fast5} -s {path_to_sequencing_summary} {path_to_basecalls}'
        LOGGER.printLog(f'nanopolish indexing command: {ANSI.RED}{command}{ANSI.END}')
        
        ret = os.system(command)
        
        if ret != 0:
            LOGGER.error('Error in nanopolish indexing!')

    else:

        LOGGER.printLog(f"Using existing {os.path.exists(path_to_basecalls + '.index')}")

    
    LOGGER.printLog(f'Run nanopolish eventalign')
    command = f'nanopolish eventalign --reads {path_to_basecalls} --bam {alignment_bam} --genome {path_to_reference} --summary={summary_csv} --scale-events --signal-index -t {threads} --progress > {result_csv}'
    LOGGER.printLog(f'Nanopolish eventalign command: {ANSI.RED}{command}{ANSI.END}')
    
    ret = os.system(command)
    
    if ret != 0:
        LOGGER.error('Error in nanopolish eventalign!')
        
    return summary_csv, result_csv, force_rebuild

def read_line(line):
    
    line_dic = {}
    line = line.strip().split()

    line_dic['contig'] = line[0]
    line_dic['position'] = int(line[1]) + 2
    line_dic['ref_kmer'] = line[2]
    line_dic['read_index'] = int(line[3])
    line_dic['strand'] = line[4]
    line_dic['event_index'] = int(line[5])
    line_dic['event_mean'] = float(line[6])
    line_dic['event_std'] = float(line[7])
    line_dic['event_len'] = float(line[8])
    line_dic['model_kmer'] = line[9]
    line_dic['model_mean'] = float(line[10])
    line_dic['model_std'] = float(line[11])
    line_dic['standardized_level'] = float(line[12])
    line_dic['start_idx'] = int(line[13])
    line_dic['end_idx'] = int(line[14])

    return line_dic

def getZNormSignal(fast5_read) -> np.ndarray:
        
    dacs = fast5_read['Raw/Signal'][:]
    r_offset = fast5_read['channel_id'].attrs['offset']
    r_range = fast5_read['channel_id'].attrs['range']
    r_digit = fast5_read['channel_id'].attrs['digitisation']
    
    # https://github.com/nanoporetech/taiyaki/blob/master/docs/FILE_FORMATS.md
    # Dacs is the pA signal transformed to integers
    # using offset, range and digit to transform it back
    read_pA = (dacs + r_offset) * (r_range / r_digit)
    
    # https://github.com/nanoporetech/taiyaki/blob/9672c66d654ad8101c4faf481646493082e090c8/taiyaki/maths.py#L8
    # shift and scale are calculated as described by nanoporetech, see the links above
    r_shift = np.median(read_pA)
    r_scale = 1.4826 * np.median(abs(read_pA - r_shift))
    
    norm_signal = (read_pA - r_shift) / r_scale
    
    return norm_signal

def aggregate_events(nanopolish_result_csv : str, nanopolish_summary_csv: str, path_to_fast5 : str, path_to_reference : str, working_dir : str, sample_label : str, force_rebuild : bool, sequencing_summary : str, calculate_data_density : bool) -> str:
    
    if not os.path.exists(os.path.join(working_dir, 'magnipore', sample_label)):
        os.makedirs(os.path.join(working_dir, 'magnipore', sample_label))

    red = os.path.join(working_dir, 'magnipore', sample_label, f'{sample_label}.red')

    if os.path.exists(red) and not force_rebuild:
        return red, force_rebuild

    else:
        force_rebuild = True

    readid2fast5 = readFast5(path_to_fast5, sequencing_summary)
    nano2readid = readNanoSum(nanopolish_summary_csv)
    sequences = createSeqDict(path_to_reference)
    
    LOGGER.printLog(f'Found sequences: {list(sequences.keys())}')
    LOGGER.printLog(f'Start getting event signal distributions per reference position and write red file for {sample_label}')
    buildModels(sequences, nano2readid, readid2fast5, nanopolish_result_csv, calculate_data_density)

    LOGGER.printLog('Writing output files')
    writeOutput(red, sequences, working_dir, sample_label, calculate_data_density)
        
    return red, force_rebuild

def createSeqDict(path_to_reference : str) -> dict:

    sequences = {}
    fasta = SeqIO.parse(open(path_to_reference), 'fasta')

    for seq in fasta:        
        sequences[seq.id] = {}
        
        for pos, base in enumerate(seq.seq):
            sequences[seq.id][pos] = {
                
                '+': {
                    'base' : base,
                    'motif' : np.nan,
                    'mean' : np.nan,
                    'std' : np.nan,
                    'omv' : omv.LocShift(100),
                    'data_density' : 0,
                    'contained_datapoints' : 0,
                    'contained_segments' : 0,
                    'n_datapoints' : 0,
                    'n_segments' : 0,
                    'n_reads': 0,
                    'nanopolish_model_mismatch': {b : 0 for b in 'ACGT'}
                    },
                
                '-': {
                    'base': COMPLEMENT.get(base, 'N'),
                    'motif' : np.nan,
                    'mean' : np.nan,
                    'std' : np.nan,
                    'omv' : omv.LocShift(100),
                    'data_density' : 0,
                    'contained_datapoints' : 0,
                    'contained_segments' : 0,
                    'n_datapoints' : 0,
                    'n_segments' : 0,
                    'n_reads': 0,
                    'nanopolish_model_mismatch': {b : 0 for b in 'ACGT'}
                    }
                }

    return sequences

def buildModels(sequences : dict, nano2readid : dict, readid2fast5 : dict, nanopolish_result_csv : str, calculate_data_density : bool):

    max_lines = 0

    for loop in ['buildModel', 'checkData']:
        LOGGER.printLog(f'Starting {loop}')
            
        with open(nanopolish_result_csv, 'r') as nano_result:
            nano_result.readline()

            opened_readid = ''
            opened_fast5 = ''
            last_position = -1

            for lidx, line in enumerate(nano_result):

                if lidx + 1 > max_lines:
                    max_lines = lidx + 1

                if (lidx + 1) % 100000 == 0:
                    LOGGER.printLog(f'Line {lidx + 1}/{max_lines}, {loop}\r', newline_after=False)

                event = read_line(line)
                strand = '+' if event['ref_kmer'] == event['model_kmer'] else '-'
                red = sequences[event['contig']][event['position']][strand]

                if event['model_kmer'] == 'NNNNN':
                    continue
                    
                readid = nano2readid[event['read_index']]
                
                if readid != opened_readid:
                    
                    if readid2fast5[readid] != opened_fast5:
                    
                        if opened_fast5 != '':
                            fast5.close()

                        fast5 = h5py.File(readid2fast5[readid], 'r')
                        opened_fast5 = readid2fast5[readid]
                        last_position = -1
                        
                    fast5_read = fast5['read_' + readid]
                    norm_signal = getZNormSignal(fast5_read)
                    opened_readid = readid
                segment = norm_signal[event['start_idx']:event['end_idx']]
                
                if loop == 'buildModel':

                    if event['ref_kmer'][2] not in IUPAC[sequences[event['contig']][event['position']][strand]['base']]:
                        red['nanopolish_model_mismatch'][event['ref_kmer'][2]] += 1

                    if red['n_segments'] <= 0:
                    
                        motif = sequences[event["contig"]][event["position"] - 2][strand]["base"] + \
                                sequences[event["contig"]][event["position"] - 1][strand]["base"] + \
                                sequences[event["contig"]][event["position"]][strand]["base"] + \
                                sequences[event["contig"]][event["position"] + 1][strand]["base"] + \
                                sequences[event["contig"]][event["position"] + 2][strand]["base"]

                        try:
                            motif = sequences[event["contig"]][event["position"] - 3][strand]["base"] + motif + sequences[event["contig"]][event["position"] + 3][strand]["base"]
                        except:
                            pass

                        red['motif'] = motif

                    red['omv'].append(segment)
                    red['n_datapoints'] += event['end_idx'] - event['start_idx']
                    red['n_segments'] += 1

                    if event['position'] != last_position:
                        if last_position != -1:
                            sequences[event['contig']][last_position][strand]['n_reads'] += 1
                        last_position = event['position']
                
                elif loop == 'checkData':

                    red['mean'], red['std'] = red['omv'].meanStdev()

                    r = 3 * red['std'] # ~99% density of normal distribution
                    l = red['mean'] - r # lower bound
                    u = red['mean'] + r # upper bound

                    contained = ((l <= segment) & (segment <= u)).sum()
                    if contained == len(segment):
                        red['contained_segments'] += 1
                    red['contained_datapoints'] += contained
                    
                    if calculate_data_density:
                        red['data_density'] += np.mean(stats.norm.pdf(
                            segment,
                            loc = red['mean'],
                            scale = red['std']))

        LOGGER.printLog(f'Line {lidx + 1}/{max_lines}, {loop}\r', newline_after=False)
        print()

    if calculate_data_density: 
        for contig in sequences:
            for positions in sequences[contig]:
                for strand in sequences[contig][positions]:
                    if sequences[contig][positions][strand]['n_segments'] != 0:
                        sequences[contig][positions][strand]['data_density'] = sequences[contig][positions][strand]['data_density'] / sequences[contig][positions][strand]['n_segments']

def writeOutput(red : str, sequences : dict, working_dir : str, sample_label : str, calculate_data_density : bool):

    nans = 0
    plotting_data = pd.DataFrame(columns=['contig', 'position', 'strand', 'Fraction of contained datapoints', 'Fraction of contained segments', 'Density difference', 'n_datapoints', 'n_segments', 'n_reads'])
    plotting_data = plotting_data.astype(
        {
            'contig' : 'str',
            'position' : 'int',
            'strand' : 'str',
            'Fraction of contained datapoints': 'float',
            'Fraction of contained segments' : 'float',
            'Density difference': 'float',
            'n_datapoints' : 'int',
            'n_segments' : 'int',
            'n_reads' : 'int'
        })
    
    with open(red, 'w') as w:

        w.write('reference\tposition\tstrand\tbase\tsignal_mean\tsignal_std\tmotif\tdata_density\texpected_model_density\tn_datapoints\tcontained_datapoints\tn_segments\tcontained_segments\tn_reads\tnanomA\tnanomC\tnanomG\tnanomT\n')

        for sequence in sequences:
            for position in sequences[sequence]:
                for strand in sequences[sequence][position]:

                    event = sequences[sequence][position][strand]

                    if not np.isnan(event['std']):
                        expected_model_density = 1 / (2 * np.sqrt(np.pi) * event['std'])
                    else:
                        nans += 1
                        expected_model_density = np.NAN

                    w.write(f'{sequence}\t{position}\t{strand}\t')
                    w.write(f'{event["base"]}\t{event["mean"]}\t{event["std"]}\t{event["motif"]}\t{event["data_density"]}\t{expected_model_density}\t{event["n_datapoints"]}\t{event["contained_datapoints"]}\t{event["n_segments"]}\t{event["contained_segments"]}\t{event["n_reads"]}')

                    for base in event['nanopolish_model_mismatch']:
                        w.write(f"\t{event['nanopolish_model_mismatch'][base]}")
            
                    w.write('\n')

                    new_entry = pd.DataFrame({
                            'contig': [sequence],
                            'position': [position],
                            'strand' : [strand],
                            'Fraction of contained datapoints' : [event["contained_datapoints"] / event["n_datapoints"] if event["n_datapoints"] != 0 else 0],
                            'Fraction of contained segments' : [event["contained_segments"] / event["n_segments"] if event["n_segments"] != 0 else 0],
                            'Density difference': [expected_model_density - event["data_density"] if not np.isnan(expected_model_density) else np.NAN],
                            'n_datapoints' : [event["n_datapoints"]],
                            'n_segments' : [event["n_segments"]],
                            'n_reads' : [event["n_reads"]]

                    })
                    plotting_data = pd.concat([plotting_data, new_entry], ignore_index=True)

    LOGGER.printLog(f'Positions without information: {nans}')

    plotStatistics(plotting_data, working_dir, sample_label, calculate_data_density)

def plotStatistics(dataFrame : pd.DataFrame, working_dir : str, sample_label : str, calculate_data_density : bool) -> None:
    LOGGER.printLog('Plotting data ...')
    working_dir = os.path.join(working_dir, 'magnipore', sample_label, 'plots')
    
    if not os.path.exists(working_dir):
        os.mkdir(working_dir)

    figure(figsize = (12,8), dpi=2000)
    g = sns.histplot(data = dataFrame, x = 'Fraction of contained datapoints', kde = True, hue = 'contig', stat = 'density')
    plt.title('Distribution of contained datapoints fractions of all position')
    plt.setp(g.get_legend().get_texts(), fontsize='6') # for legend text
    plt.setp(g.get_legend().get_title(), fontsize='6') # for legend title
    plt.tight_layout()
    plt.savefig(os.path.join(working_dir, f'{sample_label}_datapointsHist.png'))
    plt.savefig(os.path.join(working_dir, f'{sample_label}_datapointsHist.pdf'))
    plt.close()

    figure(figsize = (12,8), dpi=2000)
    g = sns.histplot(data = dataFrame, x = 'Fraction of contained segments', kde = True, hue = 'contig', stat = 'density')
    plt.title('Distribution of contained segments fractions of all position')
    plt.setp(g.get_legend().get_texts(), fontsize='6') # for legend text
    plt.setp(g.get_legend().get_title(), fontsize='6') # for legend title
    plt.tight_layout()
    plt.savefig(os.path.join(working_dir, f'{sample_label}_segmentsHist.png'))
    plt.savefig(os.path.join(working_dir, f'{sample_label}_segmentsHist.pdf'))
    plt.close()

    figure(figsize = (12,8), dpi=2000)
    g = sns.histplot(data = dataFrame, x = 'Fraction of contained datapoints', y = 'Fraction of contained segments', cbar = True)
    plt.title('Contained datapoints vs segments fractions of all position')
    plt.tight_layout()
    plt.savefig(os.path.join(working_dir, f'{sample_label}_datasegHist.png'))
    plt.savefig(os.path.join(working_dir, f'{sample_label}_datasegHist.pdf'))
    plt.close()

    figure(figsize = (12,8), dpi=2000)
    g = sns.barplot(data = dataFrame, x = 'position', y = 'n_reads', hue = 'contig', linewidth = 0)
    g.set_xticks(range(0, max(dataFrame['position']), max(dataFrame['position'])//10))
    plt.title(f'Read coverage of segmented signals for sample {sample_label}')
    plt.setp(g.get_legend().get_texts(), fontsize='6') # for legend text
    plt.setp(g.get_legend().get_title(), fontsize='6') # for legend title
    g.set_yscale("log")
    plt.tight_layout()
    plt.savefig(os.path.join(working_dir, f'{sample_label}_readCoverage.png'))
    plt.savefig(os.path.join(working_dir, f'{sample_label}_readCoverage.pdf'))
    plt.close()

    figure(figsize = (12,8), dpi=2000)
    g = sns.barplot(data = dataFrame, x = 'position', y = 'n_segments', hue = 'contig', linewidth = 0)
    g.set_xticks(range(0, max(dataFrame['position']), max(dataFrame['position'])//10))
    plt.title(f'Segment coverage of segmented signals for sample {sample_label}')
    plt.setp(g.get_legend().get_texts(), fontsize='6') # for legend text
    plt.setp(g.get_legend().get_title(), fontsize='6') # for legend title
    g.set_yscale("log")
    plt.tight_layout()
    plt.savefig(os.path.join(working_dir, f'{sample_label}_segmentCoverage.png'))
    plt.savefig(os.path.join(working_dir, f'{sample_label}_segmentCoverage.pdf'))
    plt.close()

    figure(figsize = (12,8), dpi=2000)
    g = sns.barplot(data = dataFrame, x = 'position', y = 'n_datapoints', hue = 'contig', linewidth = 0)
    g.set_xticks(range(0, max(dataFrame['position']), max(dataFrame['position'])//10))
    plt.title(f'Signal coverage of segmented signals for sample {sample_label}')
    plt.setp(g.get_legend().get_texts(), fontsize='6') # for legend text
    plt.setp(g.get_legend().get_title(), fontsize='6') # for legend title
    g.set_yscale("log")
    plt.tight_layout()
    plt.savefig(os.path.join(working_dir, f'{sample_label}_signalCoverage.png'))
    plt.savefig(os.path.join(working_dir, f'{sample_label}_signalCoverage.pdf'))
    plt.close()

    figure(figsize = (12,8), dpi=2000)
    g = sns.histplot(data = dataFrame, x = 'Density difference', kde = True, hue = 'contig', stat = 'density')
    plt.title('Distribution of model density vs data density difference of all position.' if calculate_data_density else 'Model density per position.')
    plt.setp(g.get_legend().get_texts(), fontsize='6') # for legend text
    plt.setp(g.get_legend().get_title(), fontsize='6') # for legend title
    plt.tight_layout()
    plt.savefig(os.path.join(working_dir, f'{sample_label}_densityDiff.png' if calculate_data_density else f'{sample_label}_modeldensity.png'))
    plt.savefig(os.path.join(working_dir, f'{sample_label}_densityDiff.pdf' if calculate_data_density else f'{sample_label}_modeldensity.pdf'))
    plt.close()

def main():
    
    args = parse()
    
    working_dir = args.working_dir
    path_to_fast5 = args.path_to_fast5
    path_to_reference = args.path_to_reference
    guppy_bin = args.guppy_bin
    guppy_model = args.guppy_model
    sample_label = args.sample_label
    fast5_out = args.fast5_out
    threads = args.threads
    force_rebuild = args.force_rebuild
    con_caller = args.consensus_tool
    path_to_basecalls = args.path_to_basecalls
    mx = args.minimap2x
    mk = args.minimap2k
    medaka_model = args.medaka_model
    calculate_data_density = args.calculate_data_density
        
    log_file = os.path.join(working_dir, 'log', f'{sample_label}_magnipore_{datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}.log')
    if not os.path.exists(os.path.join(working_dir, 'log')):
        os.makedirs(os.path.join(working_dir, 'log'))
    global LOGGER
    LOGGER = Logger(open(log_file, 'w'))

    if not os.path.exists(working_dir):
        
        os.mkdir(working_dir)
        LOGGER.printLog(f'Creating working directory {working_dir}')
            
    LOGGER.printLog(f'Starting magnipore pipeline')
    LOGGER.printLog(f'Writing log to {log_file}')

    sequencing_summary = None

    if path_to_basecalls is not None:

        sequencing_summary = os.path.join(path_to_basecalls, 'sequencing_summary.txt')

        path_to_basecalls = os.path.join(path_to_basecalls, sample_label + '.fastq')

        if not os.path.exists(path_to_basecalls):
            
            LOGGER.error(f'{path_to_basecalls} NOT FOUND!')
            
        if not os.path.exists(sequencing_summary):

            LOGGER.error(f'{sequencing_summary} NOT FOUND!')

    else:

        path_to_basecalls, sequencing_summary, force_rebuild = guppy_basecalling(guppy_bin, guppy_model, path_to_fast5, working_dir, sample_label, fast5_out, force_rebuild, sequencing_summary)
    
    if con_caller == 'ococo':

        alignment_bam, force_rebuild = minimap(path_to_reference, path_to_basecalls, working_dir, sample_label, threads, force_rebuild, mx, mk)
        ococo_files, force_rebuild = ococo(alignment_bam, path_to_reference, working_dir, sample_label, force_rebuild)
        consensus = ococo_files + '_consensus.fa'

    elif con_caller == 'medaka':

        consensus, force_rebuild = medaka(path_to_basecalls, path_to_reference, working_dir, sample_label, force_rebuild, threads, medaka_model)

    elif con_caller == 'skip':

        consensus = path_to_reference

    alignment_bam, force_rebuild = minimap(consensus, path_to_basecalls, working_dir, sample_label, threads, force_rebuild, mx, mk)    
    nanopolish_summary_csv, nanopolish_result_csv, force_rebuild = nanopolish(path_to_fast5, sequencing_summary, path_to_basecalls, consensus, alignment_bam, working_dir, sample_label, threads, force_rebuild)
    red_file, force_rebuild = aggregate_events(nanopolish_result_csv, nanopolish_summary_csv, path_to_fast5, consensus, working_dir, sample_label, force_rebuild, sequencing_summary, calculate_data_density)
    LOGGER.printLog(f'Done with {sample_label}')
    LOGGER.printLog(f'Aggregated reference event distributions are stored in {red_file}')

    return red_file

if __name__ == '__main__':
    main()
