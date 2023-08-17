import os

import numpy as np
import pandas as pd
from Bio import SeqIO

from magnipore import nanosherlock
from magnipore.Helper import REDENCODER
from magnipore.magnipore import align, getMapping, readRedFile
from magnipore import magnipore
from magnipore.nanosherlock import (aggregate_events, mapping,
                                    signalSegmentation)

magnipore.initLogger(None)
nanosherlock.initLogger(None)
this_file_dir = os.path.dirname(__file__)
lab1 = 'sample1'
ref1 = os.path.join(this_file_dir, 'sample1.fa')
red1 = os.path.join(this_file_dir, 'sample1.red')
lab2 = 'sample2'
ref2 = os.path.join(this_file_dir, 'sample2.fa')
red2 = os.path.join(this_file_dir, 'sample2.red')
red1DF = pd.read_csv(red1, sep='\t') # assuming pandas is tested and working correctly
red2DF = pd.read_csv(red2, sep='\t') # assuming pandas is tested and working correctly
seq_dict = {**SeqIO.to_dict(SeqIO.parse(ref1, format='fasta')), **SeqIO.to_dict(SeqIO.parse(ref2, format='fasta'))}
# print(seq_dict)
s=['+','-']
processes = 8
# ===================================
alignment_path = align(ref1, ref2, lab1, lab2, this_file_dir, 1, False)
mapping_dict, unaligned_dict, aln_dict = getMapping(alignment_path, this_file_dir, lab1, lab2)
red1_dict = readRedFile(red1, seq_dict)
red2_dict = readRedFile(red2, seq_dict)
sign_positions = [(2,2), (4,4), (10,5), (16,11)]
tdscores = [1.0801634950514278, 1.998636430295431, 2.201571772538712, 1.9101836858234442, 9.5]
bayesian_ps = [0.58914066,0.31764056,0.27098986,0.33953125,0.00000178]
classification = ['mut','mut','mut','mod', 'mod']
strands = ['+','+','+','+','-']
motifs_sample1 = ['ATCAA', 'CAATT', 'TTGGA', 'CAACA', 'TGTTG']
motifs_sample2 = ['ACCAA', 'CAAGG', 'AAGGA', 'CAACA', 'TGTTG']
bases_sample1 = 'CAGAT'
bases_sample2 = 'CAGAT'
# ===================================
targets = [38, 366]
target_mean = [-0.2045173406607994, 2.584417871497208]
target_std = [0.6050989516056662, 0.42849883560181956]
def test_ReadREDFile():

    for sequence in red1_dict:
        for i,pos in enumerate(red1_dict[sequence]):
            for j,strand in enumerate(pos):
                assert strand[REDENCODER['mean']] == red1DF['signal_mean'].iloc[[(2*i)+j]][red1DF['strand'] == s[j]].item()
                assert strand[REDENCODER['std']] == red1DF['signal_std'].iloc[[(2*i)+j]][red1DF['strand'] == s[j]].item()
                assert strand[REDENCODER['data_density']] == red1DF['data_density'].iloc[[(2*i)+j]][red1DF['strand'] == s[j]].item()
                if np.isnan(strand[REDENCODER['expected_model_density']]):
                    assert np.isnan(red1DF['expected_model_density'].iloc[[(2*i)+j]][red1DF['strand'] == s[j]].item())
                else:
                    assert strand[REDENCODER['expected_model_density']], red1DF['expected_model_density'].iloc[[(2*i)+j]][red1DF['strand'] == s[j]].item()
                assert strand[REDENCODER['n_datapoints']] == red1DF['n_datapoints'].iloc[[(2*i)+j]][red1DF['strand'] == s[j]].item()
                assert strand[REDENCODER['contained_datapoints']] == red1DF['contained_datapoints'].iloc[[(2*i)+j]][red1DF['strand'] == s[j]].item()
                assert strand[REDENCODER['n_segments']] == red1DF['n_segments'].iloc[[(2*i)+j]][red1DF['strand'] == s[j]].item()
                assert strand[REDENCODER['contained_segments']] == red1DF['contained_segments'].iloc[[(2*i)+j]][red1DF['strand'] == s[j]].item()
                assert strand[REDENCODER['n_reads']] == red1DF['n_reads'].iloc[[(2*i)+j]][red1DF['strand'] == s[j]].item()

def test_mapping():
    alignment = {0: (0, 0), 1: (1, 1), 2: (2, 2), 3: (3, 3), 4: (4, 4), 10: (5, 10), 11: (6, 11), 12: (7, 12), 13: (8, 13), 14: (9, 14), 15: (10, 15), 16: (11, 16), 17: (12, 17), 18: (13, 18), 19: (14, 19), 20: (
15, 20), 21: (16, 21)}
    assert mapping_dict == alignment
    unaligned = {'sample1': [(5, 'T'), (6, 'T'), (7, 'T'), (8, 'T'), (9, 'T'), (22, 'T')], 'sample2': []}
    assert unaligned_dict == unaligned
    assert aln_dict == {lab1:'atcaatttttggaccaacacttt', lab2:'accaa-----ggaccaacactt-'}

def test_magnipore():
    seq_dict = {key:seq.replace('-', '') for key,seq in aln_dict.items()}
    magnipore.magnipore(mapping_dict, unaligned_dict, seq_dict, aln_dict, alignment_path, red1_dict, red2_dict, lab1, lab2, this_file_dir, 'r9', processes)
    magn = pd.read_csv(os.path.join(this_file_dir, 'magnipore', 'sample1_sample2', 'sample1_sample2.magnipore'), sep='\t')
    for i, row in magn[magn['strand'] == 0].iterrows():
        assert np.isclose(row['td_score'], tdscores[i])
        assert np.isclose(row['bayesian_p'], bayesian_ps[i])
        assert (row[['ref_1', 'ref_2']].values == (lab1, lab2)).all()
        assert (row[['pos_1', 'pos_2']].values == sign_positions[i]).all()
        assert row['signal_type'] == classification[i]
        assert (row[['motif_1', 'motif_2']].values == (motifs_sample1[i], motifs_sample2[i])).all()
        assert (row[['base_1', 'base_2']].values == (bases_sample1[i], bases_sample2[i])).all()
        assert (row[['signal_mean_1', 'signal_std_1']].values == red1DF[['signal_mean', 'signal_std']].iloc[[row['pos_1']*2]].values.squeeze()).all()
        assert (row[['signal_mean_2', 'signal_std_2']].values == red2DF[['signal_mean', 'signal_std']].iloc[[row['pos_2']*2]].values.squeeze()).all()
        assert (row[['n_datapoints_1','contained_datapoints_1','n_segments_1','contained_segments_1','n_reads_1']].values == red1DF[['n_datapoints','contained_datapoints','n_segments','contained_segments','n_reads']].iloc[[row['pos_1']*2]].values.squeeze()).all()
        assert (row[['n_datapoints_2','contained_datapoints_2','n_segments_2','contained_segments_2','n_reads_2']].values == red2DF[['n_datapoints','contained_datapoints','n_segments','contained_segments','n_reads']].iloc[[row['pos_2']*2]].values.squeeze()).all()
    row = magn[magn['strand'] == '-']
    assert row['td_score'].item() == tdscores[-1]
    assert row['bayesian_p'].item() == bayesian_ps[-1]
    for item in row.items():
        print(item)
    assert row['motif_1'].item() == row['motif_2'].item() and row['motif_2'].item() == motifs_sample1[-1]
    assert row['base_1'].item() == row['base_2'].item() and row['base_2'].item() == bases_sample1[-1]

def test_stockholm():
    stk_file = os.path.join(this_file_dir, 'magnipore/sample1_sample2/sample1_sample2_marked.stk')
    with open(stk_file, 'r') as r:
        for line in r:
            line = line.strip().split(' ')
            if line[0] == 'magnipore_marked_sample1':
                assert line[-1].count('X') == 4
                assert line[-1].count('A') == 0
                assert line[-1].count('C') == 0
                assert line[-1].count('G') == 0
                assert line[-1].count('T') == 0
                assert line[-1].count('U') == 0
                assert line[-1].count('a') == 0
                assert line[-1].count('c') == 0
                assert line[-1].count('g') == 0
                assert line[-1].count('t') == 0
                assert line[-1].count('u') == 0
                assert line[-1].count('.') == 19
                assert line[-1].count('-') == 0
            if line[0] == 'magnipore_marked_sample2':
                assert line[-1].count('X') == 4
                assert line[-1].count('A') == 0
                assert line[-1].count('C') == 0
                assert line[-1].count('G') == 0
                assert line[-1].count('T') == 0
                assert line[-1].count('U') == 0
                assert line[-1].count('a') == 0
                assert line[-1].count('c') == 0
                assert line[-1].count('g') == 0
                assert line[-1].count('t') == 0
                assert line[-1].count('u') == 0
                assert line[-1].count('.') == 13
                assert line[-1].count('-') == 6
    
def test_event_aggregation_fast5():
    raw_data_path = os.path.join(this_file_dir, 'raw_data')
    basecalls = os.path.join(this_file_dir, 'basecalls', 'test.fastq')
    seq_sum = os.path.join(this_file_dir, 'basecalls', 'sequencing_summary.txt')
    map_file = os.path.join(this_file_dir, 'mapping', 'test', 'test.bam')
    reference = os.path.join(this_file_dir, 'reference', 'test.fasta')
    segmentation = os.path.join(this_file_dir, 'segmentation', 'test', 'eventalign_result.csv')
    assert os.path.exists(basecalls)
    mapping(reference, basecalls, this_file_dir, 'test', 1, True, 'splice', 14)
    assert os.path.exists(map_file)
    summary, segmentation, _ = signalSegmentation(raw_data_path, '.fast5', basecalls, reference, map_file, this_file_dir, 'test_fast5', 1, True, True, False, None)
    assert os.path.exists(segmentation)
    assert os.path.exists(summary)
    red_file = aggregate_events(segmentation, summary, raw_data_path, '.fast5', reference, this_file_dir, 'test_fast5', True, seq_sum, False, False)
    assert os.path.exists(red_file)
    # 38: 36 0 54767	54773   0008609d-0d3e-46e5-9b69-25f7ab4b194e
    # 38: 36 2 56368	56421   00425ffc-17d7-4ba0-87ae-9c01215661ca
    # 38: 36 3 15561	15586   00118376-02d0-40a7-88db-5b450adebe13
    # 366: 364 0    41422:41445   0008609d-0d3e-46e5-9b69-25f7ab4b194e
    # 366: 364 2    42862:42909   00425ffc-17d7-4ba0-87ae-9c01215661ca
    # 366: 364 5    25219:25226   003deea8-84e6-4161-9659-12a9fee2cfd4
    checkRed(red_file)

def test_event_aggregation_slow5():
    raw_data_path = os.path.join(this_file_dir, 'raw_data', 'test.slow5')
    basecalls = os.path.join(this_file_dir, 'basecalls', 'test.fastq')
    seq_sum = os.path.join(this_file_dir, 'basecalls', 'sequencing_summary.txt')
    map_file = os.path.join(this_file_dir, 'mapping', 'test', 'test.bam')
    reference = os.path.join(this_file_dir, 'reference', 'test.fasta')
    segmentation = os.path.join(this_file_dir, 'segmentation', 'test', 'eventalign_result.csv')
    assert os.path.exists(basecalls)
    mapping(reference, basecalls, this_file_dir, 'test', 1, True, 'splice', 14)
    assert os.path.exists(map_file)
    summary, segmentation, _ = signalSegmentation(raw_data_path, '.slow5', basecalls, reference, map_file, this_file_dir, 'test_slow5', 1, True, True, False, None)
    assert os.path.exists(segmentation)
    assert os.path.exists(summary)
    red_file = aggregate_events(segmentation, summary, raw_data_path, '.slow5', reference, this_file_dir, 'test_slow5', True, seq_sum, False, False)
    assert os.path.exists(red_file)
    checkRed(red_file)

def test_event_aggregation_blow5():
    raw_data_path = os.path.join(this_file_dir, 'raw_data', 'test.blow5')
    basecalls = os.path.join(this_file_dir, 'basecalls', 'test.fastq')
    seq_sum = os.path.join(this_file_dir, 'basecalls', 'sequencing_summary.txt')
    map_file = os.path.join(this_file_dir, 'mapping', 'test', 'test.bam')
    reference = os.path.join(this_file_dir, 'reference', 'test.fasta')
    segmentation = os.path.join(this_file_dir, 'segmentation', 'test', 'eventalign_result.csv')
    assert os.path.exists(basecalls)
    mapping(reference, basecalls, this_file_dir, 'test', 1, True, 'splice', 14)
    assert os.path.exists(map_file)
    summary, segmentation, _ = signalSegmentation(raw_data_path, '.slow5', basecalls, reference, map_file, this_file_dir, 'test_blow5', 1, True, True, False, None)
    assert os.path.exists(segmentation)
    assert os.path.exists(summary)
    red_file = aggregate_events(segmentation, summary, raw_data_path, '.slow5', reference, this_file_dir, 'test_blow5', True, seq_sum, False, False)
    assert os.path.exists(red_file)
    checkRed(red_file)

def checkRed(red_file : str) -> None:
    relative_tolerance = 0.05 # allow for error as the online update of meanVar is heuristic and inaccurate for low number of values in case of stdev
    with open(red_file, 'r') as red:
        next(red) # skip header
        for line in red:
            line = line.strip().split('\t')
            if int(line[1]) == 38 and line[2] == '+':
                assert np.isclose(float(line[3]), target_mean[0], rtol = relative_tolerance)
                assert np.isclose(float(line[4]), target_std[0], rtol = relative_tolerance)
            elif int(line[1]) == 366 and line[2] == '+':
                assert np.isclose(float(line[3]), target_mean[1], rtol = relative_tolerance)
                assert np.isclose(float(line[4]), target_std[1], rtol = relative_tolerance)