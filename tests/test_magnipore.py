from magnipore.magnipore import readRedFile, initLogger, mafft, getMapping, magnipore
from magnipore.Helper import REDENCODER
from Bio import SeqIO
import os
import numpy as np
import pandas as pd

initLogger(None)
this_file_dir = os.path.dirname(__file__)
lab1 = 'sample1'
ref1 = os.path.join(this_file_dir, 'sample1.fa')
red1 = os.path.join(this_file_dir, 'sample1.red')
lab2 = 'sample2'
ref2 = os.path.join(this_file_dir, 'sample2.fa')
red2 = os.path.join(this_file_dir, 'sample2.red')
red1DF = pd.read_csv(red1, sep='\t') # assuming pandas is tested and working correctly
red2DF = pd.read_csv(red2, sep='\t') # assuming pandas is tested and working correctly
seq_dict = SeqIO.to_dict(SeqIO.parse(ref1, format='fasta')) | SeqIO.to_dict(SeqIO.parse(ref2, format='fasta'))
# print(seq_dict)
s=['+','-']
# ===================================
alignment_path = mafft(ref1, ref2, lab1, lab2, this_file_dir, 1, False)
mapping_dict, unaligned_dict, aln_dict = getMapping(alignment_path, this_file_dir, lab1, lab2)
red1_dict = readRedFile(red1, seq_dict)
red2_dict = readRedFile(red2, seq_dict)
sign_positions = [(2,2), (4,4), (10,5), (16,11)]
tdscores = [1.0801634950514278, 1.998636430295431, 2.201571772538712, 1.9101836858234442, 9.5]
bayesian_ps = [0.58914066,0.31764056,0.27098986,0.33953125,0.00000178]
classification = ['mut','mut','mut','mod', 'mod']
strands = ['+','+','+','+','-']
motifs_sample1 = ['ATCAA', 'TCAATTT', 'TTTGGAC', 'CCAACAC', 'GTGTTGG']
motifs_sample2 = ['ACCAA', 'CCAAGGA', 'CAAGGAC', 'CCAACAC', 'GTGTTGG']
bases_sample1 = 'CAGAT'
bases_sample2 = 'CAGAT'

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
    plotting_data, magnipore_strings = magnipore(mapping_dict, unaligned_dict, seq_dict, aln_dict, red1_dict, red2_dict, lab1, lab2, this_file_dir)
    magn = pd.read_csv(os.path.join(this_file_dir, 'magnipore', 'sample1_sample2', 'sample1_sample2.magnipore'), sep='\t')
    for i, row in magn[magn['strand'] == '+'].iterrows():
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
    assert row['motif_1'].item() == row['motif_2'].item() == motifs_sample1[-1]
    assert row['base_1'].item() == row['base_2'].item() == bases_sample1[-1]
    for l in magnipore_strings:
        assert l.count('X') == 4
        assert l.count('A') == 0
        assert l.count('C') == 0
        assert l.count('G') == 0
        assert l.count('T') == 0
        assert l.count('U') == 0
        assert l.count('a') == 0
        assert l.count('c') == 0
        assert l.count('g') == 0
        assert l.count('t') == 0
        assert l.count('u') == 0
    assert magnipore_strings[0].count('.') == 19
    assert magnipore_strings[1].count('.') == 13
    assert magnipore_strings[1].count('-') == 6
    