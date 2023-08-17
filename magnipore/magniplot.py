#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import multiprocessing as mp
import os
import random
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from magnipore.__init__ import __version__


def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description='Creating plots (MeDAS, etc) for a given .magnipore file.',
        prog='Magnipore',
    )
    parser.add_argument('magnipore', type=str, help='Magnipore-style output')
    parser.add_argument('outdir', type=str, help='Path to write plots.')
    parser.add_argument("label_first_sample", type = str, help='Name of the sample or pipeline run')
    parser.add_argument("label_sec_sample", type = str, help='Name of the sample or pipeline run')
    parser.add_argument('-f', '--fontsize', type=int, default=18, help='Fontsize for plots')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of processes to use to create plots')
    parser.add_argument('-nl', '--num_lines', type=int, default=None, help='Providing the number of lines in file speeds up the process.')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s' + f' {__version__}')
    return parser.parse_args()

def loadPandas(magnipore_file : str, num_lines : int = None) -> pd.DataFrame:
    # sample data for given size, if it got too large
    # reduces runtime and prevent the kernel from killing the process
    if num_lines is None:
        with open(magnipore_file, "rb") as f:
            num_lines = sum(1 for _ in f)
    s = min(1000000, num_lines) - 1 # -1 because we substract the header line
    skip = sorted(random.sample(range(1, num_lines), num_lines-s))
    columns = ['strand', 'td_score', 'kl_divergence', 'signal_type', 'signal_mean_1', 'signal_std_1', 'n_reads_1', 'signal_mean_2', 'signal_std_2', 'n_reads_2']
    data = pd.read_csv(magnipore_file, sep='\t', usecols=columns, header=0, skiprows=skip)
    # prepare columns for plots
    data['Mean Distance'] = abs(data['signal_mean_1'] - data['signal_mean_2'])
    data['Avg Stdev'] = (data['signal_std_1'] + data['signal_std_2'])/2
    data['Significant'] = np.where(data['td_score'] >= 1.0, True, False)
    data['Low Coverage (<10)'] = np.where((data['n_reads_1'] < 10) | (data['n_reads_2'] < 10), True, False)
    data = data.rename(columns={'strand' : 'Strand', 'signal_type' : 'Sequence Context', 'td_score' : 'TD Score', 'kl_divergence' : 'KL Divergence'})
    data['Context, Significant'] = pd.Series(data.reindex(['Sequence Context', 'Significant'], axis='columns').astype('str').values.tolist()).str.join(', ')
    # remove unused columns
    data.drop(columns=['signal_mean_1', 'signal_std_1', 'n_reads_1', 'signal_mean_2', 'signal_std_2', 'n_reads_2'], inplace=True)
    # drop NANs
    data.dropna(how='any', inplace=True, ignore_index = True) # ignore_index essential for plots in method `plotScores`
    return data

def callbackError(error):
    print(f'Error in multiprocessing signal comparison: {error}')

def plotStatistics(data : pd.DataFrame, outdir : str, label_first_sample : str, label_sec_sample : str, threads : int, fontsize : int) -> None:
    pool = mp.Pool(threads)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    # Mean Dist vs Std Avg plot
    print(f'Plotting Mean vs Stdev of {len(data.index)} positions')
    pool.apply_async(plotMeanDistAvgStd, args=(data, outdir, label_first_sample, label_sec_sample, fontsize), error_callback=callbackError)

    print(f'Plotting Mean vs Stdev of {len(data.index)} positions excluding low coverage positions')
    pool.apply_async(plotMeanDistAvgStd, args=(data[data['Low Coverage (<10)'] == False], outdir, label_first_sample, label_sec_sample, fontsize, 'c10'), error_callback=callbackError)

    # plot MeanDistStdAvg with coverage
    print(f'Plotting Mean vs Stdev with coverage markers if {len(data.index)} positions')
    pool.apply_async(plotMeanDistAvgStdCov, args=(data, outdir, label_first_sample, label_sec_sample, fontsize), error_callback=callbackError)

    # plot scores
    print(f'Plotting TD score and KL divergence')
    pool.apply_async(plotScores, args=(data, outdir, label_first_sample, label_sec_sample), error_callback=callbackError)

    pool.close()
    pool.join()

def plotScores(data : pd.DataFrame, working_dir : str, label_first_sample : str, label_sec_sample : str) -> None:
    
    colors = {
        'mod, False':'wheat',
        'mod, True':'darkorange',
        'mut, False':'skyblue',
        'mut, True':'darkblue'}

    plt.figure(figsize = (12,8), dpi=300)
    plt.title(f'TD score for all positions\n{label_first_sample} vs. {label_sec_sample}')
    # otherwise logscale range is infinite with lower bound: -infinity
    sns.histplot(data=data[data['TD Score']>0], x='TD Score', hue=data['Context, Significant'], log_scale=(True, True), multiple="stack", palette=colors)
    plt.grid(True,  'both', 'both', alpha=0.6, linestyle='--')
    plt.tight_layout()
    plt.savefig(os.path.join(working_dir, f'{label_first_sample}_{label_sec_sample}_td_score.png'))
    plt.savefig(os.path.join(working_dir, f'{label_first_sample}_{label_sec_sample}_td_score.pdf'))
    plt.close()

    plt.figure(figsize = (12,8), dpi=300)
    plt.title(f'Kullback-Leibler divergence for all positions\n{label_first_sample} vs. {label_sec_sample}')
    # otherwise logscale range is infinite with lower bound: -infinity
    sns.histplot(data=data[data['KL Divergence']>0], x='KL Divergence', hue=data['Context, Significant'], log_scale=(True, True), multiple="stack", palette=colors)
    plt.grid(True,  'both', 'both', alpha=0.6, linestyle='--')
    plt.tight_layout()
    plt.savefig(os.path.join(working_dir, f'{label_first_sample}_{label_sec_sample}_kl_div.png'))
    plt.savefig(os.path.join(working_dir, f'{label_first_sample}_{label_sec_sample}_kl_div.pdf'))
    plt.close()
    
def plotMeanDistAvgStd(data : pd.DataFrame, working_dir : str, label_first_sample : str, label_sec_sample : str, fontsize : int, suffix : str = None) -> None:
    
    marker = lambda mut_context: 'D' if mut_context == 'mut' else 'o'
    color = lambda mut_context: 'blue' if mut_context == 'mut' else '#d95f02' 

    ### Mean Dist vs Std Avg plot
    plt.figure(figsize = (12,12), dpi=300)
    plt.rcParams.update({
        'font.size': fontsize,
        })
    label1 = label_first_sample.replace("_", " ")
    label2 = label_sec_sample.replace("_", " ")

    g = sns.JointGrid(x='Mean Distance', y='Avg Stdev', data=data, hue='Sequence Context', marginal_ticks=True, palette=['blue', '#d95f02'], hue_order=['mut', 'mod'], height = 10)
    g.plot_joint(func=sns.scatterplot, s = 8)
    g.ax_joint.cla()
    for _, row in data.iterrows():
        g.ax_joint.plot(row['Mean Distance'], row['Avg Stdev'], color = color(row['Sequence Context']), marker = marker(row['Sequence Context']), markersize=3, alpha = 0.6)
    g.fig.suptitle(f'{len(data.index)} compared bases mean distance against\naverage standard deviation\n{label1} and {label2}', y=0.98)
    g.ax_joint.grid(True, 'both', 'both', alpha = 0.4, linestyle = '--', linewidth = 0.5)
    g.plot_marginals(sns.histplot, binwidth = 0.005, kde = True, linewidth = 0)
    g.ax_marg_x.grid(True, 'both', 'both', alpha = 0.4, linestyle = '-', linewidth = 0.5)
    g.ax_marg_y.grid(True, 'both', 'both', alpha = 0.4, linestyle = '-', linewidth = 0.5)

    lims = np.array([
        [-.02, max(data['Mean Distance']) + 0.1],
        [-.02, max(data['Avg Stdev']) + 0.1]
    ])
    y1 = np.arange(min(lims[:, 0]), max(lims[:, 1]) + 0.01, 0.01)
    y2 = np.repeat(max(lims[:, 1]), len(y1))
    g.ax_joint.fill_between(y1, lims[1,0], y1, color='#1b9e77', alpha=0.15)
    g.ax_joint.fill_between(y1, y1, y2, color='#7570b3', alpha=0.15)
    sign = mlines.Line2D([], [], color='#1b9e77', marker='s', linestyle='None', markersize=10, label=r'significant, TD$\geq$1')
    insign = mlines.Line2D([], [], color='#7570b3', marker='s', linestyle='None', markersize=10, label=r'insignificant, TD$<$1')

    g.ax_joint.set_xlim(tuple(lims[0]))
    g.ax_joint.set_ylim(tuple(lims[1]))
    g.ax_joint.set_ylabel('Average standard deviation')
    g.ax_joint.set_xlabel('Mean distance')
    legend_mut = mlines.Line2D([], [], color='blue', marker='D', linestyle='None', markersize=10, label='mutation')
    legend_mod = mlines.Line2D([], [], color='#d95f02', marker='o', linestyle='None', markersize=10, label='matching reference')
    handles = [legend_mut, legend_mod, sign, insign]
    g.ax_joint.legend(handles = handles, fontsize = 'small', framealpha = 0.3)
    g.fig.subplots_adjust(top=0.95)

    plt.savefig(os.path.join(working_dir, f'{label_first_sample}_{label_sec_sample}_{suffix + "_" if suffix is not None else ""}MeDAS.png'))
    plt.savefig(os.path.join(working_dir, f'{label_first_sample}_{label_sec_sample}_{suffix + "_" if suffix is not None else ""}MeDAS.pdf'))

    plt.close()

def plotMeanDistAvgStdCov(data : pd.DataFrame, working_dir : str, label_first_sample : str, label_sec_sample : str, fontsize : int) -> None:
    ### Mean Dist vs Std Avg plot
    plt.figure(figsize = (12,12), dpi=300)
    plt.rcParams.update({
        'font.size': fontsize,
        })
    label1 = label_first_sample.replace("_", " ")
    label2 = label_sec_sample.replace("_", " ")
    sns.relplot(data = data, x='Mean Distance', y='Avg Stdev', hue='Sequence Context', style='Low Coverage (<10)', palette=['#d95f02', 'blue'], height=10)
    plt.title(f'{len(data.index)} compared bases mean distance against\naverage standard deviation\n{label1} and {label2}', y=0.98)
    plt.grid(True, 'both', 'both', alpha = 0.4, linestyle = '--', linewidth = 0.5)

    lims = np.array([
        [-.02, max(data['Mean Distance']) + 0.1],
        [-.02, max(data['Avg Stdev']) + 0.1]
    ])

    y1 = np.arange(min(lims[:, 0]), max(lims[:, 1]) + 0.01, 0.01)
    y2 = np.repeat(max(lims[:, 1]), len(y1))

    plt.fill_between(y1, lims[1,0], y1, color = '#1b9e77', alpha = 0.15)
    plt.fill_between(y1, y1, y2, color = '#7570b3', alpha = 0.15)

    plt.xlim(tuple(lims[0]))
    plt.ylim(tuple(lims[1]))

    plt.xlabel('Mean distance')
    plt.ylabel('Average standard deviation')
    plt.tight_layout()
    plt.savefig(os.path.join(working_dir, f'{label_first_sample}_{label_sec_sample}_MeDAS_coverage.png'))
    plt.savefig(os.path.join(working_dir, f'{label_first_sample}_{label_sec_sample}_MeDAS_coverage.pdf'))

    plt.ylim(bottom = min(data['Avg Stdev']))
    plt.yscale('log')
    plt.savefig(os.path.join(working_dir, f'{label_first_sample}_{label_sec_sample}_MeDAS_coverage_log.png'))
    plt.savefig(os.path.join(working_dir, f'{label_first_sample}_{label_sec_sample}_MeDAS_coverage_log.pdf'))

    plt.close()

def main() -> None:
    args = parse()
    plotStatistics(
        loadPandas(args.magnipore, args.num_lines),
        args.outdir,
        args.label_first_sample,
        args.label_sec_sample,
        args.threads,
        args.fontsize
        )

if __name__ == '__main__':
    main()