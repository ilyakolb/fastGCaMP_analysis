# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pickle, os
import pandas as pd
import seaborn as sns
from utils import df_from_percent_dict
import scipy.stats as ss
import scikit_posthocs as sp

# plot average FRAP traces and percent traces of all sensors
# takes in roi_trace dict pkls from ingest_FRAP_data_exp2.py
# show stats
# export recovery percents as csvs for analysis in prism

save_figs = True
print_and_save_stats = True


def load_pkl(name):
    with open(name, 'rb') as f:
        return pickle.load(f)

def get_mean_std(trace):
    return (np.mean(trace, axis=0), np.std(trace, axis=0))


def rescale(_trace, _max):
    # _trace = _trace+1
    _t = _trace * 1/(_max-np.min(_trace))
    return _t - _t[0]# - np.min(_t)

def plot_frap_curve(_t, _mean, _std, _color, _ax):
    # _ax.plot(t, mean, colors[i])
    plot_every = 10
    _ax.errorbar(_t[::plot_every], _mean[::plot_every],yerr=_std[::plot_every], fmt='o-', color=_color, markeredgecolor='white', markersize=4)
    
def get_trace_to_plot(trace_array, name):
    '''

    Parameters
    ----------
    trace_array : [nFrames x nTrials] ndarray
        
    name : construct name for saving

    Returns
    -------
        traces_norm_mean: mean FRAP trace
        traces_norm_std: std of FRAP trace
        t: time vector (for plotting)

    '''
    
    # scale to percent
    trace_array*=100
    
    # find min point of trace array
    (traces_mean, traces_std) = get_mean_std(trace_array)
    min_point_idx = traces_mean.argmin()
    
    
    (traces_norm_mean, traces_norm_std) = get_mean_std(trace_array[:,min_point_idx:])
    
    t = np.arange(len(traces_norm_mean))/s_rate
    pd.DataFrame(traces_norm_mean.T).to_csv(r'./analysis/normalized-csvs/' + name + '.csv')
    
    '''
    trace_array_cutoff = trace_array[:,min_point_idx:]
    trace_array_norm = np.zeros_like(trace_array_cutoff)
    for i,t in enumerate(trace_array_cutoff):
        # t = t-t[0]
        prev_mean = np.mean(trace_array[i,:3])# t[0]# np.mean(t[0])
        # stretch trace to fit between 0 and 1
        trace_array_norm[i,:] = rescale(t, prev_mean)# np.interp(t, (0, prev_mean), (0, 1))
        # set baseline to 0 only
        # trace_array_norm[i,:] = t
    
    # plot from lowest point to end
    (traces_norm_mean, traces_norm_std) = get_mean_std(trace_array_norm)
    t = np.arange(len(traces_norm_mean))/s_rate
    
    pd.DataFrame(trace_array_norm.T).to_csv(r'./analysis/normalized-csvs/' + name + '.csv')
    '''
    return(traces_norm_mean, traces_norm_std, t)
    
one_peak = False # load single-peak only (for FRAP curve)

peak_str = '_1peak' if one_peak else '_5peak'

# load 405-bleached data
traces_405 = load_pkl(r'./analysis/roi_traces_norm_regular_405' + peak_str + '.pkl')
percents_405 = load_pkl(r'./analysis/plateau_data_norm_regular_405' + peak_str + '.pkl')

# load 488-bleached data
# traces_488 = load_pkl(r'./analysis/roi_traces_norm_regular_488' + peak_str + '.pkl')
# percents_488 = load_pkl(r'./analysis/plateau_data_norm_regular_488' + peak_str + '.pkl')

s_rate = 50


plt.close('all')
colors = {'10.641': 'gray', '500.686': 'darkred', 
          '500.688': 'blue', 'EGFP.B-actin': 'darkgreen', 'mEm.Cyto':'violet'}
construct_legend = []
percent_fig, axs = plt.subplots(2,1)
percent_fig.set_size_inches([5.62, 6.85])
# remove 604.2 (405 stim)
traces_405.pop('604.2')
percents_405.pop('604.2')

n_constructs = len(traces_405)
for i,construct in enumerate(traces_405.keys()):


    traces_array = np.array(traces_405[construct])

    # plotting
    
    (mean, std, t) = get_trace_to_plot(traces_array, construct)
    
    construct_legend.append(construct + '(n={})'.format(len(percents_405[construct])))
    
    # get rid of spontaneous Ca transient
    if construct == '10.641':
        # print('m ' + str(mean.argmax()))
        mean[mean.argmax() + np.arange(0,5)] =np.nan

    plot_frap_curve(t, mean, std, colors[construct], axs[0])
    
    percents_construct = np.array(percents_405[construct])
    (percent_mean, percent_std) = get_mean_std(percents_construct)
    # axs[1].bar(i, 100*percent_mean, yerr = 100*percent_std, facecolor = 'white', edgecolor = colors[i])
    # axs[1] = sns.swarmplot(x = i*np.ones_like(percents_construct.squeeze()), y = 100*percents_construct.squeeze(), color = colors[i])



df_melted = df_from_percent_dict(percents_405)


sns.swarmplot(x = 'index', 
              y='mean percent', 
              color='gray',
              hue = 'index', 
              palette=colors, 
              data=df_melted,
              order=colors.keys())
# sns.barplot(x='index', y='mean percent', data=df_melted, facecolor='white', edgecolor='black', linewidth=1.5)
sns.boxplot(x='index', 
            y='mean percent', 
            # hue='index', 
            data=df_melted, 
            # palette=colors,
            color='white',
            whis=1.5,
            showfliers=False,
            dodge=False,
            order=colors.keys(),
            width=0.5)
# add 604.2 bleached with 488 traces

axs[1].get_legend().remove()
'''
(mean, std, t) = get_trace_to_plot(np.array(traces_488['604.2']), '604.2_488stim')
(percent_mean, percent_std) = get_mean_std(np.array(percents_488['604.2']))
# axs[0].plot(t, mean, colors[i+1])
# axs[0].fill_between(t, mean + std, mean - std, facecolor=colors[i+1], color=colors[i+1], alpha=0.2)
plot_frap_curve(t, mean, std, colors[i+1], axs[0])
axs[1].bar(i+1, 100*percent_mean, yerr = 100*percent_std, color='black')
construct_legend.append('604.2 (488 bleach) (n={})'.format(len(percents_488['604.2'])))
'''


axs[0].set_ylabel('Recovery (%)')
axs[1].set_ylabel('Immobile fraction (%)')
axs[1].set_xticks(np.arange(0,i+1))
axs[1].set_xticklabels(colors.keys(), rotation = -45, ha="left")

# axs[0].legend(construct_legend)

# unity lines
axs[0].plot(t, 100*np.ones_like(t), 'k--')
axs[1].plot([0, i+1], [0,0], 'k--')

axs[0].set_xlim([-0.2, 2])
axs[0].set_xlabel('Time (s)')

plt.tight_layout()


# save figures
if save_figs:
    percent_fig.savefig(os.path.join('./analysis/normalized', 'percent_change.pdf'))
    

# print stats
if print_and_save_stats:
    
    # stats -- do not use, use prism instead
    '''
    [_,p] = ss.kruskal(*df_melted.groupby('index')['mean percent'].agg(list), nan_policy='omit')
    print('MEANS')
    df_melted.groupby('index')['mean percent'].mean()
    print('STDs')
    df_melted.groupby('index')['mean percent'].std()
    print('KW test: p = {:4f}'.format(p))
    pval_table = sp.posthoc_dunn(df_melted, val_col = 'mean percent', group_col='index')
    print(pval_table)
    '''
    
    # write csvs to be read by prism
    for c_id in df_melted['index'].unique():
        percent_array = df_melted[df_melted['index'] == c_id]['mean percent']# np.array(percents_405[key])
        np.savetxt('./analysis/normalized/recovery_percent_' + c_id + '.csv', percent_array, delimiter='\n')
        # print(key + '= ' + str(percent_array.mean()) + ' +/- ' + str(percent_array.std()))
    