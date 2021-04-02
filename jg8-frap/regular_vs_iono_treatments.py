import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import pickle, os
from scipy import stats
from plot_averages_2 import plot_frap_curve
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd
from utils import df_from_percent_dict
import seaborn as sns

matplotlib.rcParams['pdf.fonttype'] = 42
# plot average FRAP traces of all sensors
# takes in roi_trace dict pkls from ingest_FRAP_data_exp2.py

plt.close('all')

def load_pkl(name):
    with open(name, 'rb') as f:
        return pickle.load(f)

def get_mean_std(trace):
    return (np.mean(trace, axis=0), np.std(trace, axis=0))

def rescale(_trace, _max):
    # _trace = _trace+1
    _t = _trace * 1/(_max-np.min(_trace))
    return _t - _t[0]# - np.min(_t)


def get_trace_to_plot(trace_array):
    (traces_mean, traces_std) = get_mean_std(trace_array)
    min_point_idx = traces_mean.argmin()
    
    trace_array_cutoff = trace_array[:,min_point_idx:]

    (traces_norm_mean, traces_norm_std) = get_mean_std(trace_array_cutoff)
    t = np.arange(len(traces_norm_mean))/s_rate
    return(traces_norm_mean, traces_norm_std, t)
    
    '''
     trace_array_norm = np.zeros_like(trace_array_cutoff)
     for i,t in enumerate(trace_array_cutoff):
        prev_mean = np.mean(trace_array[i,:3])# t[0]# np.mean(t[0])
        # stretch trace to fit between 0 and 1
        trace_array_norm[i,:] = rescale(t, prev_mean)# np.interp(t, (0, prev_mean), (0, 1))
    '''



save_figs = True
one_peak = False
peak_str = '_1peak' if one_peak else '_10peak'

traces_regular = load_pkl(r'./analysis/roi_traces_norm_regular_405' + peak_str + '.pkl')
traces_iono = load_pkl(r'./analysis/roi_traces_norm_iono_405' + peak_str + '.pkl')

percents_regular = load_pkl(r'./analysis/plateau_data_norm_regular_405' + peak_str + '.pkl')
percents_iono = load_pkl(r'./analysis/plateau_data_norm_iono_405' + peak_str + '.pkl')

# load 488-bleached data
traces_regular_488 = load_pkl(r'./analysis/roi_traces_norm_regular_488' + peak_str + '.pkl')
percents_regular_488 = load_pkl(r'./analysis/plateau_data_norm_regular_488' + peak_str + '.pkl')

s_rate = 50

all_constructs = ['10.641','500.688','500.686'] # '500.688', '604.2', '500.686'

colors = {'regular': 'gray', 'iono': 'red'}

df_percents_regular = df_from_percent_dict(percents_regular)
df_percents_regular['exp'] = 'regular'
df_percents_iono = df_from_percent_dict(percents_iono)
df_percents_iono['exp'] = 'iono'
df_percents = pd.concat([df_percents_regular, df_percents_iono])
for construct in all_constructs:
    
    if construct == '604.2': # if fungal GCaMP, use 488 data
        traces_regular_construct = traces_regular_488[construct]
        
    else:
        traces_regular_construct = traces_regular[construct]

    
    traces_iono_construct = traces_iono[construct]
    
    
    construct_legend = ['regular (488 bleach)' if construct == '604.2' else 'regular', 'iono']
    
    construct_legend[0] += (' (n={})'.format(len(traces_regular_construct)))
    construct_legend[1] += (' (n={})'.format(len(traces_iono_construct)))
    
    f = plt.figure(figsize = [4.5, 4.5])
    f.suptitle(construct)
    ax1 = plt.subplot()
    t = np.arange(traces_regular_construct[0].shape[0])/s_rate
    (reg_mean, reg_std, t_reg) = get_trace_to_plot(np.array(traces_regular_construct))
    
    # ax1.plot(t, reg_mean, colors[0])
    # ax1.fill_between(t, reg_mean + reg_std, reg_mean - reg_std, facecolor=colors[0], color=colors[0], alpha=0.2)
    
    (iono_mean, iono_std, t_iono) = get_trace_to_plot(np.array(traces_iono_construct))
    
    # array sizes can be different, shorten
    shortest_length = np.min([t_reg.shape[0], t_iono.shape[0]])
    reg_mean = reg_mean[:shortest_length]
    reg_std = reg_std[:shortest_length]
    t_reg = t_reg[:shortest_length]
    iono_mean = iono_mean[:shortest_length]
    iono_std = iono_std[:shortest_length]
    t_iono = t_iono[:shortest_length]
    
    # remove Ca transient
    if construct == '10.641':
        reg_mean[reg_mean.argmax() + np.arange(0,5)] =np.nan
        reg_mean[reg_std.argmax() + np.arange(0,5)] =np.nan
        
    plot_frap_curve(t_reg, 100*reg_mean, 100*reg_std, colors['regular'], ax1)
    plot_frap_curve(t_iono, 100*iono_mean, 100*iono_std, colors['iono'], ax1)
    ax1.legend(construct_legend)
    ax1.plot(t, 100*np.ones_like(t), 'k--')
    
    plt.ylim([80, 105])
    ax1.set_xlabel('Time (s)')
    
    
    # percent change bar plots    
    ax2 = inset_axes(ax1, width="30%", height="40%", loc=4, borderpad=3)
    
    df_barplot = df_percents[df_percents['index']==construct]

    sns.swarmplot(x = 'exp', 
              y='mean percent', 
              color='black',
              data=df_barplot,
              order=colors.keys())
    
    sns.boxplot(x='exp', 
                y='mean percent', 
                data=df_barplot, 
                # palette=colors,
                color='white',
                whis=False,
                showfliers=False,
                dodge=False,
                hue='exp',
                palette=colors,
                order=colors.keys(),
                width=0.5)
    
    
    
    ax2.set_xlabel(None)
    ax2.set_ylim([-1,5])
    ax2.set_ylabel('Immobile fraction (%)')
    # plt.tight_layout()
    ax2.get_legend().remove()
    # stats
    (_,pval) = stats.ttest_ind(df_barplot[df_barplot['exp'] == 'regular']['mean percent'], df_barplot[df_barplot['exp'] == 'iono']['mean percent'], equal_var = False, nan_policy='omit')
    print(construct + ': regular vs iono: p = ' + str(pval))
    if save_figs:
        f.savefig(os.path.join('./analysis/normalized', 'iono_' + construct + '.pdf'))

