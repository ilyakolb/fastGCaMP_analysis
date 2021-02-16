import numpy as np
import matplotlib.pyplot as plt
import pickle, os
import pandas as pd


# plot average FRAP traces and percent traces of all sensors
# takes in roi_trace dict pkls from ingest_FRAP_data_exp2.py
# show stats
# export recovery percents as csvs for analysis in prism

save_figs = 0

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
    
one_peak = True # load single-peak only (for FRAP curve)

peak_str = '_1peak' if one_peak else ''

# load 405-bleached data
traces_405 = load_pkl(r'./analysis/roi_traces_norm_regular_405' + peak_str + '.pkl')
percents_405 = load_pkl(r'./analysis/plateau_data_norm_regular_405' + peak_str + '.pkl')

# load 488-bleached data
traces_488 = load_pkl(r'./analysis/roi_traces_norm_regular_488' + peak_str + '.pkl')
percents_488 = load_pkl(r'./analysis/plateau_data_norm_regular_488' + peak_str + '.pkl')

s_rate = 50


plt.close('all')
colors = ['gray', 'darkred', 'blue', 'darkgreen', 'violet']
construct_legend = []
percent_fig, axs = plt.subplots(2,1)
percent_fig.set_size_inches([5.62, 6.85])
# remove 604.2 (405 stim) from traces_405
traces_405.pop('604.2')

n_constructs = len(traces_405)
for i,construct in enumerate(traces_405.keys()):
    traces_array = np.array(traces_405[construct])

    (mean, std, t) = get_trace_to_plot(traces_array, construct)
    
    construct_legend.append(construct + '(n={})'.format(len(traces_array)))
    
    plot_frap_curve(t, mean, std, colors[i], axs[0])

    percents_construct = percents_405[construct]
    x = np.arange(1,percents_construct[0].shape[0]+1)
    (percent_mean, percent_std) = get_mean_std(np.array(percents_construct))
    axs[1].bar(i, 100*percent_mean, yerr = 100*percent_std, facecolor = 'white', edgecolor = colors[i])


# add 604.2 bleached with 488 traces

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
axs[1].set_ylabel('Recovery (%)')
axs[1].set_xticks(np.arange(0,i+1))
axs[1].set_xticklabels(construct_legend, rotation = -45, ha="left")

axs[0].legend(construct_legend)

# unity lines
axs[0].plot(t, 100*np.ones_like(t), 'k--')
axs[1].plot([0, i+1], [100,100], 'k--')

axs[0].set_xlim([-0.2, 8])
axs[0].set_xlabel('Time (s)')

plt.tight_layout()




### 405 vs 488 comparison of 10.641
f, axs = plt.subplots(2,1)
f.set_size_inches([4.99, 6.56])

(mean, std, t) = get_trace_to_plot(np.array(traces_405['10.641']), '10.641')
(percent_mean_405, percent_std_405) = get_mean_std(np.array(percents_405['10.641']))
plot_frap_curve(t, mean, std, 'darkviolet', axs[0])
# axs[0].plot(t, mean, color = 'darkviolet')
# axs[0].fill_between(t, mean + std, mean - std, facecolor='darkviolet', color='darkviolet', alpha=0.2)
axs[1].bar(0, 100*percent_mean_405, yerr=100*percent_std_405, color='darkviolet')
construct_legend = ['10.641 (405 bleach) (n={})'.format(len(percent_std_405))]


(mean, std, t) = get_trace_to_plot(np.array(traces_488['10.641']), '10.641_488stim')
(percent_mean_488, percent_std_488) = get_mean_std(np.array(percents_488['10.641']))
plot_frap_curve(t, mean, std, 'cyan', axs[0])
# axs[0].plot(t, mean, color = 'cyan')
# axs[0].fill_between(t, mean + std, mean - std, facecolor='cyan', color='cyan', alpha=0.2)
axs[1].bar(1, 100*percent_mean_488, yerr=100*percent_std_488, color='cyan')
axs[1].set_xticks([0,1])
construct_legend = ['10.641 (405 bleach) (n={})'.format(len(traces_405['10.641'])), '10.641 (488 bleach) (n={})'.format(len(traces_488['10.641']))]
axs[0].legend(construct_legend)
axs[1].set_xticklabels(construct_legend, rotation = -45, ha="left")
axs[0].set_ylabel('Recovery (%)')
axs[0].set_xlabel('Time (s)')
axs[1].set_ylabel('Recovery (%)')
plt.tight_layout()
plt.show()

# save figures
if save_figs:
    percent_fig.savefig(os.path.join('./analysis/normalized', 'percent_change.pdf'))
    f.savefig(os.path.join('./analysis/normalized', '405_vs_488_6s.pdf'))
    

# print stats
print('recovery percents (mean +/- std')
for key in percents_405.keys():
    percent_array = np.array(percents_405[key])
    np.savetxt('./analysis/normalized/recovery_percent_' + key + '.csv', percent_array, delimiter='\n')
    print(key + '= ' + str(percent_array.mean()) + ' +/- ' + str(percent_array.std()))
