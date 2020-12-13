import numpy as np
import matplotlib.pyplot as plt
import pickle


# plot average FRAP traces and percent traces of all sensors
# takes in roi_trace dict pkls from ingest_FRAP_data_exp2.py


def load_pkl(name):
    with open(name, 'rb') as f:
        return pickle.load(f)

def get_mean_std(trace):
    return (np.mean(trace, axis=0), np.std(trace, axis=0))
    
def get_trace_to_plot(trace_array):
    (traces_mean, traces_std) = get_mean_std(trace_array)
    
    # plot from lowest point to end
    min_point_idx = traces_mean.argmin()
    traces_mean = traces_mean[min_point_idx:]
    traces_std = traces_std[min_point_idx:]
    t = np.arange(len(traces_mean))/s_rate
    return(traces_mean, traces_std, t)
    
# load 405-bleached data
traces_405 = load_pkl(r'./analysis/roi_traces_norm_regular_405.pkl')
percents_405 = load_pkl(r'./analysis/plateau_data_norm_regular_405.pkl')

# load 488-bleached data
traces_488 = load_pkl(r'./analysis/roi_traces_norm_regular_488.pkl')
percents_488 = load_pkl(r'./analysis/plateau_data_norm_regular_488.pkl')

s_rate = 50


plt.close('all')
colors = ['gray', 'red', 'blue', 'cyan', 'green']
construct_legend = []
_, axs = plt.subplots(2,1)

n_constructs = len(traces_405)
for i,construct in enumerate(traces_405.keys()):
    traces_array = np.array(traces_405[construct])

    (mean, std, t) = get_trace_to_plot(traces_array)
    
    construct_legend.append(construct)
    
    axs[0].plot(t, mean, colors[i])
    axs[0].fill_between(t, mean + std, mean - std, facecolor=colors[i], color=colors[i], alpha=0.2)
    axs[0].set_ylim([-0.4, 0.4])
    
    percents_construct = percents_405[construct]
    x = np.arange(1,percents_construct[0].shape[0]+1)
    (percent_mean, percent_std) = get_mean_std(np.array(percents_construct))
    axs[1].errorbar(x-n_constructs/1000+i/500, 100*percent_mean, 100*percent_std, color=colors[i], capsize=5, marker='o')


# add 604.2 bleached with 488 traces
(mean, std, t) = get_trace_to_plot(traces_488['604.2'])
(percent_mean, percent_std) = get_mean_std(np.array(percents_488['604.2']))
axs[0].plot(t, mean, colors[i+1])
axs[0].fill_between(t, mean + std, mean - std, facecolor=colors[i+1], color=colors[i+1], alpha=0.2)
axs[0].set_ylim([-0.4, 0.4])
axs[1].errorbar(x-n_constructs/1000+(i+1)/500, 100*percent_mean, 100*percent_std, color=colors[i+1], capsize=5, marker='o')

construct_legend.append('604.2 (488 bleach)')

# 0 lines
axs[0].plot(t, np.zeros_like(t), 'k--')
axs[1].plot(x, np.zeros_like(x), 'k--')


axs[0].set_xlabel('time (s)')
axs[0].set_ylabel('F (norm.)')
axs[1].set_xlabel('stim number')
axs[1].set_ylabel('Percent change (norm.)')
axs[0].legend(construct_legend)


### 405 vs 488 comparison of 10.641
_, axs = plt.subplots(2,1)

(mean, std, t) = get_trace_to_plot(np.array(traces_405['10.641']))
(percent_mean_405, percent_std_405) = get_mean_std(np.array(percents_405['10.641']))
axs[0].plot(t, mean, color = 'darkviolet')
axs[0].fill_between(t, mean + std, mean - std, facecolor='darkviolet', color='darkviolet', alpha=0.2)
axs[1].errorbar(x-n_constructs/1000+(i+1)/500, 100*percent_mean_405, 100*percent_std_405, color='darkviolet', capsize=5, marker='o')

(mean, std, t) = get_trace_to_plot(np.array(traces_488['10.641']))
(percent_mean_488, percent_std_488) = get_mean_std(np.array(percents_488['10.641']))
axs[0].plot(t, mean, color = 'cyan')
axs[0].fill_between(t, mean + std, mean - std, facecolor='cyan', color='cyan', alpha=0.2)
axs[0].set_ylim([-0.4, 0.4])


axs[1].errorbar(x-n_constructs/1000+(i+1)/500, 100*percent_mean_488, 100*percent_std_488, color='cyan', capsize=5, marker='o')
axs[0].legend(['10.641 (405 bleach)', '10.641 (488 bleach)'])

plt.show()