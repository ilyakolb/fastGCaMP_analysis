import numpy as np
import matplotlib.pyplot as plt
import pickle


# plot average FRAP traces of all sensors
# takes in roi_trace dict pkls from ingest_FRAP_data_exp2.py


def load_pkl(name):
    with open(name, 'rb') as f:
        return pickle.load(f)

def get_mean_std(trace):
    return (np.mean(trace, axis=0), np.std(trace, axis=0))
    
traces = load_pkl(r'./analysis/roi_traces_norm_regular_405.pkl')
percents = load_pkl(r'./analysis/plateau_data_norm_regular_405.pkl')
s_rate = 50


plt.close('all')
colors = ['gray', 'red', 'blue', 'cyan']
construct_legend = []
_, axs = plt.subplots(2,1)

n_constructs = len(traces)
for i,construct in enumerate(traces.keys()):
    traces_array = np.array(traces[construct])

    (traces_mean, traces_std) = get_mean_std(traces_array)
    
    t = np.arange(traces_array.shape[1])/s_rate
    
    # [plt.plot(t, c, color = colors[i]) for c in traces[construct] ]
    
    construct_legend.append(construct)
    
    axs[0].plot(t, traces_mean, colors[i])
    axs[0].fill_between(t, traces_mean + traces_std, traces_mean - traces_std, facecolor=colors[i], color=colors[i], alpha=0.2)
    axs[0].set_ylim([-0.4, 0.4])
    
    percents_construct = percents[construct]
    x = np.arange(1,percents_construct[0].shape[0]+1)
    (percent_mean, percent_std) = get_mean_std(np.array(percents_construct))
    axs[1].errorbar(x-n_constructs/1000+i/500, 100*percent_mean, 100*percent_std, color=colors[i], capsize=5, marker='o')
    axs[1].set_xlabel('stim number')
    axs[1].set_ylabel('Percent change (norm.)')
    
axs[0].plot(t, np.zeros_like(t), 'k--')
axs[0].legend(construct_legend)
# plt.xlim([.3, 2])
# plt.ylim([-.5, .1])
# plt.plot(t_combo, np.zeros_like(t_combo), 'k--')
# 
plt.show()
