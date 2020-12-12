import numpy as np
import matplotlib.pyplot as plt
import pickle


# plot average FRAP traces of all sensors
# takes in roi_trace dict pkls from ingest_FRAP_data_exp2.py


def load_pkl(name):
    with open(name, 'rb') as f:
        return pickle.load(f)
    
traces = load_pkl(r'./analysis/roi_traces_norm_regular_405.pkl')
s_rate = 50


# plt.close('all')
colors = ['gray', 'red', 'blue', 'cyan']
construct_legend = []
_, ax = plt.subplots()
for i,construct in enumerate(traces.keys()):
    traces_array = np.array(traces[construct])
    
    traces_mean = np.mean(traces_array, axis=0)
    traces_std = np.std(traces_array, axis=0)
    
    t = np.arange(traces_array.shape[1])/s_rate
    
    # [plt.plot(t, c, color = colors[i]) for c in traces[construct] ]
    
    construct_legend.append(construct)
    plt.plot(t, traces_mean, colors[i])
    plt.fill_between(t, traces_mean + traces_std, traces_mean - traces_std, facecolor=colors[i], color=colors[i], alpha=0.2)

ax.legend(construct_legend)
# plt.xlim([.3, 2])
# plt.ylim([-.5, .1])
# plt.plot(t_combo, np.zeros_like(t_combo), 'k--')
# 
plt.show()
