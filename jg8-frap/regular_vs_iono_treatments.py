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

def get_trace_to_plot(trace_array):
    (traces_mean, traces_std) = get_mean_std(trace_array)
    
    # plot from lowest point to end
    min_point_idx = traces_mean.argmin()
    traces_mean = traces_mean[min_point_idx:]
    traces_std = traces_std[min_point_idx:]
    t = np.arange(len(traces_mean))/s_rate
    return(traces_mean, traces_std, t)


plt.close('all')
traces_regular = load_pkl(r'./analysis/roi_traces_norm_regular_405.pkl')
traces_iono = load_pkl(r'./analysis/roi_traces_norm_iono_405.pkl')

percents_regular = load_pkl(r'./analysis/plateau_data_norm_regular_405.pkl')
percents_iono = load_pkl(r'./analysis/plateau_data_norm_iono_405.pkl')

s_rate = 50

all_constructs = ['10.641', '604.2','500.688','500.686'] # '500.688', '604.2', '500.686'

for construct in all_constructs:
    
    traces_regular_construct = traces_regular[construct]
    traces_iono_construct = traces_iono[construct]
    
    colors = ['gray', 'red']
    
    plt.figure()
    plt.subplot(2,1,1)
    t = np.arange(traces_regular_construct[0].shape[0])/s_rate
    plt.legend(['regular', 'iono'])
    (reg_mean, reg_std, t) = get_trace_to_plot(traces_regular_construct)
    plt.plot(t, reg_mean, colors[0])
    plt.fill_between(t, reg_mean + reg_std, reg_mean - reg_std, facecolor=colors[0], color=colors[0], alpha=0.2)
    
    (iono_mean, iono_std, t) = get_trace_to_plot(traces_iono_construct)
    plt.plot(t, iono_mean, colors[1])
    plt.fill_between(t, iono_mean + iono_std, iono_mean - iono_std, facecolor=colors[1], color=colors[1], alpha=0.2)
    
    plt.plot(t, np.zeros_like(t), 'k--')
    plt.legend(['regular', 'iono'])
    plt.show()
    plt.ylim([-0.4, 0.4])
    plt.title(construct)
    plt.xlabel('Time (s)')
    
    # percent change plot
    percents_regular_construct = percents_regular[construct]
    percents_iono_construct = percents_iono[construct]
    plt.subplot(2,1,2)
    x = np.arange(1,percents_regular_construct[0].shape[0]+1)
    (percent_reg_mean, percent_reg_std) = get_mean_std(np.array(percents_regular_construct))
    (percent_iono_mean, percent_iono_std) = get_mean_std(np.array(percents_iono_construct))
    plt.errorbar(x, 100*percent_reg_mean, 100*percent_reg_std, color=colors[0])
    plt.errorbar(x, 100*percent_iono_mean, 100*percent_iono_std, color=colors[1])
    plt.legend(['regular', 'iono'])
    plt.xlabel('stim number')
    plt.ylabel('Percent change (norm.)')