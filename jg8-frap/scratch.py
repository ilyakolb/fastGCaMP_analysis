import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# num_repeats = 4 # number of repeats
# ith_repeat_to_plot = 0 # plot 0th repeat

avg = np.load(r'./analysis/all_roi_avg_data_norm.pkl', allow_pickle=True)
s_rate = 37.037

plt.close('all')
colors = ['gray', 'red', 'blue']
for i,construct in enumerate(avg.keys()):
    avg_to_plot = avg[construct]
    t_combo = np.arange(0,len(avg_to_plot[0]))/s_rate
    [plt.plot(t_combo, a, color=colors[i]) for a in avg[construct]]

plt.xlim([.3, 2])
plt.ylim([-.5, .1])
plt.plot(t_combo, np.zeros_like(t_combo), 'k--')
plt.show()