# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 12:59:13 2020

@author: kolbi

NOTES:
    Second FRAP experiment (12/10/20)
    
USAGE:
    env: conda activate base2
    run from D:\pythonTesting\jg8-frap
    For un-normalized, plot N=20 pulses to get fraction change plot
    For normalized, plot N=5 pulses to get uncontaminated curves
    
NOTES:
        # roi1: frap pixel
        # roi2: frap roi
        # roi3: control roi within cell
        # roi4: neighboring cell roi <<< NOT PRESENT IN second experiment
        
        # @todo: check sampling rate
        
@todo: inspect video where resistant fraction is visible?

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import os, pickle

def decode_filename(fname):
    '''
    decode csv filename like: 10.641.20201210.006.regular.stim405.csv

    Parameters
    ----------
    fname : str
        csv str in the format {sensor}.{exp date}.{cellnum}.{regular/iono}.{stim405/stim488}.csv

    Returns
    -------
        construct: MUST HAVE . IN CONSTRUCT NAME E.G. '10.641'
        expdate: e.g. 20201211
        cellnum: e.g. 005
        imgbuffer: 'iono' or 'regular'
        stim_laser: '405' or '488'

    '''
    fname_split = fname.split('.')
    construct = fname_split[0] + '.' + fname_split[1]
    expdate = fname_split[2]
    cellnum = fname_split[3]
    imgbuffer = fname_split[4]
    stimlaser = fname_split[5].lstrip('stim')
    return (construct, expdate, cellnum, imgbuffer, stimlaser)
    

plt.close('all')
num_peaks_to_plot = 10 # 40 # 'all' to plot all
length_to_plot = 200
samples_pre_stim = 10 # was 10
peak_thresh = 20000
plateau_start_idx = 100
plateau_end_idx = 20 # was 10
s_rate = 50 # DOUBLE-CHECK

save_figs       = False
save_data       = True
normalize_roi   = True
keep_figs_open  = False # True to keep all generated figures open. Memory errors if too many open

bleachlaser_condition = 'stim405' # or 'stim488'
solution_condition = 'regular' # or 'iono'
all_constructs = ['604.2', '10.641']# [ '10.641', '604.2','500.688','500.686'] # '500.688', '604.2', '500.686'
# all_cell_num = ['001.iono']# ['001', '002', '003', '004', '005']


# directory of combined data
combo_dir = r'Z:\ilya\code\fastGCaMP_analysis\jg8-frap\data\combined' 
csv_filenames_all = os.listdir(combo_dir)

findmatches = lambda x: any([1 for c in all_constructs if c in x]) and solution_condition in x and bleachlaser_condition in x

csv_filenames = list(filter(findmatches, csv_filenames_all))


plateau_data = {} # dict for storing all plateau data
all_roi_avg_data = {} # 

exp_folder = 'exp3_20201211'# 'exp3_20201211'

# for construct in all_constructs:

for csv_filename in csv_filenames:

    (construct, expdate, cellnum, imgbuffer, stim_laser) = decode_filename(csv_filename)    
    # plateau_data[construct] = []
    # all_roi_avg_data[construct] = []
    data = pd.read_csv(os.path.join(combo_dir, csv_filename)) 
    t = data['Time [s]']

    roi1 = data['#1 (CSU (488))'].values
    roi2 = data['#2 (CSU (488))'].values
    roi3 = data['#3 (CSU (488))'].values
    # roi4 = data['#4 ( 1)'].values
    
    if normalize_roi:
        roi2 = roi2 / roi3
        roi3 = roi3 / roi3
        # roi4 = roi4 / roi3
        folder_name = 'normalized'
    else:
        folder_name = 'unnormalized'
        
    
    all_rois = (roi1, roi2, roi3) #, roi4)
    idx_peaks, _ = find_peaks(roi1, height=peak_thresh, distance=20)
    
    print('total number of stims: {}'.format(len(idx_peaks)))
    if num_peaks_to_plot != 'all':
        idx_peaks = idx_peaks[:num_peaks_to_plot]
    
    # adjust peak indices to find rightmost peak
    for i,_ in enumerate(idx_peaks):
        idx_peaks[i] += np.where(roi1[idx_peaks[i]:idx_peaks[i] +10] > peak_thresh)[0][-1]
    
    num_stims = len(idx_peaks)
    
    fig, raw_ax = plt.subplots()
    raw_ax.plot(t, roi2)
    raw_ax.plot(t[idx_peaks], roi2[idx_peaks], 'ro')
    raw_ax.set_xlabel("s")
    raw_ax.set_ylabel("roi2")

    # funcs to get plateau region
    plateau = lambda idx, r: r[idx-plateau_start_idx:idx-plateau_end_idx]
    plateau_t = lambda idx: t[idx-plateau_start_idx:idx-plateau_end_idx]
    
    # change in plateaus over time
    plateaus_roi1 = [plateau(i, roi1).mean() for i in idx_peaks] # plateaus from roi 2 (FRAP zone)
    plateaus_roi2 = [plateau(i, roi2).mean() for i in idx_peaks] # plateaus from roi 2 (FRAP zone)
    plateaus_roi3 = [plateau(i, roi3).mean() for i in idx_peaks] # plateaus from roi 3 (non-FRAP zone in same cell)
    # plateaus_roi4 = [plateau(i, roi4).mean() for i in idx_peaks] # plateaus from roi 2 (FRAP zone)
    
    # plot stim-triggered averages
    fig_stimavg, ax = plt.subplots(1,4, sharex='all', sharey='all')
    fig.set_figwidth(8)
    for i,roi in enumerate(all_rois):
    
        
        t_combo = np.arange(-1*samples_pre_stim, length_to_plot)/s_rate
        # print(t_combo)
        roi_avg = np.zeros_like(t_combo)
        for j in range(num_stims):
            start_idx = idx_peaks[j]
            current_roi = roi[start_idx + np.arange(-1*samples_pre_stim,length_to_plot)]
            
            if i == 0:
                current_roi = current_roi - plateaus_roi1[j]
            elif i == 1:
                current_roi = current_roi - plateaus_roi2[j]
            elif i == 2:
                current_roi = current_roi - plateaus_roi3[j]
#                elif i == 3:
#                    current_roi = current_roi - plateaus_roi4[j]
            roi_avg += current_roi
            ax[i].plot(t_combo, current_roi, 'gray')
        
        roi_avg = roi_avg/(j+1)
        ax[i].plot(t_combo, np.zeros_like(t_combo), 'k--')
        ax[i].plot(t_combo, roi_avg, 'r-')
        ax[i].set_title('roi ' + str(i+1))
        
        if i == 1: # roi 2 only (FRAP point)
            if construct not in all_roi_avg_data.keys():
                all_roi_avg_data[construct] = []
            all_roi_avg_data[construct].append(roi_avg)
        
        if normalize_roi:
            ax[i].set_ylim([-1.5, 1])
        else:
            ax[i].set_ylim([-1000, 500])
    
    # plot plateaus on entire timeseries and save (roi2 only)
    for i in idx_peaks:
        current_plateau = plateau(i, roi2)
        raw_ax.plot(plateau_t(i), current_plateau, 'k-')
            
    # percent change in plateaus (next plateau - current plateau ) / current plateau
    percent_change = np.diff(plateaus_roi2) / plateaus_roi2[0:-1]
    
    if construct not in plateau_data:
        plateau_data[construct] = []
    plateau_data[construct].append(percent_change)
    
    
    # temporarily commenting out plotting stuff
    '''
    fig_plateaus, ax = plt.subplots()
    ax.plot(plateaus_roi2, 'ko-')
    ax.plot(plateaus_roi3, 'ro-')
    ax.legend(['FRAP zone (roi2)', 'non-FRAP zone (roi3)'])
    ax.set_title("plateaus")
    ax.set_xlabel("stim num")
    ax.set_ylabel("frap roi")
    
    
    fig_prcnt_chng = plt.figure()
    plt.plot(percent_change, 'ko-')
    plt.title("fraction change plateau")
    plt.xlabel("stim num")
    plt.ylabel("fraction change")
    
    
    
    if save_figs:
        # save plots
        fig_stimavg.savefig(os.path.join('./analysis', folder_name, filename.replace('.csv', '.png')))
        fig_prcnt_chng.savefig(os.path.join('./analysis', folder_name, 'percent_change_' + filename.replace('.csv', '.png')))
        fig_plateaus.savefig(os.path.join('./analysis', folder_name, 'plateaus_' + filename.replace('.csv', '.png')))
    '''
  
    if not keep_figs_open:
        plt.close('all')    


if save_data:
    
    plateau_name = 'plateau_data_{}_{}_{}'.format('norm' if normalize_roi else '', solution_condition, stim_laser)
    with open(r'./analysis/'+ plateau_name + '.pkl', 'wb') as f:
        pickle.dump(plateau_data, f, pickle.HIGHEST_PROTOCOL)
    
    roi_traces_name = 'roi_traces_{}_{}_{}'.format('norm' if normalize_roi else '', solution_condition, stim_laser)
    with open(r'./analysis/'+ roi_traces_name + '.pkl', 'wb') as f:
        pickle.dump(all_roi_avg_data, f, pickle.HIGHEST_PROTOCOL)
    
    # convert plateau to dataframe, save
    '''
    df = pd.DataFrame.from_dict(plateau_data)
    df.to_csv(r'./analysis/plateau_data{}.csv'.format('_norm' if normalize_roi else ''))
    df.to_pickle(r'./analysis/plateau_data{}.pkl'.format('_norm' if normalize_roi else ''))
    
    # df_roi = pd.DataFrame.from_dict(all_roi_avg_data)
    
    # save all_roi_avg_data
    if normalize_roi:
        df_roi.to_pickle(r'./analysis/all_roi_avg_data_norm.pkl')
    else:
        df_roi.to_pickle(r'./analysis/all_roi_avg_data.pkl')
    '''