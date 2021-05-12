# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 12:59:13 2020

@author: kolbi

NOTES:
    combo analysis of FRAP experiments
    
USAGE:
    env: conda activate base2
    run from D:\pythonTesting\jg8-frap
    Obtain data for plot_average_2.py and regular_vs_iono_treatments.py:
        
        --- all 405-stimmed constructs
        bleachlaser_condition = 'stim405'
        solution_condition = 'regular'
        all_constructs =   ['604.2', '10.641', '500.688','500.686']
        
        --- all 488-stimmed constructs
        bleachlaser_condition = 'stim405'
        solution_condition = 'regular'
        all_constructs =   ['604.2', '10.641']
               
        --- all 405-stimmed iono constructs
        bleachlaser_condition = 'stim405' # 'stim405' or 'stim488'
        solution_condition = 'iono' # 'regular' or 'iono'
        all_constructs =   ['604.2', '10.641', '500.688','500.686']

    
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

def rescale_0to1(_trace, _max):
    # rescale so that _max value (plateau) = 1 and min value of trace = 0
    _t = _trace * 1/(_max-np.min(_trace))
    return _t - np.min(_t)

def rescale_to1(_trace, _max):
    # rescale so that _max value (plateau) = 1
    _t = _trace / _max
    return _t

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

save_figs       = True
save_data       = True
normalize_roi   = True
keep_figs_open  = False # True to keep all generated figures open. Memory errors if too many open

bleachlaser_condition = 'stim405' # 'stim405' or 'stim488'
solution_condition = 'regular' # 'regular' or 'iono'

# ['604.2', '10.641'] # ['604.2', '10.641', '500.688','500.686']
all_constructs = ['EGFP.B-actin' , 'mEm.Cyto', '604.2', '10.641', '500.688','500.686']

length_to_plot_s = 6.5
seconds_pre_stim = 0.4
plateau_start_sec = 2
plateau_end_sec = 1.5

num_peaks_to_plot = 5 # or 10
peak_thresh = 1 # diff > peak_thresh used to detect stimuli
# s_rate = 25# was 50

# directory of combined data
combo_dir = r'Z:\ilya\code\fastGCaMP_analysis\jg8-frap\data\combined' # r'Z:\ilya\code\fastGCaMP_analysis\jg8-frap\data\exp7_20210511' #  #
csv_filenames_all = os.listdir(combo_dir)

findmatches = lambda x: any([1 for c in all_constructs if c in x]) and solution_condition in x and bleachlaser_condition in x and (x.endswith('.csv'))

csv_filenames = list(filter(findmatches, csv_filenames_all))

print('Evaluating')
[print(p) for p in csv_filenames]

plateau_data = {} # dict for storing all plateau data
all_roi_avg_data = {} # 

# for construct in all_constructs:

for csv_filename in csv_filenames:

    (construct, expdate, cellnum, imgbuffer, stim_laser) = decode_filename(csv_filename)    
    # plateau_data[construct] = []
    # all_roi_avg_data[construct] = []
    data = pd.read_csv(os.path.join(combo_dir, csv_filename)) 
    t = data['Time [s]']
    
    # get constants
    s_rate = int(1/(t[1]-t[0]))
    length_to_plot = int(length_to_plot_s*s_rate) # was 325
    samples_pre_stim = int(seconds_pre_stim*s_rate if bleachlaser_condition == 'stim405' else 3/2*seconds_pre_stim*s_rate) # was 20
    plateau_start_idx = int(plateau_start_sec*s_rate) # was 100
    plateau_end_idx = int(plateau_end_sec*s_rate)
    
    # exp7_20210511 data does not have #1 channel 
    if '#1 (CSU (488))' not in data.keys():
        roi1=np.zeros_like(t)
    else:
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
    roi2_diff = np.diff(roi2, append=roi2[-1])
    idx_peaks, _ = find_peaks(roi2_diff, height=peak_thresh, distance=20)
    
    # fudge factor to adjust for stim duration
    idx_peaks = idx_peaks + 2
    print('total number of stims: {}'.format(len(idx_peaks)))
    if num_peaks_to_plot != 'all':

        idx_peaks = idx_peaks[:num_peaks_to_plot]
            
    # adjust peak indices to find rightmost peak
    for i,_ in enumerate(idx_peaks):
        # find rightmost edge of stimulus to align to
        idx_peaks[i] += np.where(roi2[idx_peaks[i]:idx_peaks[i] +10] > roi2[idx_peaks[i]] / 2)[0][-1]
    
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
    fig_stimavg, ax = plt.subplots(1,1)
    fig.set_figwidth(8)
    percent_change = np.zeros(num_stims)
    
    # for i,roi in enumerate(all_rois):
    
        
    t_combo = np.arange(-1*samples_pre_stim, length_to_plot)/s_rate
    # print(t_combo)
    roi_avg = np.zeros_like(t_combo)
    
    for j in range(num_stims):
        start_idx = idx_peaks[j]
        current_roi = roi2[start_idx + np.arange(-1*samples_pre_stim,length_to_plot)]
        
        current_roi = rescale_to1(current_roi, plateaus_roi2[j])
        percent_change[j] = np.mean(current_roi[-1*plateau_end_idx:])
        roi_avg += current_roi

    roi_avg = roi_avg/(j+1)
    
    ax.plot(t_combo, np.ones_like(t_combo), 'k--')
    ax.plot(t_combo, roi_avg, 'r-')
    ax.set_title('roi2')
    ax.set_ylim([-0.5, 1.5])
    
    if construct not in all_roi_avg_data.keys():
        all_roi_avg_data[construct] = []
    
    '''
    a subset of EGFP.B-actin recordings are recorded at 25 Hz (~1/2 of normal 50 Hz)
    for those, interpolate to match 50 Hz rate
    '''
    if roi_avg.size == 172:
        roi_avg = np.interp(np.arange(0,345), np.arange(0,344,2), roi_avg)
        print('Encountered low-sampling rate recording: upsampling...')
    all_roi_avg_data[construct].append(roi_avg)

    # percent_change[j] = (np.mean(current_roi[-1*plateau_end_idx:]) - plateaus_roi2[j])/plateaus_roi2[j]
    print(percent_change)

    # plot plateaus on entire timeseries and save (roi2 only)
    for idx in idx_peaks:
        current_plateau = plateau(idx, roi2)
        raw_ax.plot(plateau_t(idx), current_plateau, 'k-')
            
    # percent change in plateaus (next plateau - current plateau ) / current plateau
    # percent_change = np.diff(plateaus_roi2) / plateaus_roi2[0:-1]
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
    '''
    
    
    if save_figs:
        # save plots
        fig_stimavg.savefig(os.path.join('./analysis', folder_name, csv_filename.replace('.csv', '.png')))
        fig_stimavg.savefig(os.path.join('./analysis', folder_name, csv_filename.replace('.csv', '.pdf')))
        
        # fig_prcnt_chng.savefig(os.path.join('./analysis', folder_name, 'percent_change_' + csv_filename.replace('.csv', '.png')))
        # fig_plateaus.savefig(os.path.join('./analysis', folder_name, 'plateaus_' + csv_filename.replace('.csv', '.png')))
    
  
    if not keep_figs_open:
        plt.close('all')    


if save_data:
    
    plateau_name = 'plateau_data_{}_{}_{}_{}peak'.format('norm' if normalize_roi else '', solution_condition, stim_laser, num_peaks_to_plot)
    with open(r'./analysis/'+ plateau_name + '.pkl', 'wb') as f:
        pickle.dump(plateau_data, f, pickle.HIGHEST_PROTOCOL)
    
    roi_traces_name = 'roi_traces_{}_{}_{}_{}peak'.format('norm' if normalize_roi else '', solution_condition, stim_laser, num_peaks_to_plot)
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