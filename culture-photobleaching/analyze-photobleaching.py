# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 19:03:47 2020

@author: kolbi

photobleaching of gcamp variants
exp performed by Daniel Reep

images taken every 10 s
imaging power same as during screening

protocol:
    background subtraction, normalize to first point
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')
T_s = 10 # sampling period (s)

def get_column_name(table, plate_id, construct_name, cell_or_bg_num):
    '''
    Parameters
    ----------
    table : pandas dataframe
        table of photobleaching results
    plate_id : str
        e.g. P3a or P4a.
    construct_name : str
        e.g. '500.456'.
    cell_or_bg_num : int
        0 for background. 1-3 for cell num

    Returns
    -------
    Corresponding pandas column for data table. If not found or > 1 entry found, error

    '''
    cols = data.columns
    
    if cell_or_bg_num == 0:
        # background trace
        col_name = [m for m in cols if (construct_name in m and 'background' in m and plate_id in m)]
    else:
        col_name = [m for m in cols if (construct_name in m and 'cell'+str(cell_or_bg_num) in m and plate_id in m)]
    
    # check to make sure column is unique
    if len(col_name) > 1:
        print(col_name)
        raise Exception('More than 1 column found that matches plate, construct name, and cell num!')
    elif len(col_name) < 1:
        print(col_name)
        raise Exception('No column found matching plate, construct name, and cell num')
    return col_name[0]
    
data = pd.read_csv('photobleaching-DR.csv')
hits = ['500dot456', '500dot686', '500dot688', '500dot712', '500dot543', '500dot707', '500dot640', '500dot455', '10dot921', '10dot1473', '10dot1513', '10dot1561', '538dot1']

plate_names = ['P3a', 'P4a']

t = np.arange(0,len(data))*T_s

df_plt = pd.DataFrame(index = hits, columns=['bleach'])
for h in hits:
    # plt.figure()
    
    h_array = list()
    for i,plate in enumerate(plate_names):
        bg_trace = data[get_column_name(data, plate, h, 0)]
        
        # ax = plt.subplot(2,1,i+1)
        for c in range(1,4):
            cell_trace = data[get_column_name(data, plate, h, c)]
            cell_trace_bg = cell_trace - bg_trace
            cell_trace_bg = cell_trace_bg / cell_trace_bg[0]
            # plt.plot(cell_trace_bg)
            
            # add array to list
            h_array.append(cell_trace_bg.to_numpy())
        # ax.set_xlabel('Time (s)')
        # ax.set_title(h)
    
    df_plt.loc[h] = [h_array]

# remove bad cells
df_plt.loc['10dot1473']['bleach'].pop(0) # cell is all over the place
df_plt.loc['10dot1513']['bleach'].pop(4) # cell fluorescence goes up

plt.figure()

i=1
# plot mean +/- std for each sensor
for index, row in df_plt.iterrows():
    
    ax = plt.subplot(4,4,i)
    i+=1
    bleach_nd = np.array(row['bleach']) # convert to 2d np array
    bleach_mean = np.mean(bleach_nd, axis = 0)
    bleach_std = np.std(bleach_nd, axis = 0)
    # _,ax= plt.subplots()
    plt.fill_between(t, bleach_mean+bleach_std, bleach_mean-bleach_std, 
                     facecolor="blue", # The fill color
                     color='blue',       # The outline color
                     alpha=0.2)          # Transparency of the fill)
    ax.set_title(index)
    # ax.set_xlabel('Time (s)')
    # plt.tight_layout()
    ax.set_ylim([0, 1.2])
    
plt.tight_layout()