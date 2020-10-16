# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 09:19:08 2020

@author: kolbi
"""

import numpy as np
from scipy.io import loadmat
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
import plotly.express as px

pio.templates.default = "plotly_white"

def add_shaded_trace(t, y_mean, y_sterr, stim_subplot_idx, color):
    '''
    plot AP mean trace +/- s.e.m

    Parameters
    ----------
    t : 1d numpy array
        time array (s).
    y_mean : TYPE
        DESCRIPTION.
    y_sterr : TYPE
        DESCRIPTION.
    stim_subplot_idx : TYPE
        DESCRIPTION.
    color : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    y_upper = y_mean + y_sterr
    y_lower = y_mean- y_sterr
    y_lower = np.flip(y_lower)
    t_rev = np.flip(t)
    
    # add shade
    fig.add_trace(go.Scatter(
    x=np.append(t, t_rev),
    y=np.append(y_upper, y_lower),
    fill='toself',
    fillcolor= color, # 'rgba(0,100,80,0.2)',
    line_color='rgba(255,255,255,0)',
    showlegend=False,
    name='GCaMP6s',
    opacity= 0.5
    ), row=1,
    col=stim_subplot_idx)
    
    # add line
    fig.add_trace(go.Scatter(
    x=t, y=y_mean,
    line_color=color, # 'rgb(0,100,80)',
    name='GCaMP6s'
    ), row=1,
    col=stim_subplot_idx,)

plot_mat = loadmat(r'..\plotting.mat')

control_med_med_dff = plot_mat['plot_out'][0,0]['control_med_med_dff'] # [time x nStims]
hits_med_med_dff = plot_mat['plot_out'][0,0]['hits_med_med_dff'] # [time x nStims x nHits]
time = plot_mat['plot_out'][0,0]['time'][0] # time vector
hits = plot_mat['plot_out'][0,0]['hits'] # array of hit name strings
control = plot_mat['plot_out'][0,0]['control'] # control name string
hits_med_med_dff_sterr = np.random.random_sample(hits_med_med_dff.shape)/10
control_med_med_dff_sterr = np.random.random_sample(control_med_med_dff.shape)/10

n_stims = control_med_med_dff.shape[1]
n_hits = hits_med_med_dff.shape[2]

stim_names = {'1', '3', '10', '160'}
colorscheme = px.colors.qualitative.Plotly # 0th: control

fig = make_subplots(rows=1, cols=n_stims) # go.Figure()

# cycle through stim number
for i in range(n_stims):
    
    # plot control
    add_shaded_trace(time, control_med_med_dff[:,i], control_med_med_dff_sterr[:,i], i+1, colorscheme[0])
    
    # plot hits
    for j in range(n_hits):
        pass

fig.update_traces(mode='lines')
fig.update_layout(hovermode="closest", width=1400, height=800)
fig.show()
fig.write_html('first_figure.html', auto_open=True)

