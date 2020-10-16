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

def add_shaded_trace(t, y_mean, y_sterr, stim_subplot_idx, color, construct_label):
    '''
    plot AP mean trace +/- s.e.m

    Parameters
    ----------
    t : 1d numpy array
        time array (s).
    y_mean : 1d numpy array
        mean response.
    y_sterr : 1d numpy array
        st.err around mean of response.
    stim_subplot_idx : int
        index of subplot (e.g. 1 for 1AP, 4 for 160AP).
    color : str
        color of trace (e.g. '#FEAF16').
    variant_label : str
        construct name for legend.

    Returns
    -------
    None.

    '''
    y_upper = y_mean + y_sterr
    y_lower = y_mean- y_sterr
    y_lower = np.flip(y_lower)
    t_rev = np.flip(t)
    
    # draw shaded area
    fig.add_trace(go.Scatter(
    x=np.append(t, t_rev),
    y=np.append(y_upper, y_lower),
    fill='toself',
    fillcolor= color, # 'rgba(0,100,80,0.2)',
    line_color='rgba(255,255,255,0)',
    showlegend=False,
    name=construct_label,
    opacity= 0.5
    ), row=1,
    col=stim_subplot_idx)
    
    
    # draw mean line
    fig.add_trace(go.Scatter(
    x=t, y=y_mean,
    line_color=color, # 'rgb(0,100,80)',
    name=construct_label,
    showlegend=(stim_subplot_idx == 1),
    ), row=1,
    col=stim_subplot_idx,)

plot_mat = loadmat(r'..\plotting.mat')

control_med_med_dff = plot_mat['plot_out'][0,0]['control_med_med_dff'] # [time x nStims]
hits_med_med_dff = plot_mat['plot_out'][0,0]['hits_med_med_dff'] # [time x nStims x nHits]
time = plot_mat['plot_out'][0,0]['time'][0] # time vector
hits = plot_mat['plot_out'][0,0]['hits'] # array of hit name strings
control = plot_mat['plot_out'][0,0]['control'][0] # control name string
hits_med_med_dff_sterr = np.random.random_sample(hits_med_med_dff.shape)/100
control_med_med_dff_sterr = np.random.random_sample(control_med_med_dff.shape)/100

hits_label = [h[0] for h in hits[0]]
n_stims = control_med_med_dff.shape[1]
n_hits = hits_med_med_dff.shape[2]

stim_names = {'1', '3', '10', '160'}
colorscheme = px.colors.qualitative.Alphabet # 0th: control

fig = make_subplots(rows=1, cols=n_stims) # go.Figure()

# cycle through stim number
for i in range(n_stims):
    
    # plot control
    add_shaded_trace(time, control_med_med_dff[:,i], control_med_med_dff_sterr[:,i], i+1, colorscheme[0], control)
    
    # plot hits
    for j in range(n_hits):
        add_shaded_trace(time, hits_med_med_dff[:,i,j], hits_med_med_dff_sterr[:,i,j], i+1, colorscheme[j+1], hits_label[j])

fig.update_traces(mode='lines')
fig.update_layout(hovermode="closest", width=1400, height=800)
fig.show()
fig.write_html('first_figure.html', auto_open=True)

