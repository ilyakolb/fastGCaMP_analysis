# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 09:19:08 2020

@author: kolbi

@todo:
    why is response occuring before 1 s?
    
"""

import numpy as np
from scipy.io import loadmat
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
import plotly.express as px


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
    legendgroup=construct_label,
    name=construct_label,
    opacity= 0.5
    ), row=1,
    col=stim_subplot_idx)
    
    
    # draw mean line
    fig.add_trace(go.Scatter(
    x=t, y=y_mean,
    line_color=color, # 'rgb(0,100,80)',
    name=construct_label,
    legendgroup=construct_label,
    showlegend=(stim_subplot_idx == 1),
    ), row=1,
    col=stim_subplot_idx,)

t_to_ignore_s = 0.5 # initial seconds to remove to get rid of bleaching artifacts
pio.templates.default = "plotly_white"
plot_mat = loadmat(r'..\plotting.mat')

control_med_med_dff = plot_mat['plot_out'][0,0]['control_med_med_dff'] # [time x nStims]
hits_med_med_dff = plot_mat['plot_out'][0,0]['hits_med_med_dff'] # [time x nStims x nHits]
time = plot_mat['plot_out'][0,0]['time'][0] # time vector
hits = plot_mat['plot_out'][0,0]['hits'] # array of hit name strings
control = plot_mat['plot_out'][0,0]['control'][0] # control name string
hits_med_med_dff_sterr =  plot_mat['plot_out'][0,0]['hits_med_med_dff_sterr']
control_med_med_dff_sterr = plot_mat['plot_out'][0,0]['control_med_med_dff_sterr']

fs = 1/(time[1]-time[0])
t_to_ignore_samples = int(t_to_ignore_s * fs)

# ignore samples in beginning
time = time[t_to_ignore_samples:]
control_med_med_dff = control_med_med_dff[t_to_ignore_samples:,]
hits_med_med_dff = hits_med_med_dff[t_to_ignore_samples:,]
control_med_med_dff_sterr = control_med_med_dff_sterr[t_to_ignore_samples:,]
hits_med_med_dff_sterr = hits_med_med_dff_sterr[t_to_ignore_samples:,]

hits_label = [h[0] for h in hits[0]]
n_stims = control_med_med_dff.shape[1]
n_hits = hits_med_med_dff.shape[2]

stim_names = ('1 AP', '3 AP', '10 AP', '160 AP')
colorscheme = px.colors.qualitative.Alphabet # 0th: control

fig = make_subplots(rows=1, cols=n_stims, subplot_titles=stim_names, x_title='Time (s)',   y_title='dF/F') # go.Figure()

# cycle through stim number
for i in range(n_stims):
    
    # plot control
    add_shaded_trace(time, control_med_med_dff[:,i], control_med_med_dff_sterr[:,i], i+1, colorscheme[0], control)
    
    # plot hits
    for j in range(n_hits):
        add_shaded_trace(time, hits_med_med_dff[:,i,j], hits_med_med_dff_sterr[:,i,j], i+1, colorscheme[j+1], hits_label[j])

fig.update_traces(mode='lines')
fig.update_layout(hovermode="closest", #, width=800, height=400,
    title="AP responses in cultured neuron screen",
    title_font_size=18,
    #xaxis_title="time (s)",
    legend_title="variants",
    font=dict(
        family="Arial",
        size=14,
    )
)
fig.show()
fig.write_html('interactive_AP_traces.html', auto_open=True)

