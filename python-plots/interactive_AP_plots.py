# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 09:19:08 2020

@author: kolbi
"""

import numpy as np
from scipy.io import loadmat
import plotly.graph_objects as go
import plotly.io as pio

pio.templates.default = "plotly_white"

def add_shaded_trace(t, y_mean, y_sterr):
    y_upper = y_mean + y_sterr
    y_lower = y_mean- y_sterr
    y_lower = np.flip(y_lower)
    t_rev = np.flip(t)
    
    # add shade
    fig.add_trace(go.Scatter(
    x=np.append(t, t_rev),
    y=np.append(y_upper, y_lower),
    fill='toself',
    fillcolor='rgba(0,100,80,0.2)',
    line_color='rgba(255,255,255,0)',
    showlegend=False,
    name='GCaMP6s',
    ))
    
    # add line
    fig.add_trace(go.Scatter(
    x=t, y=y_mean,
    line_color='rgb(0,100,80)',
    name='GCaMP6s',
    ))

plot_mat = loadmat(r'..\plotting.mat')

control_med_med_dff = plot_mat['plot_out'][0,0]['control_med_med_dff'] # [time x nStims]
hits_med_med_dff = plot_mat['plot_out'][0,0]['hits_med_med_dff'] # [time x nStims x nHits]
time = plot_mat['plot_out'][0,0]['time'][0] # time vector
hits = plot_mat['plot_out'][0,0]['hits'] # array of hit name strings
control = plot_mat['plot_out'][0,0]['control'] # control name string
hits_med_med_dff_sterr = np.random.random_sample(hits_med_med_dff.shape)/10
control_med_med_dff_sterr = np.random.random_sample(control_med_med_dff.shape)/10

fig = go.Figure()

add_shaded_trace(time, control_med_med_dff[:,0], control_med_med_dff_sterr[:,0])

fig.update_traces(mode='lines')
fig.update_layout(hovermode="closest", width=1400, height=800)
fig.show()
fig.write_html('first_figure.html', auto_open=True)

