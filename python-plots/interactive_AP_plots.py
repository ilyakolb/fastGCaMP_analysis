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


plot_subset_for_paper = True # 0: plot entire dataset, save as html. 1: plot subset, save as pdf for paper
html_write_dir = r'./figs/interactive_AP_traces.html'
pdf_dir = r"./figs/AP_plots.pdf"

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
    ), row=(2 if stim_subplot_idx > 2 else 1),
    col=(1 if stim_subplot_idx%2 == 1 else 2))
    
    
    # draw mean line
    fig.add_trace(go.Scatter(
    x=t, y=y_mean,
    line_color=color, # 'rgb(0,100,80)',
    name=construct_label,
    legendgroup=construct_label,
    showlegend=(stim_subplot_idx == 1),
    ), row=(2 if stim_subplot_idx > 2 else 1),
    col=(1 if stim_subplot_idx%2 == 1 else 2),)

t_to_ignore_s = 0.5 # initial seconds to remove to get rid of bleaching artifacts
pio.templates.default = "plotly_white"
plot_mat = loadmat(r'data/plotly_AP_traces_WEBSITE.mat')

control_med_med_dff = plot_mat['plot_out'][0,0]['control_med_med_dff'] # [time x nStims]
hits_med_med_dff = plot_mat['plot_out'][0,0]['hits_med_med_dff'] # [time x nStims x nHits]
time = plot_mat['plot_out'][0,0]['time'][0] # time vector
hits = plot_mat['plot_out'][0,0]['hits'] # array of hit name strings
hits = [h[0] for h in hits[0]]


control = plot_mat['plot_out'][0,0]['control'][0] # control name string
hits_med_med_dff_sterr =  plot_mat['plot_out'][0,0]['hits_med_med_dff_sterr']
control_med_med_dff_sterr = plot_mat['plot_out'][0,0]['control_med_med_dff_sterr']

fs = 1/(time[1]-time[0])
t_to_ignore_samples = int(t_to_ignore_s * fs)

# ignore samples in beginning and end (due to alignment)
time = time[t_to_ignore_samples:-1*t_to_ignore_samples]
control_med_med_dff = control_med_med_dff[t_to_ignore_samples:-1*t_to_ignore_samples,]
hits_med_med_dff = hits_med_med_dff[t_to_ignore_samples:-1*t_to_ignore_samples,]
control_med_med_dff_sterr = control_med_med_dff_sterr[t_to_ignore_samples:-1*t_to_ignore_samples,]
hits_med_med_dff_sterr = hits_med_med_dff_sterr[t_to_ignore_samples:-1*t_to_ignore_samples,]

# mapping between construct IDs and names
mapping = {'GCaMP6s': '10.641' , 'jGCaMP8f': '500.456', 'jGCaMP8m': '500.686', 'jGCaMP8s': '500.688','jGCaMP8.712': '500.712', 'GCaMP6f': '10.693', 'jGCaMP7f': '10.921', 'jGCaMP7s': '10.1473', 'jGCaMP7c': '10.1513', 'jGCaMP7b': '10.1561', 'XCaMP-Gf': '538.1', 'XCaMP-G': '538.2', 'XCaMP-Gf0': '538.3'} # [h[0] for h in hits[0]]

control_label= 'GCaMP6s'
# hits_label = ['GCaMP6f', 'jGCaMP7f', 'jGCaMP8f', 'jGCaMP8m', 'jGCaMP8s', 'jGCaMP8.712', 
#              'jGCaMP8.543', 'jGCaMP8.707', 'jGCaMP8.455', 'jGCaMP7s', 'jGCaMP7c', 'jGCaMP7b', 'XCaMP-Gf', 'XCaMP-G', 'XCaMP-Gf0']# [h[0] for h in hits[0]]

# labels and order of legend
if plot_subset_for_paper:
    hits_label = ['jGCaMP8f', 'jGCaMP8m', 'jGCaMP8s', 'GCaMP6s', 'jGCaMP7f', 'jGCaMP7s', 'XCaMP-Gf']# [h[0] for h in hits[0]]
else:
    hits_label = ['jGCaMP8f', 'jGCaMP8m', 'jGCaMP8s','jGCaMP8.712', 'GCaMP6s', 'GCaMP6f', 'jGCaMP7f', 'jGCaMP7s', 
                'jGCaMP7c', 'jGCaMP7b', 'XCaMP-Gf', 'XCaMP-G', 'XCaMP-Gf0']# [h[0] for h in hits[0]]

n_stims = control_med_med_dff.shape[1]
n_hits = hits_med_med_dff.shape[2]

stim_names = ('1 AP', '3 AP', '1 AP zoomed', '3 AP zoomed')

if plot_subset_for_paper:
    colorscheme = ['#336699', '#CC3333', '#666666', '#FFCC99', '#66CC33', '#CC99CC', '#3399CC']
else:
    colorscheme = px.colors.qualitative.Alphabet # 0th: control

stim_iter = [0, 1, 0, 1] # [1AP, 3AP, 1AP, 3AP]

fig = make_subplots(rows=2, cols=2, subplot_titles=stim_names, x_title='Time (s)',   y_title='dF/F') # go.Figure()

# cycle through stim number
for plot_i, i in enumerate(stim_iter):

    # plot hits
    for j,hit in enumerate(hits_label): # for j in range(n_hits):
        
        if hit == control_label:
            # plot control
            add_shaded_trace(time, control_med_med_dff[:,i], control_med_med_dff_sterr[:,i], plot_i+1, colorscheme[j], control_label)
        else:
            hit_id = mapping.get(hit)
            hit_idx = hits.index(hit_id)
            add_shaded_trace(time, hits_med_med_dff[:,i,hit_idx], hits_med_med_dff_sterr[:,i,hit_idx], plot_i+1, colorscheme[j], hit)

fig.update_traces(mode='lines')
fig.update_layout(hovermode="closest", #, width=800, height=400,
    title="<b>AP responses in cultured neuron screen</b>",
    title_font_size=18,
    #xaxis_title="time (s)",
    legend_title="variants",
    font=dict(
        family="Arial",
        size=14,
    ),
    height=300
)


fig.update_yaxes(showline=True, linewidth=1, linecolor='black', ticks='inside', showgrid=False)

# for 1AP, 3AP case, zoom in
fig.update_xaxes(range=[0.8, 1.6], row=1, col=1, ticks='inside', showline=True, linewidth=1, linecolor='black', nticks=6, showgrid=False)
fig.update_xaxes(range=[0.8, 1.6], row=1, col=2, ticks='inside', showline=True, linewidth=1, linecolor='black', nticks=6, showgrid=False)

# for bottom row, zoom in more
fig.update_xaxes(range=[0.95, 1.2], row=2, col=1, ticks='inside', showline=True, linewidth=1, linecolor='black', nticks=6, showgrid=False)
fig.update_xaxes(range=[0.95, 1.2], row=2, col=2, ticks='inside', showline=True, linewidth=1, linecolor='black', nticks=6, showgrid=False)

fig.show()

# save subset version to pdf directory, save webpage version to site directory
if plot_subset_for_paper:
    fig.write_image(pdf_dir)
else:
    fig.write_html(html_write_dir, auto_open=True)
