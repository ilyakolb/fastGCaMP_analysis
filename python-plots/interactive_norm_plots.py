# -*- coding: utf-8 -*-
'''
load and plot data generated by normPlots.m MATLAB file using plotly

normPlots.mat: contains normPlots_struct and nAP

normPlots_struct:
dtype([('construct', 'O'), ('dff_mean', 'O'), ('dff_sterr', 'O'), ('SNR_mean', 'O'), ('SNR_sterr', 'O'), ('halfrise_mean', 'O'), ('halfrise_sterr', 'O'), ('timetopeak_mean', 'O'), ('timetopeak_sterr', 'O'), ('halfdecay_mean', 'O'), ('halfdecay_sterr', 'O')])

total 5 plots to generate

@todo:
    multi row subplots
    
'''
from scipy.io import loadmat
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
import plotly.express as px

plot_subset_for_paper = 0 # 0: plot entire dataset, save as html. 1: plot subset, save as pdf for paper
html_write_dir = r'D:\site\ilyakolb.github.io\interactive_norm_plots.html'
pdf_dir = r"D:\ufgcamp_paper_data\culture-screen-figs/norm_plots.pdf"
plot_mat = loadmat(r'data/plotly_normPlots_WEBSITE.mat')

n_subplots = 6
n_rows = 2
n_cols = 3

pio.templates.default = "plotly_white"

subplot_titles = ['peak dF/F', 'SNR (norm.)', 'half-rise time (norm.)', 'full rise time (norm.)', 'half-decay time (norm.)']
fig = make_subplots(rows=n_rows, cols=n_cols, subplot_titles=subplot_titles, x_title='number of action potentials') 


colorscheme = px.colors.qualitative.Alphabet # 0th: control

nAPs = plot_mat['nAPs'][0]

# labels and order of legend
if plot_subset_for_paper:
    hits_label = ['jGCaMP8f', 'jGCaMP8m', 'jGCaMP8s', 'GCaMP6s', 'jGCaMP7f', 'jGCaMP7s', 'XCaMP-Gf']# [h[0] for h in hits[0]]
else:
    hits_label = ['jGCaMP8f', 'jGCaMP8m', 'jGCaMP8s','jGCaMP8.712', 'GCaMP6s', 'GCaMP6f', 'jGCaMP7f', 'jGCaMP7s', 
                'jGCaMP7c', 'jGCaMP7b', 'XCaMP-Gf', 'XCaMP-G', 'XCaMP-Gf0']# [h[0] for h in hits[0]]

all_norm_plots = plot_mat['normPlots_struct'][0]

all_norm_plots_constructs = [c[0] for c in all_norm_plots['construct']]

for i,construct_name in enumerate(hits_label): # enumerate(all_norm_plots):
    # construct_name = c['construct'][0]
    c = all_norm_plots[all_norm_plots_constructs.index(construct_name)]
    
    for j in range(1,n_subplots):
        x_plot = nAPs
        y_plot = c[2*j-1].squeeze()
        
        # if decay plot, remove the 160AP data (unreliable)
        if all_norm_plots.dtype.names[2*j-1] == 'halfdecay_mean_norm':
            x_plot = x_plot[:-1]
            y_plot = y_plot[:-1]
        
        fig.add_trace(go.Scatter(
            x=x_plot,
            y=y_plot,
            error_y=dict(
                type='data', # value of error bar given in data coordinates
                array=c[2*j].squeeze(),
                visible=True),
            showlegend = (j == 1),
            legendgroup = construct_name,
            line_color = colorscheme[i],
            name = construct_name,
        ), row=(1 if j < 4 else 2),
        col= (1 if (j == 1 or j == 4) else (2 if (j == 2 or j == 5) else 3))
        )
  
fig.update_traces(mode='lines')
fig.update_layout(hovermode="closest", #, width=800, height=400,
    title="<b>jGCaMP8 performance in cultured neurons</b>",
    title_font_size=18,
    #xaxis_title="time (s)",
    legend_title="variants",
    font=dict(
        family="Arial",
        size=14,
    )
)
fig.update_xaxes(type="log",
        tickmode = 'array',
        tickvals = nAPs
    )

fig.show()

if plot_subset_for_paper:
    fig.write_image(pdf_dir)
else:
    fig.write_html(html_write_dir, auto_open=True)



