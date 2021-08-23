load('F:\ufgcamp_paper_data\sample_culture_movies\jG8m\96Well28-C04\para_array_cherry.mat')

savedir = 'C:\Users\labadmin\Dropbox (HHMI)\jGCaMP8_paper_for_coauthors\neuron_culture_figures';
figure

f = para_array(1,36).df_f;
t = (1:length(f))/200;
plot(t, f, 'k-')
% savefig([savedir filesep 'sample_jG8_well_trace.fig'])
exportgraphics(gcf,[savedir filesep 'sample_jG8_well_trace.pdf'],'ContentType','vector')
