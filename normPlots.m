%% Run this after running compare_constructs_GCaMP96uf
clc
close all

% set to 1 to save figures
saveOn = 0;

saveFolder = 'D:\Dropbox (HHMI)\janelia\writing\jGCaMP8 patent\figures\all\';
% nAP vs parameter plots (like jGCaMP7 paper)
nAPs = [1 3 10 160];
relevantAPs = [1 2 3 4];

% filter only mutants that are in 'hits'
mutant_hits = controlMutant;
for i = 1:length(hits)
    mutant_hits = [mutant_hits mutant(contains({mutant.construct}, hits{i}))];
end
% addendum: string to add on to construct names for clarity
% original
% addendum = {' GCaMP6s', '', '', '', '', '', '', '', ' jGCaMP7f', ' jGCaMP7s', 'jGCaMP7c', 'jGCaMP7b', 'XCaMP-Gf', 'XCaMP-G', 'XCaMP-Gf0'};
% new
addendum = {'GCaMP6s', 'GCaMP6f', 'jGCaMP7f', 'jGCaMP8f', 'jGCaMP8m', 'jGCaMP8s', 'jGCaMP8.712', 'jGCaMP8.543', 'jGCaMP8.707', 'jGCaMP8.455', 'jGCaMP7s', 'jGCaMP7c', 'jGCaMP7b', 'XCaMP-Gf', 'XCaMP-G', 'XCaMP-Gf0'};
% hits = {'10.693', '10.921','500.456', '500.686', '500.688', '500.712', '500.543', '500.707', '500.455', '10.1473', '10.1513', '10.1561', '538.1', '538.2', '538.3'};

% addendum = {' GCaMP6s', ' jGCaMP7f', ' jGCaMP7s', 'jGCaMP7c', 'jGCaMP7b', '', '', 'XCaMP-Gf', 'XCaMP-G', 'XCaMP-Gf0'};

normPlots_struct = [];
control_dff_mean = nanmean(controlMutant.df_fpeak_med,2);
control_dff_sterr = std(controlMutant.df_fpeak_med,0,2)/sqrt(size(controlMutant.df_fpeak_med,2));
control_SNR_mean = nanmean(controlMutant.SNR,2);
control_SNR_sterr = normalized_error(controlMutant.SNR, controlMutant.SNR, 2);
control_halfrise_mean = nanmean(controlMutant.rise_half_med,2);
control_halfrise_sterr = std(controlMutant.rise_half_med,0,2)/sqrt(size(controlMutant.rise_half_med,2));
control_timetopeak_mean = nanmean(controlMutant.timetopeak_med,2);
control_timetopeak_sterr = std(controlMutant.timetopeak_med,0,2)/sqrt(size(controlMutant.timetopeak_med,2));
control_halfdecay_mean = nanmean(controlMutant.decay_half_med,2);
control_halfdecay_sterr = std(controlMutant.decay_half_med,0,2)/sqrt(size(controlMutant.decay_half_med,2));

% add control to normPlots_struct
normPlots_struct(1).construct = addendum{1};
normPlots_struct(1).dff_mean = control_dff_mean;
normPlots_struct(1).dff_sterr = control_dff_sterr;
normPlots_struct(1).SNR_mean = control_SNR_mean;
normPlots_struct(1).SNR_sterr = control_SNR_sterr;
normPlots_struct(1).halfrise_mean = control_halfrise_mean;
normPlots_struct(1).halfrise_sterr = control_halfrise_sterr;
normPlots_struct(1).timetopeak_mean = control_timetopeak_mean;
normPlots_struct(1).timetopeak_sterr = control_timetopeak_sterr;
normPlots_struct(1).halfdecay_mean = control_halfdecay_mean;
normPlots_struct(1).halfdecay_sterr = control_halfdecay_sterr;


dff_fig = figure; errorbar(nAPs, control_dff_mean, control_dff_sterr, 'k-', 'linewidth', 2);
SNR_fig = figure; errorbar(nAPs, ones(length(nAPs),1), control_SNR_sterr, 'k-', 'linewidth', 2);
halfrise_fig = figure; errorbar(nAPs, ones(length(nAPs),1), control_halfrise_sterr, 'k-', 'linewidth', 2);
timetopeak_fig = figure; errorbar(nAPs, ones(length(nAPs),1), control_timetopeak_sterr, 'k-', 'linewidth', 2);
halfdecay_fig = figure; errorbar(nAPs, ones(length(nAPs),1), control_halfdecay_sterr, 'k-', 'linewidth', 2);

% plotLegend = {mutant_hits.construct};
% plotLegend = [plotLegend' addendum'];
plotLegend = addendum;
plotLegend = join(plotLegend);

control_f0 = controlMutant.f0';

% text of mean +/- sterr
disp('GCaMP6s DFF')
join([string(control_dff_mean(relevantAPs))  string(control_dff_sterr(relevantAPs))], '�')
disp('GCaMP6s F0')
join([string(mean(control_f0))  string(std(control_f0) / length(control_f0))], '�')
disp('GCaMP6s half-rise')
join([string(control_halfrise_mean(relevantAPs))  string(control_halfrise_sterr(relevantAPs))], '�')
disp('GCaMP6s time to peak')
join([string(control_timetopeak_mean(relevantAPs))  string(control_timetopeak_sterr(relevantAPs))], '�')
disp('GCaMP6s half-decay')
join([string(control_halfdecay_mean(relevantAPs))  string(control_halfdecay_sterr(relevantAPs))], '�')
disp('GCaMP6s SNR')
join([string(control_SNR_mean(relevantAPs))  string(control_SNR_sterr(relevantAPs))], '�')

% colors = {'m', 'b', 'r', 'c', 'y'};
for i = 2:length(mutant_hits)
    plotMutant = mutant_hits(i);
    mutant_dff_mean = nanmean(plotMutant.df_fpeak_med,2);
    mutant_dff_sterr = std(plotMutant.df_fpeak_med,0,2)/sqrt(size(plotMutant.df_fpeak_med,2));
    figure(dff_fig); hold on, errorbar(nAPs, mutant_dff_mean, mutant_dff_sterr, 'color', cMap(i-1,:), 'linewidth', 2);
    
    mutant_SNR_mean = nanmean(plotMutant.SNR,2);
    mutant_SNR_sterr = nanstd(plotMutant.SNR,0,2) / size(plotMutant.SNR,2);
    mutant_SNR_mean_norm = mutant_SNR_mean ./ control_SNR_mean;
    mutant_SNR_sterr_norm = normalized_error(plotMutant.SNR, controlMutant.SNR, 2);
    
    % SNR looks weird for 1AP, remove
    mutant_SNR_sterr(1) = 0;
    SNR = plotMutant.SNR;
    % figure, hist(SNR(8,:))
    figure(SNR_fig); hold on, errorbar(nAPs, mutant_SNR_mean_norm, mutant_SNR_sterr_norm, 'color', cMap(i-1,:), 'linewidth', 2);
    
    mutant_halfrise_mean = nanmean(plotMutant.rise_half_med,2);
    mutant_halfrise_sterr = nanstd(plotMutant.rise_half_med,0,2) / size(plotMutant.rise_half_med,2);
    mutant_halfrise_mean_norm = mutant_halfrise_mean ./ control_halfrise_mean;
    mutant_halfrise_sterr_norm = normalized_error(plotMutant.rise_half_med, controlMutant.rise_half_med, 2);
    figure(halfrise_fig); hold on, errorbar(nAPs, mutant_halfrise_mean_norm, mutant_halfrise_sterr_norm,'color', cMap(i-1,:), 'linewidth', 2);
    
    mutant_timetopeak_mean = nanmean(plotMutant.timetopeak_med,2);
    mutant_timetopeak_sterr = nanstd(plotMutant.timetopeak_med,0,2) / size(plotMutant.timetopeak_med,2);
    mutant_timetopeak_mean_norm = mutant_timetopeak_mean ./ control_timetopeak_mean;
    mutant_timetopeak_sterr_norm = normalized_error(plotMutant.timetopeak_med, controlMutant.timetopeak_med, 2);
    figure(timetopeak_fig); hold on, errorbar(nAPs, mutant_timetopeak_mean_norm, mutant_timetopeak_sterr_norm,'color', cMap(i-1,:), 'linewidth', 2);
    
    mutant_halfdecay_mean = nanmean(plotMutant.decay_half_med,2);
    mutant_halfdecay_sterr = nanstd(plotMutant.decay_half_med,0,2) / size(plotMutant.decay_half_med,2);
    mutant_halfdecay_mean_norm = mutant_halfdecay_mean ./ control_halfdecay_mean;
    mutant_halfdecay_sterr_norm = normalized_error(plotMutant.decay_half_med, controlMutant.decay_half_med, 2);
    figure(halfdecay_fig); hold on, errorbar(nAPs, mutant_halfdecay_mean_norm, mutant_halfdecay_sterr_norm,'color', cMap(i-1,:), 'linewidth', 2);
    
    mutant_f0_mean = mean(plotMutant.f0);
    mutant_f0_sterr = std(plotMutant.f0)./length(plotMutant.f0);
    
    % text of mean +/- sterr
    disp(plotMutant.construct)
    disp([plotMutant.construct ' DFF'])
    join([string(mutant_dff_mean(relevantAPs))  string(mutant_dff_sterr(relevantAPs))], '�')
    disp([plotMutant.construct ' F0'])
    join([string(mutant_f0_mean)  string(mutant_f0_sterr)], '�')
    disp([plotMutant.construct ' half-rise'])
    join([string(mutant_halfrise_mean(relevantAPs))  string(mutant_halfrise_sterr(relevantAPs))], '�')
    disp([plotMutant.construct ' time to peak'])
    join([string(mutant_timetopeak_mean(relevantAPs))  string(mutant_timetopeak_sterr(relevantAPs))], '�')
    disp([plotMutant.construct ' half-decay'])
    join([string(mutant_halfdecay_mean(relevantAPs))  string(mutant_halfdecay_sterr(relevantAPs))], '�')
    disp([plotMutant.construct ' SNR'])
    join([string(mutant_SNR_mean(relevantAPs))  string(mutant_SNR_sterr(relevantAPs))], '�')
    
    % add mutant to normPlots_struct
    normPlots_struct(i).construct = addendum{i};
    normPlots_struct(i).dff_mean = mutant_dff_mean;
    normPlots_struct(i).dff_sterr = mutant_dff_sterr;
    normPlots_struct(i).SNR_mean = mutant_SNR_mean;
    normPlots_struct(i).SNR_sterr = mutant_SNR_sterr;
    normPlots_struct(i).halfrise_mean = mutant_halfrise_mean;
    normPlots_struct(i).halfrise_sterr = mutant_halfrise_sterr;
    normPlots_struct(i).timetopeak_mean = mutant_timetopeak_mean;
    normPlots_struct(i).timetopeak_sterr = mutant_timetopeak_sterr;
    normPlots_struct(i).halfdecay_mean = mutant_halfdecay_mean;
    normPlots_struct(i).halfdecay_sterr = mutant_halfdecay_sterr;
end

figure(dff_fig); 
% legend(plotLegend)
xticks(nAPs); xlim([.8 190])
set(gca, 'XScale', 'log'); box off
ylabel('DFF')
figure(SNR_fig); 
xticks(nAPs); xlim([.8 190])
ylabel('SNR')
legend(plotLegend)
set(gca, 'XScale', 'log'); box off
figure(halfrise_fig); 
xticks(nAPs); xlim([.8 190])
ylabel('half-rise')
% legend(plotLegend)
set(gca, 'XScale', 'log'); box off
figure(timetopeak_fig); 
xticks(nAPs); xlim([.8 190])
ylabel('time to peak')
% legend(plotLegend)
set(gca, 'XScale', 'log'); box off
figure(halfdecay_fig); 
xticks(nAPs); xlim([.8 190])
ylabel('half-decay')
% legend(plotLegend)
set(gca, 'XScale', 'log'); box off

% apply 'default' style
sdf(dff_fig, 'default')
sdf(SNR_fig, 'default')
sdf(halfrise_fig, 'default')
sdf(halfdecay_fig, 'default')

% save figs and pdfs
if saveOn
    saveas(dff_fig, fullfile(saveFolder, 'dff.fig'));
    saveas(SNR_fig, fullfile(saveFolder, 'SNR.fig'));
    saveas(halfrise_fig, fullfile(saveFolder, 'halfrise.fig'));
    saveas(timetopeak_fig, fullfile(saveFolder, 'timetopeak.fig'));
    saveas(halfdecay_fig, fullfile(saveFolder, 'halfdecay.fig'));
    
    saveas(dff_fig, fullfile(saveFolder, 'dff.pdf'));
    saveas(SNR_fig, fullfile(saveFolder, 'SNR.pdf'));
    saveas(halfrise_fig, fullfile(saveFolder, 'halfrise.pdf'));
    saveas(timetopeak_fig, fullfile(saveFolder, 'timetopeak.pdf'));
    saveas(halfdecay_fig, fullfile(saveFolder, 'halfdecay.pdf'));
    
    % save DFF inset figure
    figure(dff_fig)
    xlim([0.870   3.709])
    ylim([-0.28   3.24])
    saveas(dff_fig, fullfile(saveFolder, 'dff_inset.fig'))
    saveas(dff_fig, fullfile(saveFolder, 'dff_inset.pdf'))
    
    
    % f0 values to import into prism
    f0_ID = fopen(fullfile(saveFolder, 'f0.txt'),'w');
    for i = 1:length(plotLegend)
        fprintf(f0_ID, plotLegend{i});
        fprintf(f0_ID, '\n');
        fprintf(f0_ID, '%f\n', mutant_hits(i).f0');
        fprintf(f0_ID, '\n');
    end
    fclose(f0_ID);
end

% save struct for plotting in plotly
save('normPlots.mat', 'normPlots_struct', 'nAPs')
%% testing f0
% mngGECO 1374
% 6s: 1302.4355�25.5401
% 7b: 3673