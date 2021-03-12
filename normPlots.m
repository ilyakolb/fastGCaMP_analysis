%% Run this after running compare_constructs_GCaMP96uf
clc
close all

% set to 1 to save figures
saveOn = 0;

saveFolder = 'D:\Dropbox (HHMI)\janelia\writing\jGCaMP8 patent\figures\all\';
% nAP vs parameter plots (like jGCaMP7 paper)
nAPs = [1 2 3 5 10 40];
relevantAPs = [1 2 3 4 5 6];

% filter only mutants that are in 'hits'
mutant_hits = controlMutant;
for i = 1:length(hits)
    mutant_hits = [mutant_hits mutant(contains({mutant.construct}, hits{i}))];
end
% addendum: string to add on to construct names for clarity
% original
% addendum = {' GCaMP6s', '', '', '', '', '', '', '', ' jGCaMP7f', ' jGCaMP7s', 'jGCaMP7c', 'jGCaMP7b', 'XCaMP-Gf', 'XCaMP-G', 'XCaMP-Gf0'};
% new
addendum = {'GCaMP6s', 'jGCaMP7f', 'jGCaMP8f', 'jGCaMP8m', 'jGCaMP8s', 'jGCaMP8.712', 'jGCaMP8.543', 'jGCaMP8.707', 'jGCaMP8.455', 'jGCaMP7s', 'jGCaMP7c', 'jGCaMP7b', 'XCaMP-Gf', 'XCaMP-G', 'XCaMP-Gf0'};
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
decays_not_nan = sum(~isnan(controlMutant.decay_half_med),2);
control_halfdecay_sterr = nanstd(controlMutant.decay_half_med,0,2)./sqrt(decays_not_nan);

% add control to unnormPlots_singleWells_struct
unnormPlots_singleWells_struct(1).construct = addendum{1};
unnormPlots_singleWells_struct(1).df_fpeak_med = controlMutant.df_fpeak_med;
unnormPlots_singleWells_struct(1).SNR = controlMutant.SNR;
unnormPlots_singleWells_struct(1).rise_half_med = controlMutant.rise_half_med;
unnormPlots_singleWells_struct(1).timetopeak_med = controlMutant.timetopeak_med;
unnormPlots_singleWells_struct(1).decay_half_med_comp = controlMutant.decay_half_med;

% add control to normPlots_struct
% for control, all means = 1
normPlots_struct(1).construct = addendum{1};
normPlots_struct(1).dff_mean = control_dff_mean;
normPlots_struct(1).dff_sterr = control_dff_sterr;
normPlots_struct(1).SNR_mean_norm = ones(length(nAPs),1);
normPlots_struct(1).SNR_sterr_norm = control_SNR_sterr;
normPlots_struct(1).halfrise_mean_norm = ones(length(nAPs),1);
normPlots_struct(1).halfrise_sterr_norm = control_halfrise_sterr;
normPlots_struct(1).timetopeak_mean_norm = ones(length(nAPs),1);
normPlots_struct(1).timetopeak_sterr_norm = control_timetopeak_sterr;
normPlots_struct(1).halfdecay_mean_norm = ones(length(nAPs),1);
normPlots_struct(1).halfdecay_sterr_norm = control_halfdecay_sterr;


dff_fig = figure; errorbar(nAPs, control_dff_mean, control_dff_sterr, 'k-', 'linewidth', 2);
SNR_fig = figure; errorbar(nAPs, ones(length(nAPs),1), control_SNR_sterr, 'k-', 'linewidth', 2);
halfrise_fig = figure; errorbar(nAPs, ones(length(nAPs),1), control_halfrise_sterr, 'k-', 'linewidth', 2);
timetopeak_fig = figure; errorbar(nAPs, ones(length(nAPs),1), control_timetopeak_sterr, 'k-', 'linewidth', 2);
halfdecay_fig = figure; errorbar(nAPs, ones(length(nAPs),1), control_halfdecay_sterr, 'k-', 'linewidth', 2);

plotLegend = addendum;
% plotLegend = join(plotLegend);

control_f0 = controlMutant.f0';

% text of mean +/- sterr
disp('GCaMP6s DFF')
join([string(control_dff_mean(relevantAPs))  string(control_dff_sterr(relevantAPs))], '±')
disp('GCaMP6s F0')
join([string(mean(control_f0))  string(std(control_f0) / length(control_f0))], '±')
disp('GCaMP6s half-rise')
join([string(control_halfrise_mean(relevantAPs))  string(control_halfrise_sterr(relevantAPs))], '±')
disp('GCaMP6s time to peak')
join([string(control_timetopeak_mean(relevantAPs))  string(control_timetopeak_sterr(relevantAPs))], '±')
disp('GCaMP6s half-decay')
join([string(control_halfdecay_mean(relevantAPs))  string(control_halfdecay_sterr(relevantAPs))], '±')
disp('GCaMP6s SNR')
join([string(control_SNR_mean(relevantAPs))  string(control_SNR_sterr(relevantAPs))], '±')

% colors = {'m', 'b', 'r', 'c', 'y'};
for i = 2:length(mutant_hits)
    plotMutant = mutant_hits(i);
    mutant_dff_mean = nanmean(plotMutant.df_fpeak_med,2);
    mutant_dff_sterr = std(plotMutant.df_fpeak_med,0,2)/sqrt(size(plotMutant.df_fpeak_med,2));
    figure(dff_fig); hold on, errorbar(nAPs, mutant_dff_mean, mutant_dff_sterr, 'color', cMap(i-1,:), 'linewidth', 2);
    
    mutant_SNR_mean = double(nanmean(plotMutant.SNR,2));
    mutant_SNR_sterr = double(nanstd(plotMutant.SNR,0,2)) / size(plotMutant.SNR,2);
    mutant_SNR_mean_norm = mutant_SNR_mean ./ control_SNR_mean;
    mutant_SNR_sterr_norm = normalized_error(double(plotMutant.SNR), controlMutant.SNR, 2);
    
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
    
    decays_not_nan = sum(~isnan(plotMutant.decay_half_med),2);
    mutant_halfdecay_sterr_norm = normalized_error(plotMutant.decay_half_med, controlMutant.decay_half_med, 2);
    figure(halfdecay_fig); hold on, errorbar(nAPs, mutant_halfdecay_mean_norm, mutant_halfdecay_sterr_norm,'color', cMap(i-1,:), 'linewidth', 2);
    
    mutant_f0_mean = mean(plotMutant.f0);
    mutant_f0_sterr = std(plotMutant.f0)./length(plotMutant.f0);
    
    % text of mean +/- sterr
    disp(plotMutant.construct)
    disp([plotMutant.construct ' DFF'])
    join([string(mutant_dff_mean(relevantAPs))  string(mutant_dff_sterr(relevantAPs))], '±')
    disp([plotMutant.construct ' F0'])
    join([string(mutant_f0_mean)  string(mutant_f0_sterr)], '±')
    disp([plotMutant.construct ' half-rise'])
    join([string(mutant_halfrise_mean(relevantAPs))  string(mutant_halfrise_sterr(relevantAPs))], '±')
    disp([plotMutant.construct ' time to peak'])
    join([string(mutant_timetopeak_mean(relevantAPs))  string(mutant_timetopeak_sterr(relevantAPs))], '±')
    disp([plotMutant.construct ' half-decay'])
    join([string(mutant_halfdecay_mean(relevantAPs))  string(mutant_halfdecay_sterr(relevantAPs))], '±')
    disp([plotMutant.construct ' SNR'])
    join([string(mutant_SNR_mean(relevantAPs))  string(mutant_SNR_sterr(relevantAPs))], '±')
    
    % add mutant to unnormPlots_singleWells_struct
    unnormPlots_singleWells_struct(i).construct = addendum{i};
    unnormPlots_singleWells_struct(i).df_fpeak_med = plotMutant.df_fpeak_med;
    unnormPlots_singleWells_struct(i).SNR = plotMutant.SNR;
    unnormPlots_singleWells_struct(i).rise_half_med = plotMutant.rise_half_med;
    unnormPlots_singleWells_struct(i).timetopeak_med = plotMutant.timetopeak_med;
    unnormPlots_singleWells_struct(i).decay_half_med_comp = plotMutant.decay_half_med;
    
    % add mutant to normPlots_struct
    normPlots_struct(i).construct = addendum{i};
    normPlots_struct(i).dff_mean = mutant_dff_mean;
    normPlots_struct(i).dff_sterr = mutant_dff_sterr;
    normPlots_struct(i).SNR_mean_norm = mutant_SNR_mean_norm;
    normPlots_struct(i).SNR_sterr_norm = mutant_SNR_sterr_norm;
    normPlots_struct(i).halfrise_mean_norm = mutant_halfrise_mean_norm;
    normPlots_struct(i).halfrise_sterr_norm = mutant_halfrise_sterr_norm;
    normPlots_struct(i).timetopeak_mean_norm = mutant_timetopeak_mean_norm;
    normPlots_struct(i).timetopeak_sterr_norm = mutant_timetopeak_sterr_norm;
    normPlots_struct(i).halfdecay_mean_norm = mutant_halfdecay_mean_norm;
    normPlots_struct(i).halfdecay_sterr_norm = mutant_halfdecay_sterr_norm;
    

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


% save figs and pdfs
if saveOn
    % apply 'default' style
    sdf(dff_fig, 'default')
    sdf(SNR_fig, 'default')
    sdf(halfrise_fig, 'default')
    sdf(halfdecay_fig, 'default')

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

save('plotly_normPlots.mat', 'normPlots_struct', 'nAPs')
save('unnormPlots_singleWells_struct.mat', 'unnormPlots_singleWells_struct')
%% testing f0
% mngGECO 1374
% 6s: 1302.4355±25.5401
% 7b: 3673