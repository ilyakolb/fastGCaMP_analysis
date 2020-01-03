%% Run this after running compare_constructs_GCaMP96uf
close all
clc

% nAP vs parameter plots (like jGCaMP7 paper)
nAPs = [1 2 3 5 10 20 40 160];
relevantAPs = [1 3 5 8];

control_dff_mean = nanmean(controlMutant.df_fpeak_med,2);
control_dff_sterr = std(controlMutant.df_fpeak_med,0,2)/sqrt(size(controlMutant.df_fpeak_med,2));
control_SNR_mean = nanmean(controlMutant.SNR,2);
% control_SNR_sterr = std(controlMutant.SNR,0,2)/sqrt(size(controlMutant.SNR,2));
control_SNR_sterr = normalized_error(controlMutant.SNR, controlMutant.SNR, 2);
control_halfrise_mean = nanmean(controlMutant.rise_half_med,2);
control_halfrise_sterr = std(controlMutant.rise_half_med,0,2)/sqrt(size(controlMutant.rise_half_med,2));
control_halfdecay_mean = nanmean(controlMutant.decay_half_med,2);
control_halfdecay_sterr = std(controlMutant.decay_half_med,0,2)/sqrt(size(controlMutant.decay_half_med,2));

dff_fig = figure; errorbar(nAPs, control_dff_mean, control_dff_sterr, 'g-', 'linewidth', 2);
SNR_fig = figure; errorbar(nAPs, ones(length(nAPs),1), control_SNR_sterr, 'g-', 'linewidth', 2);
halfrise_fig = figure; errorbar(nAPs, ones(length(nAPs),1), control_halfrise_sterr, 'g-', 'linewidth', 2);
halfdecay_fig = figure; errorbar(nAPs, ones(length(nAPs),1), control_halfdecay_sterr, 'g-', 'linewidth', 2);

plotLegend = {mutant.construct};
plotLegend = plotLegend(1:end-1);
addendum = {' GCaMP6s', 'jGCaMP7f', 'jGCaMP7s', 'jGCaMP7c', 'jGCaMP7b', '545.1'};
plotLegend = [plotLegend' addendum'];
plotLegend = join(plotLegend);

control_f0 = controlMutant.f0';

% text of mean +/- sterr
disp('GCaMP6s DFF')
join([string(control_dff_mean(relevantAPs))  string(control_dff_sterr(relevantAPs))], '±')
disp('GCaMP6s F0')
join([string(mean(control_f0))  string(std(control_f0) / length(control_f0))], '±')
disp('GCaMP6s half-rise')
join([string(control_halfrise_mean(relevantAPs))  string(control_halfrise_sterr(relevantAPs))], '±')
disp('GCaMP6s half-decay')
join([string(control_halfdecay_mean(relevantAPs))  string(control_halfdecay_sterr(relevantAPs))], '±')
disp('GCaMP6s SNR')
join([string(control_SNR_mean(relevantAPs))  string(control_SNR_sterr(relevantAPs))], '±')

colors = {'m', 'b', 'r', 'c', 'y'};
for i = 2:length(mutant)-1
    plotMutant = mutant(i);
    mutant_dff_mean = nanmean(plotMutant.df_fpeak_med,2);
    mutant_dff_sterr = std(plotMutant.df_fpeak_med,0,2)/sqrt(size(plotMutant.df_fpeak_med,2));
    figure(dff_fig); hold on, errorbar(nAPs, mutant_dff_mean, mutant_dff_sterr, [colors{i-1} '-'], 'linewidth', 2);
    
    mutant_SNR_mean = nanmean(plotMutant.SNR,2);
    mutant_SNR_sterr = nanstd(plotMutant.SNR,0,2) / size(plotMutant.SNR,2);
    mutant_SNR_mean_norm = mutant_SNR_mean ./ control_SNR_mean;
    mutant_SNR_sterr_norm = normalized_error(plotMutant.SNR, controlMutant.SNR, 2);
    
    % SNR looks weird for 1AP, remove
    %mutant_SNR_sterr(1) = 0;
%     SNR = plotMutant.SNR;
%     figure, hist(SNR(8,:))
    figure(SNR_fig); hold on, errorbar(nAPs, mutant_SNR_mean_norm, mutant_SNR_sterr_norm, [colors{i-1} '-'], 'linewidth', 2);
    
    mutant_halfrise_mean = nanmean(plotMutant.rise_half_med,2);
    mutant_halfrise_sterr = nanstd(plotMutant.rise_half_med,0,2) / size(plotMutant.rise_half_med,2);
    mutant_halfrise_mean_norm = mutant_halfrise_mean ./ control_halfrise_mean;
    mutant_halfrise_sterr_norm = normalized_error(plotMutant.rise_half_med, controlMutant.rise_half_med, 2);
    figure(halfrise_fig); hold on, errorbar(nAPs, mutant_halfrise_mean_norm, mutant_halfrise_sterr_norm, [colors{i-1} '-'], 'linewidth', 2);
    
    mutant_halfdecay_mean = nanmean(plotMutant.decay_half_med,2);
    mutant_halfdecay_sterr = nanstd(plotMutant.decay_half_med,0,2) / size(plotMutant.rise_half_med,2);
    mutant_halfdecay_mean_norm = mutant_halfdecay_mean ./ control_halfdecay_mean;
    mutant_halfdecay_sterr_norm = normalized_error(plotMutant.decay_half_med, controlMutant.decay_half_med, 2);
    figure(halfdecay_fig); hold on, errorbar(nAPs, mutant_halfdecay_mean_norm, mutant_halfdecay_sterr_norm, [colors{i-1} '-'], 'linewidth', 2);
    
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
    disp([plotMutant.construct ' half-decay'])
    join([string(mutant_halfdecay_mean(relevantAPs))  string(mutant_halfdecay_sterr(relevantAPs))], '±')
    disp([plotMutant.construct ' SNR'])
    join([string(mutant_SNR_mean(relevantAPs))  string(mutant_SNR_sterr(relevantAPs))], '±')
    

end

figure(dff_fig)
legend(plotLegend)
xticks(nAPs); xlim([0 190])
set(gca, 'XScale', 'log')
ylabel('DFF')
figure(SNR_fig)
xticks(nAPs); xlim([0 190])
ylabel('SNR')
legend(plotLegend)
set(gca, 'XScale', 'log')
figure(halfrise_fig)
xticks(nAPs); xlim([0 190])
ylabel('half-rise')
legend(plotLegend)
set(gca, 'XScale', 'log')
figure(halfdecay_fig)
xticks(nAPs); xlim([0 190])
ylabel('half-decay')
legend(plotLegend)
set(gca, 'XScale', 'log')

% f0 plot into prism
% look at mutant.f0'

%% testing f0
% mngGECO 1374
% 6s: 1302.4355±25.5401
% 7b: 3673