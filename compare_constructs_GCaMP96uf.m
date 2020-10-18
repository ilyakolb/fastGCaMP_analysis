%% Ultra-fast GCaMP hits analysis
% NOTES:
% 10/2/19
% loads mutants variable from latest data in nearline analysis folder e.g. /Volumes/genie/GENIE_Pipeline/Analysis/Results/pile_all_GCaMP96uf_upto_20190509.mat
% loads figure comparing all hits
% calls plotRawCellData to show figure for each well and optionally open Fiji to show raw videos
% 17/11/19
% adapting for 8 stim pulses
% usage:
% select hits you want to see (e.g. '500.331') and relevant control (e.g. '10.921') in the user params section
% select max number of wells to plot for each construct (numSampleWells)
% 
% plate controls:
% 10.641 GCaMP6s   : pGP-SIV-(SalI)-syn-IRES-nls-mCherry-WPRE-GCaMP3 K78H T302L R303P D380Y T381R S383T R392G.10.641
% 10.693 GCaMP6f   : pGP-SIV-(SalI)-syn-IRES-nls-mCherry-WPRE-GCaMP3 T302L R303P A317E D380Y T381R S383T R392G.10.693
% 10.921 GCaMP7f  : pGP-SIV-(SalI)-syn-IRES-nls-mCherry-WPRE-GCaMP3 T302L R303P A317L D380Y.10.921
% 10.1513 jGCaMP7c: pGP-SIV-(SalI)-syn-IRES-nls-mCherry-WPRE-GCaMP3 L59Q E60P T302L R303P M378G K379S D380Y T381R R392G T412N.10.1513
% 10.1473 jGCaMP7s
% 10.1561: jGCaMP7b pGP-SIV-(SalI)-syn-IRES-nls-mCherry-WPRE-GCaMP3 T302P R303P A317L M374Y D380Y T381R S383T R392G.10.1561
% XCaMP controls:
% 538.1: XCaMP-Gf
% 538.2: XCaMP-G
% 538.3: XCaMP-Gf0

%% user parameters

addpath('accessory_funcs\')
clearvars -except mutant
clc

rng('default'); % for reproducibility

% hits = {'10.921', '500.311', '500.330', '500.333', '500.336', '500.350', '500.378'};
% hits = {};

% all hits
hits = {'10.693', '10.921','500.456', '500.686', '500.688', '500.712', '500.543', '500.707', '500.455', '10.1473', '10.1513', '10.1561', '538.1', '538.2', '538.3'};

% all variants from 3-5-20 PPT slide except 640 + best performers + xcamps + 7 series (loaner + our camera, EM gain 25)
% hits = {'500.456', '500.688', '500.712', '500.543', '500.707', '500.455', '10.921', '10.1473', '10.1513', '10.1561', '538.1', '538.2', '538.3'};
% hits = {}
% all variants from best performers + xcamps + 7 series week (pile_week_GCaMP96uf_upto_20200310_GCaMP96uf_raw)
% substituting 10.641 for variants that are not in this batch as a hacky
% way to preserve color scheme
% hits = {'500.456', '10.641', '500.688', '500.712', '500.543', '500.707', '500.455', '10.921', '10.1473', '10.1513', '10.1561', '538.1', '538.2', '538.3'};

% hits = {'10.921', '500.456', '500.686', '500.688', '500.333', '500.640', '500.712', '500.543', '500.707', '500.455'};

% 6th round hits (dff, kinetics)
%hits = {'10.921', '500.456', '500.640', '500.686', '500.675', '500.676', '500.688'};

% 6th round hits (other params)
% hits = {'10.921', '500.456', '500.666', '500.543', '500.712'};


% 4th round (best 1AP)
% hits = {'10.921', '500.333', '500.456'};

% 4th round (best other params)
% hits = {'10.921', '500.333', '500.512', '500.543'};

% XCaMP comparison
% hits = {'10.921', '10.693', '538.1', '538.2', '538.3', '500.333'};

% XCaMP comparison (compared to our best)
% hits = {'500.456', '500.640', '538.1', '538.2', '538.3'};

% 5th round (20191125)
% hits = {'10.921', '500.333', '500.378', '500.668'};

%5th round (20191209, ilastik)
% hits = {'10.921', '500.333', '500.668', '500.663', '500.649'};

% 5th round (20191118, no 333)
% hits = {'10.921', '500.668', '500.650', '500.658'};

% abhi's sensors (8 APs)
% hits = {'10.921', '10.1473', '10.1513', '545.1'}; 

% hits = {'10.641', '500.333', '500.659'};
% hits = {'10.641'};
% hits = {'10.921', '500.333', '500.378'}; % hits from newest screen
% hits = {'500.469', '500.508'  , '500.460', '500.462'}; % hits from 4th round
%hits = {'10.921', '10.1473', '10.1513', '545.1a', '545.1b', '545.1c', '545.1d'}; % abhi's sensors
% hits = {'10.641', '500.465', '500.464', '500.455'}; % fast rise, fast decay, high f0 for flies

control= '10.641';

alignControlToStimPulse = 1; % 1 to correct for stim pulse timing variability in controls. takes longer time
alignMutantToStimPulse = 1;  % 1 to correct for stim pulse timing variability in mutants. takes longer time 
bleachCorrect = 1;           % 1 to bleach correct the 1FP traces
Fs = 200;                    % sampling rate (Hz) assuming GCaMPuf
plotRaw = 0;                 % 1 to plot raw well figures
numSampleWells =10;           % number of sample wells to plot
launchFiji = 0;              % 1 to launch Fiji and show every tiff stack
apNumIdx = 1;                % AP index for  (1, 3, 10, 160) to 
% plot colors
% col=['b','r','g','m','c','k', 'b','r','g','m','c','y'];
APstimNames = {'1AP', '3AP', '10AP', '160AP'};
% APstimNames = {'1AP', '2AP', '3AP', '5AP', '10AP', '20AP', '40AP', '160AP'};


%%
if bleachCorrect && (alignControlToStimPulse || alignMutantToStimPulse)
    warning('Doing bleach correction with stim alignment: expect first 100 ms of traces to be weird!')
end

nStims = length(APstimNames);
base = 'Z:/';
if ismac
	base = fullfile('/Volumes', 'genie');
        error('Code can''t run on Mac due to old MATLAB version!')
end

if isempty(whos('mutant'))
    % load latest MAT
    
    % ALL including best performers + xcamps + 7
    load(fullfile(base,'GECIScreenData\Analysis\pile_all_GCaMP96uf_upto_20200325.mat'), 'mutant')
    
    % best performers + xcamps + 7 series (loaner + our camera, EM gain 25)
    % load(fullfile(base,'GECIScreenData\Analysis\pile_week_GCaMP96uf_upto_20200310_GCaMP96uf_analyzed.mat'), 'mutant')
    
    % 6th round ONLY with fixed jgcamp7f control
    % load(fullfile(base,'GECIScreenData\Analysis\pile_week_GCaMP96uf_upto_20200303_GCaMP96uf_analyzed.mat'), 'mutant')

    % ALL after week 2 of 6th round (updated ilastik parameters)
    % load(fullfile(base,'GECIScreenData\Analysis\pile_all_GCaMP96uf_upto_20200308.mat'), 'mutant')

    % ALL after week 2 of 6th round (old ilastik parameters)
    % load(fullfile(base,'GECIScreenData\Analysis\pile_all_GCaMP96uf_upto_20200212.mat'), 'mutant')


    % 6th round 20200211 (week 2)
    % load(fullfile(base,'GECIScreenData\Analysis\pile_week_GCaMP96uf_upto_20200211_GCaMP96uf_raw.mat'), 'mutant')
    
    % 6th round 20200205 (week 1)
    % load(fullfile(base,'GECIScreenData\Analysis\pile_week_GCaMP96uf_upto_20200205_GCaMP96uf_analyzed.mat'), 'mutant')
    
    % all up to 5th round w/ ilastik
    % load(fullfile(base,'GECIScreenData\Analysis\pile_all_GCaMP96uf_upto_20200128.mat'), 'mutant')
    
    % 5th round w/ ilastik
    % load(fullfile(base,'GECIScreenData\Analysis\pile_week_GCaMP96uf_upto_20191209_GCaMP96uf_raw.mat'), 'mutant')
    
    
    % ufGCaMPs round 1-4 with fixed pixels, Hod's cellfinder (not ilastik)
    % load(fullfile(base,'GECIScreenData\Analysis\pile_all_GCaMP96uf_upto_20191211_OLDCELLFINDER.mat'), 'mutant')
    
    % Abhi's sensor from 20191112 (8 AP stims) BUG CORRECTED
    % load(fullfile(base,'GECIScreenData\Analysis\pile_week_mngGECO_upto_20191112_GCaMP96uf_analyzed.mat'), 'mutant')
    
    % fixed pixels from Yan's 5th round
    % 20191118 data
    % load(fullfile(base,'GECIScreenData\Analysis\pile_week_GCaMP96uf_upto_20191118_GCaMP96uf_analyzed.mat'), 'mutant')
    
    % fixed pixels from Yan's 5th round
    % 20191125 data
    % load(fullfile(base,'GECIScreenData\Analysis\pile_week_GCaMP96uf_upto_20191125_GCaMP96uf_analyzed.mat'), 'mutant')
    
    
    % fixed pixels from Yan's 4th round (combined)
    % load(fullfile(base,'GECIScreenData\Analysis\pile_all_4thround.mat'), 'mutant')
    
    % fixed pixels from Yan's 4th round
    % 0910 data
    % load(fullfile(base,'GECIScreenData\Analysis\pile_week_GCaMP96uf_upto_20190910_GCaMP96uf_analyzed.mat'), 'mutant')
    
    % 0903 data
    % load(fullfile(base,'GECIScreenData\Analysis\pile_week_GCaMP96uf_upto_20190903_GCaMP96uf_analyzed.mat'), 'mutant')
    
    % 0827 data
    % load(fullfile(base,'GECIScreenData\Analysis\pile_week_GCaMP96uf_upto_20190827_GCaMP96uf_analyzed.mat'), 'mutant')

    % XCaMP from 20190904 analysis
    % load(fullfile(base,'GECIScreenData\Analysis\pile_week_GCaMP96uf_upto_20190904_GCaMP96uf_analyzed.mat'), 'mutant')
    
    % XCaMP from 10-24 analysis
    % load(fullfile(base,'GECIScreenData\Analysis\pile_week_GCaMP96uf_upto_20191005_GCaMP96uf_analyzed.mat'), 'mutant')
    
    % XCaMP from 10-28 analysis
    % load(fullfile(base,'GECIScreenData\Analysis\pile_week_GCaMP96uf_upto_20191022_GCaMP96uf_raw.mat'), 'mutant')
    
    % hits from 20190903 (where 456 variant was discovered)
    % 456 variant is in P5a-20190819_GCaMP96uf
    % load(fullfile(base,'GECIScreenData\Analysis\pile_week_GCaMP96uf_upto_20190903_GCaMP96uf_analyzed.mat'), 'mutant')
    
    % Abhi's sensor from 20191112 (8 AP stims)
    % load(fullfile(base,'GECIScreenData\Analysis\pile_week_mngGECO_upto_20191112_GCaMP96uf_analyzed.mat'), 'mutant')
    
    % Abhi's sensor from 10-24 analysis
    %load(fullfile(base,'GECIScreenData\Analysis\pile_week_GCaMP96uf_upto_20191006_GCaMP96uf_analyzed.mat'), 'mutant')

end

tableVarNames = {'construct', 'nWells', 'df_f_AP_mean', 'df_f_AP_std', 'rise_half_AP_ms_mean', 'rise_half_AP_ms_std', 'rise_full_AP_ms_mean', 'rise_full_AP_ms_std',...
    'f0_mean', 'f0_std', 'decay_half_med_mean', 'decay_half_med_std', 'SNR', 'SNR_std'};
comparisonTable = table('Size', [length(hits)+1 length(tableVarNames)], 'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', tableVarNames);

% controlMutant = mutant(find(strcmp({mutant.construct},control)));
controlMutant = mutant(find(endsWith({mutant.construct},control)));

% number of ROIs for controlMutant and all mutants
nROI = {};
nROI{end+1} = controlMutant.nSegment;

comparisonTable(1,:) = {controlMutant.construct, controlMutant.nreplicate, nanmean(controlMutant.df_fpeak_med(apNumIdx,:)), nanstd(controlMutant.df_fpeak_med(apNumIdx,:)), ...
    nanmean(controlMutant.rise_half_med(apNumIdx,:)), nanstd(controlMutant.rise_half_med(apNumIdx,:)),...
    nanmean(controlMutant.timetopeak_med(apNumIdx,:)), nanstd(controlMutant.timetopeak_med(apNumIdx,:)), nanmean(controlMutant.f0), nanstd(controlMutant.f0),...
    nanmean(controlMutant.decay_half_med(apNumIdx,:) + controlMutant.decay_half_med_comp(apNumIdx,:)), nanstd(controlMutant.decay_half_med(apNumIdx,:)+controlMutant.decay_half_med_comp(apNumIdx,:)),...
    nanmean(controlMutant.SNR(apNumIdx,:)), nanstd(controlMutant.SNR(apNumIdx,:))};


if alignControlToStimPulse
    control_df_f_med = align_responses(controlMutant, nStims, 0);
else
    control_df_f_med = controlMutant.df_f_med;
end

time=(1:size(control_df_f_med,1))/Fs;%;%1/35:1/35:249/35;
control_df_f_med_mean=nanmean(control_df_f_med,3);

% bleach correct only 1FP trace!
if bleachCorrect
     control_df_f_med_mean(:,1) = bleachCorr(time,control_df_f_med_mean(:,1));
end

control_med_med_dff_sterr=std(control_df_f_med,0,3)/sqrt(size(control_df_f_med,3));


for i=1:length(hits)
    found_hits_idx = find(strcmp({mutant.construct},hits{1,i}));
    assert(~isempty(found_hits_idx), ['Hit ' hits{i} ' not found in pile_all mat file!'])
    hits_idx(i) = found_hits_idx;
end



for i=1:length(hits)
    currentMutant = mutant(1,hits_idx(i));
    
    if alignMutantToStimPulse
        mutant_df_f_med = align_responses(currentMutant, nStims, 0);
    else
        mutant_df_f_med = currentMutant.df_f_med;
    end
    
    hits_med_med_dff(:,:,i)=mean(mutant_df_f_med,3);
    % hits_med_med_dff_unaligned(:,:,i)=mean(currentMutant.df_f_med,3);
    hits_med_med_dff_sterr(:,:,i)=std(mutant_df_f_med,0,3)/sqrt(size(mutant_df_f_med,3));
    
    if plotRaw
        plotRawCellData(currentMutant,numSampleWells, launchFiji);
    end
    % populate stats table with variant
    comparisonTable(1+i, :) = {currentMutant.construct, currentMutant.nreplicate, nanmean(currentMutant.df_fpeak_med(apNumIdx,:)), nanstd(currentMutant.df_fpeak_med(apNumIdx,:)),...
    nanmean(currentMutant.rise_half_med(apNumIdx,:)), nanstd(currentMutant.rise_half_med(apNumIdx,:)),...
    nanmean(currentMutant.timetopeak_med(apNumIdx,:)), nanstd(currentMutant.timetopeak_med(apNumIdx,:)), nanmean(currentMutant.f0), nanstd(currentMutant.f0), ...
    nanmean(currentMutant.decay_half_med(apNumIdx,:) + currentMutant.decay_half_med_comp(apNumIdx,:)), nanstd(currentMutant.decay_half_med(apNumIdx,:)+currentMutant.decay_half_med_comp(apNumIdx,:)), ...
    nanmean(currentMutant.SNR(apNumIdx,:)), nanstd(currentMutant.SNR(apNumIdx,:))};

    nROI{end+1} = currentMutant.nSegment;
end

% bleach correct 1AP mutant traces -- BELONGS HERE?
if bleachCorrect
    for i = 1:length(hits)
        hits_med_med_dff(:,1,i) = bleachCorr(time, hits_med_med_dff(:,1,i));
    end
end

f = figure('name', 'comparison', 'position', [2165         616        1238         249]);
cMap = getColorMap(length(hits));

for n=1:nStims
    subplot(1,nStims,n)
    shadedErrorBar(time,control_df_f_med_mean(:,n),control_med_med_dff_sterr(:,n), 'lineprops', 'k-', 'transparent',1); % IK mod
    
    hold on
    for i=1:length(hits)
        shadedErrorBar(time,hits_med_med_dff(:,n,i),hits_med_med_dff_sterr(:,n,i),'lineprops', {'color', cMap(i,:)}, 'transparent',1);
    end

    set(gca,'FontSize',10)
    title(APstimNames{n})

end


legendWithControl = [control hits];
legend(legendWithControl)
set(gcf,'Visible','on')

% set x and y limits
subplot(1,4,1); xlim([.23 4.1]); ylim([-0.1471    1.3917])
subplot(1,4,2); xlim([.23 4.1]); ylim([-0.3724    2.7636])
subplot(1,4,3); xlim([.23 4.1]); ylim([-1 6.5])
subplot(1,4,4); xlim([.23 6])

% Table for Prism
% Format is: {Name, mean, st.err, num samples} for each parameter
% Parameters: DF_F, Half Rise Time, Full Rise Time, Half Decay Time, F0
% open this in Variable explorer, copy, and paste to Prism
% Prism options: New Table (column), Enter values calculated elsewhere,
% Mean, 
comparisonTable_forPrism = table(comparisonTable.construct, comparisonTable.df_f_AP_mean, comparisonTable.df_f_AP_std, comparisonTable.nWells, ...
    comparisonTable.rise_half_AP_ms_mean, comparisonTable.rise_half_AP_ms_std, comparisonTable.nWells, comparisonTable.rise_full_AP_ms_mean, comparisonTable.rise_full_AP_ms_std, comparisonTable.nWells, ...
    comparisonTable.decay_half_med_mean, comparisonTable.decay_half_med_std, comparisonTable.nWells, comparisonTable.f0_mean, comparisonTable.f0_std, comparisonTable.nWells);
% calculate statistics
disp('STATISTICS')
disp([num2str(length(mutant)) ' unique constructs']);

% save data for plotting
plot_out.control_med_med_dff = control_df_f_med_mean;
plot_out.hits_med_med_dff = hits_med_med_dff;
plot_out.time = time;
plot_out.hits = hits;
plot_out.control = control;
plot_out.hits_med_med_dff_sterr = hits_med_med_dff_sterr;
plot_out.control_med_med_dff_sterr = control_med_med_dff_sterr;

save('plotting.mat', 'plot_out')

% normPlots

