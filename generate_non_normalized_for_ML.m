% generate non-normalized data for ML
% 6/26/20
% ilya kolb
% loads pile_all file and generates csv datatable
% NOTE: loading pile_all file may take long time!
clearvars -except mutant

csvFileName = 'gcamp3unnormalized.csv'; % file to save

% csv column names
colNames = {'name', 'nReplicates', 'dff_1_mean', 'dff_3_mean', 'dff_10_mean', 'dff_160_mean',...
    'dff_1_std', 'dff_3_std', 'dff_10_std', 'dff_160_std', ...
    'rise_1_mean', 'rise_3_mean', 'rise_10_mean', 'rise_160_mean', ...
    'rise_1_std', 'rise_3_std', 'rise_10_std', 'rise_160_std',...
    'decay_1_mean', 'decay_3_mean', 'decay_10_mean', 'decay_160_mean', ...
    'decay_1_std', 'decay_3_std', 'decay_10_std', 'decay_160_std', ...
    'df_1_mean', 'df_3_mean', 'df_10_mean', 'df_160_mean', ...
    'df_1_std', 'df_3_std', 'df_10_std', 'df_160_std', 'f0'};
nCols = length(colNames);

% variable types for columns
varTypes = {'string'}';
varTypes = vertcat(varTypes, repmat({'double'}, nCols-1,1));

% generate csv of unnormalized data
if isempty(whos('mutant'))
    % gcamp6 dataset
    % load('Z:\GECIScreenData\Analysis\pile_all_GCaMP96b_upto_20161214.mat')
    % gcamp3 dataset
    load('Z:\GECIScreenData\Analysis\pile_all_GCaMP_upto_20130703.mat')
end

FPidx = [1 3 5 9]; % indices of 1,3,10,160 AP stims for GCaMP3 dataset
% FPidx = [1 2 3 4]; % indices of 1,3,10,160 AP stims for GCaMP6 dataset

%%
nRows = length(mutant);

T = table('Size', [nRows, nCols], 'VariableTypes', varTypes, 'VariableNames', colNames);
for i = 1:length(mutant)
    currentMutant = mutant(i);
    fullname = currentMutant.fullname;
    nReplicates = currentMutant.nreplicate;
    dff = currentMutant.df_fpeak_med(FPidx,:) + currentMutant.df_fpeak_med_comp(FPidx,:);
    nStims = size(dff,1);
    dff_mean = nanmean(dff, 2);
    dff_std = nanstd(dff, [], 2);
    f0_mean = nanmean(currentMutant.f0);
    
    rise_half = currentMutant.rise_half_med(FPidx,:); % + currentMutant.rise_half_med_comp(FPidx,:);
    rise_half_mean = nanmean(rise_half, 2);
    rise_half_std = nanstd(rise_half, [], 2);
    
    decay_half = currentMutant.decay_half_med(FPidx,:); % + currentMutant.decay_half_med_comp(FPidx,:);
    decay_half_mean = nanmean(decay_half, 2);
    decay_half_std = nanstd(decay_half, [], 2);
    
    df = dff.*repmat(f0_mean, size(dff));
    df_mean = nanmean(df, 2);
    df_std = nanstd(df, [], 2);
    
    T(i,:) = {['"' fullname '"'], nReplicates, dff_mean(1), dff_mean(2), dff_mean(3), dff_mean(4), ...
        dff_std(1), dff_std(2), dff_std(3), dff_std(4),...
        rise_half_mean(1), rise_half_mean(2), rise_half_mean(3), rise_half_mean(4),...
        rise_half_std(1), rise_half_std(2), rise_half_std(3), rise_half_std(4),...
        decay_half_mean(1), decay_half_mean(2), decay_half_mean(3), decay_half_mean(4),...
        decay_half_std(1), decay_half_std(2), decay_half_std(3), decay_half_std(4),...
        df_mean(1), df_mean(2), df_mean(3), df_mean(4), ...
        df_std(1), df_std(2), df_std(3), df_std(4), f0_mean};
end
writetable(T, csvFileName, 'QuoteStrings',true)

