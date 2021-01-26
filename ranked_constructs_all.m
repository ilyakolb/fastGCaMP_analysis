% generate ranked construct plots (as in Dana et al 2019 1c) but rank every
% plot

clearvars -except good
clc
% close all

base = 'Z:/';

control = '10.641';
hits = {'10.921', '500.456', '500.686', '500.688', '10.1473', '10.1513', '10.1561', '538.1', '538.2', '538.3'};
% hits = {'10.921', '500.456', '500.688'};
if isempty(whos('good'))
    good = readtable(fullfile(base,'\GECIScreenData\Analysis\data_all_20200325_GCaMP96uf.xlsx'));
end

switch width(good)
    case 52
        load varNames_4AP % old xlsx files without TimeToPeak variable
    case 62
        load varNames_4AP_withTimeToPeak % new xlsx files with TimeToPeak variable
    case 84
        load varNames_8AP % 8AP data for mng-GECO analysis only
    otherwise
        error('Wrong number of columns')
end

good.Properties.VariableNames(~startsWith(good.Properties.VariableNames, 'Var'))= varNames;

cMap = getColorMap(length(hits));

% remove bad stuff here
good_filt = good;

good_filt = good_filt(good_filt.replicate_number > 2, :); % at least 2 wells
% good_filt = good_filt(good_filt.x1_fp > 1, :); % > GCaMP6s response
good_filt = good_filt(~contains(good_filt.construct, {'TE', 'none'}), :); % remove empty wells
good_filt = good_filt(good_filt.rise_1_fp > 0, :); % remove oddities
% good_filt = good(~contains(good.construct, '410.'), :);




nConstructs = size(good_filt,1);
varNamesToPlot = { 'x1_fp', 'rise_1_fp', 'decay_1_fp', 'timetopeak_1_fp'};
nVars = length(varNamesToPlot);
figure 

for i = 1:nVars
    currentVarName = varNamesToPlot{i};
    
    % if kinetics, sort by descending (low values = fast = good)
    if contains(currentVarName, 'timetopeak') || contains(currentVarName, 'decay') || contains(currentVarName, 'rise')
        ascendOrDescent = 'descend';
    else
        ascendOrDescent = 'ascend';
    end
    [rankedTable, rankedIdx] = sortrows(good_filt, currentVarName, ascendOrDescent);
    tableVar = rankedTable.(currentVarName);
    
    subplot(length(varNamesToPlot), 1, i), bar(tableVar, 1)
    title(currentVarName, 'Interpreter', 'none')
    
    hold on, plot([0 nConstructs], [1 1], 'r-', 'linewidth' ,1)
    
    % plot black marker for control
    scatter(find(contains(rankedTable.construct, control)), rankedTable(contains(rankedTable.construct, control),:).(currentVarName), [], [0 0 0]);
    text(find(contains(rankedTable.construct, control)), rankedTable(contains(rankedTable.construct, control),:).(currentVarName), control);

    % plot markers for control and relevant hits here
    for j = 1:length(hits)
        currentHit = hits{j};
        currentHitIdx = contains(rankedTable.construct, hits{j});
        scatter(find(currentHitIdx), rankedTable(currentHitIdx,:).(currentVarName), [], cMap(j,:));
        text(find(currentHitIdx), rankedTable(currentHitIdx,:).(currentVarName), currentHit)
    end
end

hits = {};
tail(rankedTable)

function [subset] = getBoundedSubset(g, series, startIdx, endIdx)
% getBoundedSubset: return 'good' array matching a particular series,
% bounded by startIndx and endIdx
% g: good table
% series (str): series should end with . i.e. '500.'
% startIdx: (int) index with which to start/stop matching
% example: getBoundedSubset(g, '500.', 16, 27) returns table with
% constructs 500.16 - 500.27

% first match the series

subset = g(contains(g.construct, series), :);

idx_array = startIdx:endIdx;
idx_array_str = arrayfun(@(x) ['.' num2str(x)], idx_array, 'UniformOutput', 0);

% then match the startIdx:endIdx array
subset = subset(contains(subset.construct, idx_array_str));

end
% comboVar = good.df_f;
% [group, id] = findgroups(good.construct);
% groupFunc = @(p) [nanmean(p) nanstd(p)];
% comboData = splitapply(groupFunc, comboVar, group);
% groupT = table(id, comboData(:,1), comboData(:,2), 'VariableNames', {'construct', 'dff_mean', 'dff_std'});
% groupT = sortrows(groupT, 'dff_mean', 'ascend');
% figure,errorbar(groupT.dff_mean, groupT.dff_std)