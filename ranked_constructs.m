% generate ranked construct plots (as in Dana et al 2019 1c)
clearvars -except good
clc
% close all

base = 'Z:/';

control = '10.641';
% hits = {'10.921', '500.456', '500.688', '10.1473', '10.1513', '10.1561', '538.1', '538.2', '538.3'};
hits = {'10.921', '500.456', '500.688'};
if isempty(whos('good'))
    good = readtable(fullfile(base,'\GECIScreenData\Analysis\data_all_20200308_GCaMP96uf.xlsx'));
end

switch width(good)
    case 52
        load varNames_4AP
    case 84
        load varNames_8AP
    otherwise
        error('Wrong number of columns')
end

good.Properties.VariableNames(~startsWith(good.Properties.VariableNames, 'Var'))= varNames;

cMap = getColorMap(length(hits));

% remove bad stuff here
good_filt = good(~contains(good.construct, '410.'), :);
good_filt = good_filt(good_filt.replicate_number > 2, :); % at least 2 wells
good_filt = good_filt(good_filt.x1_fp > 1, :); % > GCaMP6s response

% filter good table
rankVarName = 'rise_1_fp';
[rankedTable, rankedIdx] = sortrows(good_filt, rankVarName, 'descend');
nConstructs = size(rankedTable,1);
varNamesToPlot = { 'x1_fp', 'decay_1_fp', 'norm_f0'};
nOtherVars = length(varNamesToPlot);
figure('Position', [352          97        1242         875])
subplot(nOtherVars+1, 1, 1)
bar(table2array(rankedTable(:, rankVarName)), 1);
hold on, plot([0 nConstructs], [1 1], 'r-', 'linewidth' ,1)

% plot markers for hits here

for i = 1:length(hits)
    currentHit = hits{i};
    currentHitIdx = contains(rankedTable.construct, hits{i});
    scatter(find(currentHitIdx), rankedTable(currentHitIdx,:).(rankVarName), [], cMap(i,:));
    text(find(currentHitIdx), rankedTable(currentHitIdx,:).(rankVarName), currentHit)
end

% plot ordered list of all
for i = 1:nOtherVars
    currentVarName = varNamesToPlot{i};
    
    tableVar = rankedTable(:, currentVarName);
    tableVar = table2array(tableVar);
    
    subplot(length(varNamesToPlot)+1, 1, i+1), bar(tableVar, 1)
    title(currentVarName, 'Interpreter', 'none')
    
    hold on, plot([0 nConstructs], [1 1], 'r-', 'linewidth' ,1)
    
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

% comboVar = good.df_f;
% [group, id] = findgroups(good.construct);
% groupFunc = @(p) [nanmean(p) nanstd(p)];
% comboData = splitapply(groupFunc, comboVar, group);
% groupT = table(id, comboData(:,1), comboData(:,2), 'VariableNames', {'construct', 'dff_mean', 'dff_std'});
% groupT = sortrows(groupT, 'dff_mean', 'ascend');
% figure,errorbar(groupT.dff_mean, groupT.dff_std)