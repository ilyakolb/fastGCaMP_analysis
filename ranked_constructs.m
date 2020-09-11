% generate ranked construct plots (as in Dana et al 2019 1c)
clearvars -except good
clc
% close all

pareto = 0; % set to 1 to do pareto optimization
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
good_filt = good(~contains(good.construct, '410.'), :);
good_filt = good_filt(good_filt.replicate_number > 2, :); % at least 2 wells
good_filt = good_filt(good_filt.x1_fp > 1, :); % > GCaMP6s response

%% trying pareto optimization
% weights
if pareto
    wDFF = 1;
    wRise = -5;
    wTTP = -5;
    wDecay = -2;
    pareto = wDFF * good_filt.x1_fp + wRise * log(good_filt.rise_1_fp) + wTTP * log(good_filt.timetopeak_1_fp) + wDecay * log(good_filt.decay_1_fp);
    good_filt.pareto = pareto;
    %%
    % filter good table
    rankVarName = 'pareto';
    [rankedTable, rankedIdx] = sortrows(good_filt, rankVarName, 'ascend');
else
    rankVarName = 'rise_1_fp';
    [rankedTable, rankedIdx] = sortrows(good_filt, rankVarName, 'descend');
end 
nConstructs = size(rankedTable,1);
varNamesToPlot = { 'x1_fp', 'decay_1_fp', 'norm_f0', 'timetopeak_1_fp'};
nOtherVars = length(varNamesToPlot);
figure %('Position', [2034          47        1087         945])
subplot(nOtherVars+1, 1, 1)
bar(table2array(rankedTable(:, rankVarName)), 1);
hold on, plot([0 nConstructs], [1 1], 'r-', 'linewidth' ,1)

% plot black marker for control
scatter(find(contains(rankedTable.construct, control)), rankedTable(contains(rankedTable.construct, control),:).(rankVarName), [], [0 0 0]);
text(find(contains(rankedTable.construct, control)), rankedTable(contains(rankedTable.construct, control),:).(rankVarName), control);

% plot markers for hits
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

% comboVar = good.df_f;
% [group, id] = findgroups(good.construct);
% groupFunc = @(p) [nanmean(p) nanstd(p)];
% comboData = splitapply(groupFunc, comboVar, group);
% groupT = table(id, comboData(:,1), comboData(:,2), 'VariableNames', {'construct', 'dff_mean', 'dff_std'});
% groupT = sortrows(groupT, 'dff_mean', 'ascend');
% figure,errorbar(groupT.dff_mean, groupT.dff_std)