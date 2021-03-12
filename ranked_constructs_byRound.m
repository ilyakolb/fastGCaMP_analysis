% generate ranked construct plots (as in Dana et al 2019 1c) but rank every
% plot


clearvars -except good
clc
% close all

hits_as_indiv_points = 1; % if 0, hits are plotted as points. if 1, hits plotted as horizontal lines
rank_by_rounds = 0;
saveFig = 0;

min_dff = 0;
max_rise = 4;
min_rise = 0.1;
min_decay = 0.01;
timetopeak_max = 3;

if hits_as_indiv_points
    hit_plot_str = '_indivpoints';
else
    hit_plot_str = '_horizline';
end
if rank_by_rounds
    pdf_dir = ['D:\ufgcamp_paper_data\culture-screen-figs\ranked\ranked_constructs_byRound' hit_plot_str '.pdf'];
else
    pdf_dir = ['D:\ufgcamp_paper_data\culture-screen-figs\ranked\ranked_constructs_all' hit_plot_str '.pdf'];
end
base = 'Z:/';
varNamesToPlot = { 'x1_fp', 'rise_1_fp', 'decay_1_fp', 'timetopeak_1_fp'};
varNameTitles = {'\DeltaF/F (norm.)', 'Rise time (norm.)', 'Decay time (norm.)', 'Peak time (norm.)'};
hits = {'10.641', '10.693', '10.921', '500.456', '500.686', '500.688', '10.1473', '10.1513', '10.1561', '538.1', '538.2', '538.3'};
hits_labels = {'GCaMP6s', 'GCaMP6f', 'jGCaMP7f', 'jGCaMP8f', 'jGCaMP8m', 'jGCaMP8s', 'jGCaMP7s', 'jGCaMP7c', 'jGCaMP7b', 'XCaMP-Gf', 'XCaMP-G', 'XCaMP-Gf0'};
% hits = {'10.921', '500.456', '500.688'};

labelAngle = 45;
if isempty(whos('good'))
    good = readtable(fullfile(base,'\GECIScreenData\Analysis\data_all_20210203_GCaMP96uf.xlsx'));
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

cMap = getColorMap(length(hits)); % for hits

% remove bad stuff here
good_filt = good;

% good_filt = good_filt(good_filt.replicate_number > 2, :); % at least 2 wells
% good_filt = good_filt(good_filt.x1_fp > 1, :); % > GCaMP6s response
good_filt = good_filt(~contains(good_filt.construct, {'TE', 'none'}), :); % remove empty wells
good_filt = good_filt(good_filt.rise_1_fp > 0, :); % remove oddities
good_filt = good_filt(good_filt.decay_1_fp > 0.05, :); % decays < 0.05 typically erroneous
good_filt = good_filt(good_filt.timetopeak_1_fp > 0.05, :); % decays < 0.05 typically erroneous
good_filt = good_filt(good_filt.timetopeak_1_fp < 3, :); % time to peak > 3 typically erroneous
% good_filt = good(~contains(good.construct, '410.'), :);

nConstructs = size(good_filt,1);

disp(['Total constructs: ' int2str(height(good))])
disp(['Filtered constructs: ' int2str(nConstructs)])


% sort by screening rounds as follows: 
% Round 0: 410.1-410.31.         Different peptides
% Round 1: 410.33-410.96        mainly linker1 screening for two different variants with either 1NIW or 1YR5 peptide
% Round 2: 500.2-500.303.       Single mutations
% Round 3: 500.306-500.374
% Round 4: 500.375-500.646
% Round 5: 500.647-500.671
% Round 6: 500.672-500.722
good_round_0 = getBoundedSubset(good_filt, '410.', 1, 31);
good_round_1 = getBoundedSubset(good_filt, '410.', 33, 96);
good_round_2 = getBoundedSubset(good_filt, '500.', 2, 303);
good_round_3 = getBoundedSubset(good_filt, '500.', 306, 374);
good_round_4 = getBoundedSubset(good_filt, '500.', 375, 646);
good_round_5 = getBoundedSubset(good_filt, '500.', 647, 671);
good_round_6 = getBoundedSubset(good_filt, '500.', 672, 722);
good_10_controls = getBoundedSubset(good_filt, '10.', 1, 2000); % get all GCaMP controls
good_XCaMP_controls = getBoundedSubset(good_filt, '538.', 1, 3); % get all GCaMP controls + XCaMP

if rank_by_rounds
    all_good_rounds = {good_round_0,good_round_1,good_round_2,good_round_3,good_round_4,good_round_5,good_round_6, good_10_controls, good_XCaMP_controls};
else
    all_good_rounds = {good_filt};
end
nRounds = length(all_good_rounds);

nVars = length(varNamesToPlot);
f = figure('position', [29          98        1803         880]);
f.Renderer='Painters';

cMap = getColorMap(length(hits));

startIdx = 0; % x index to start the bar plot
for j = 1:nRounds
    good_subset = all_good_rounds{j};
    nConstructs_round = height(good_subset);
    for i = 1:nVars
        currentVarName = varNamesToPlot{i};
        
        % if kinetics, sort by descending (low values = fast = good)
        if contains(currentVarName, 'timetopeak') || contains(currentVarName, 'decay') || contains(currentVarName, 'rise')
            ascendOrDescent = 'descend';
        else
            ascendOrDescent = 'ascend';
        end
        [rankedTable, rankedIdx] = sortrows(good_subset, currentVarName, ascendOrDescent);
        tableVar = rankedTable.(currentVarName);
        
        bar_idx = startIdx:(startIdx + nConstructs_round-1);
        subplot(length(varNamesToPlot), 1, i)
        
        if rank_by_rounds
            barcolor = 1 - 0.1 * j * ones(3,1);
        else
            barcolor = [0 0 0];
        end
        hold on, bar(bar_idx, tableVar, 1, 'FaceColor', barcolor, 'EdgeColor', barcolor)
        title(varNameTitles{i})
        
        if j == nRounds
            hold on, plot([0 nConstructs], [1 1], 'b-', 'linewidth' ,1)
        end
        % plot markers for control and relevant hits here
        
        for k = 1:length(hits)
            currentHit = hits{k};
            currentHitIdx = contains(rankedTable.construct, hits{k});
            
            if sum(currentHitIdx > 0)
                
                if hits_as_indiv_points
                    % hits as individual points
                    scatter(startIdx-1 + find(currentHitIdx), rankedTable(currentHitIdx,:).(currentVarName), 20, [1 0 0], 'filled', 'v');
                    
                    text(startIdx-1 + find(currentHitIdx), ...
                        rankedTable(currentHitIdx,:).(currentVarName), ...
                        ['  ' hits_labels{k}], ...
                        'Rotation', labelAngle,...
                        'FontSize', 10);
                else
                    % hits as horizontal lines
                    hold on, plot([0 nConstructs], [1, 1] * rankedTable(currentHitIdx,:).(currentVarName), 'color', cMap(k,:), 'linewidth', 1)
                    text(nConstructs, ...
                        rankedTable(currentHitIdx,:).(currentVarName), ...
                        ['  ' hits_labels{k}], ...
                        'FontSize', 9,...
                        'color', cMap(k,:));
                end

            end
        end
    end
    
    % offset start index for next round
    startIdx = startIdx + nConstructs_round;
end

if saveFig
    print(pdf_dir, '-dpdf', '-fillpage')
end