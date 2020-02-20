base = 'Z:/';
if isempty(whos('mutant'))
    load(fullfile(base,'GECIScreenData\Analysis\pile_all_GCaMP96uf_upto_20200212.mat'), 'mutant')

end

newDateToken = {'2020127','20200120','20191202','20191125'};
oldDateToken = {'20190826', '20190819', '20190812', '20190813'};
APidx = 1; % 1:1AP, 2: 3AP, 3: 10AP, 4: 160AP

figure
currentControl = mutant(strcmp({mutant.construct}, '10.921'));
latest_dff_currentControl = currentControl.rise_half_med(:,startsWith(currentControl.date, newDateToken));
early_dff_currentControl = currentControl.rise_half_med(:,startsWith(currentControl.date, oldDateToken));

latest_nCells_currentControl = currentControl.nSegment(startsWith(currentControl.date, newDateToken));
early_nCells_currentControl = currentControl.nSegment(startsWith(currentControl.date, oldDateToken));


latest_dff = latest_dff_currentControl(APidx,:);
early_dff = early_dff_currentControl(APidx,:);

% subplot(1,2,1)
title('df/f')
scatter(rand(length(latest_dff),1)+.5, latest_dff)
hold on, scatter(rand(length(early_dff),1)+2.5, early_dff)
xlim([-1 5])

% subplot(1,2,2)
% title('# cells')
% scatter(rand(length(latest_nCells_currentControl),1)+.5, latest_nCells_currentControl)
% hold on, scatter(rand(length(early_nCells_currentControl),1)+2.5, early_nCells_currentControl)
% xlim([-1 5])
