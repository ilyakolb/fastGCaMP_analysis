base = 'Z:/';
if isempty(whos('mutant'))
    load(fullfile(base,'GECIScreenData\Analysis\pile_all_GCaMP96uf_upto_20191127.mat'), 'mutant')

end

dateToken = '20191111';
APidx = 1; % 1:1AP, 2: 3AP, 3: 10AP, 4: 160AP

figure
currentControl = mutant(strcmp({mutant.construct}, '10.921'));
latest_dff_currentControl = currentControl.f0(:,startsWith(currentControl.date, dateToken));
early_dff_currentControl = currentControl.f0(:,startsWith(currentControl.date, '2018'));

latest_nCells_currentControl = currentControl.nSegment(startsWith(currentControl.date, dateToken));
early_nCells_currentControl = currentControl.nSegment(startsWith(currentControl.date, '2018'));


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
