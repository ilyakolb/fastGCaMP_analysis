% Gcamp6s decay half med = 
%     0.0093
%     0.4310
%     1.2597
%     2.6416
% 
% 	
%     0.0074
%     0.0554
%     1.6672
%     0.7665

% APidx = 3; % 1: 1AP, 2: 3AP, 3: 10AP, 4: 160AP
% figure
% 
% dff_new_gcamp6s = control.df_fpeak_med(:,contains(control.plate, '201808'));
% dff_old_gcamp6s = control.df_fpeak_med(:,~contains(control.plate, '201808'));
% 
% subplot(2,1,1)
% histogram(dff_new_gcamp6s(APidx,:))
% % xlim([0 3])
% title('new')
% box off
% subplot(2,1,2)
% title('old')
% hold on, histogram(dff_old_gcamp6s(APidx,:))
% % xlim([0 3])
% 
% %% 7-30-19 GCaMP comparisons for Ron Vale presentation
% new_1fp = readtable('D:\Dropbox (HHMI)\janelia\meetings\illustrations\gcamp_redux_for_RV\500.311_1FP.csv');
% new_10fp = readtable('D:\Dropbox (HHMI)\janelia\meetings\illustrations\gcamp_redux_for_RV\500.311_10FP.csv');
% 
% new_1fp_y = (new_1fp.Y - min(new_1fp.Y)) / (mean(new_1fp.Y));
% new_10fp_y = (new_10fp.Y - min(new_10fp.Y)) / (mean(new_10fp.Y));
% 
% %
% gcamp6s_1fp = readtable('D:\Dropbox (HHMI)\janelia\meetings\illustrations\gcamp_redux_for_RV\gcamp6s_1FP.csv');
% gcamp6s_10fp = readtable('D:\Dropbox (HHMI)\janelia\meetings\illustrations\gcamp_redux_for_RV\gcamp6s_10FP.csv');
% 
% gcamp6s_1fp_y = (gcamp6s_1fp.Y - min(gcamp6s_1fp.Y)) / (mean(gcamp6s_1fp.Y));
% gcamp6s_10fp_y = (gcamp6s_10fp.Y - min(gcamp6s_10fp.Y)) / (mean(gcamp6s_10fp.Y));
% 
% 
% %
% gcamp7f_1fp = readtable('D:\Dropbox (HHMI)\janelia\meetings\illustrations\gcamp_redux_for_RV\jgcamp7f_1FP.csv');
% gcamp7f_10fp = readtable('D:\Dropbox (HHMI)\janelia\meetings\illustrations\gcamp_redux_for_RV\jgcamp7f_10FP.csv');
% 
% gcamp7f_1fp_y = (gcamp7f_1fp.Y - min(gcamp7f_1fp.Y)) / (mean(gcamp7f_1fp.Y));
% gcamp7f_10fp_y = (gcamp7f_10fp.Y - min(gcamp7f_10fp.Y)) / (mean(gcamp7f_10fp.Y));
% 
% 
% figure, plot(new_1fp.X, new_1fp_y, gcamp6s_1fp.X, gcamp6s_1fp_y, gcamp7f_1fp.X, gcamp7f_1fp_y)
% legend({'500.311', 'gcamp6s', 'gcamp7f'})
% figure, plot(new_10fp.X, new_10fp_y, gcamp6s_10fp.X, gcamp6s_10fp_y, gcamp7f_10fp.X, gcamp7f_10fp_y)
% legend({'500.311', 'gcamp6s', 'gcamp7f'})

%% plot raw cell data launchpad

if isempty(whos('mutant'))
    load('Y:\GENIE_Pipeline\Analysis\Results\pile_all_GCaMP96uf_upto_20190814.mat', 'mutant')
end
currentMutant = mutant(contains({mutant.construct}, '500.456'));
figure, plot(squeeze(currentMutant.df_f_med(:,3,:)))
% plotRawCellData(currentMutant,5, 0);

[df_f_aligned] = align_responses(currentMutant);
figure, plot(squeeze(df_f_aligned(:,3,:)))
%% stim frame testing
% plateFolder = 'Y:\GENIE_Pipeline\GECI Imaging Data\20190806_GCaMP96uf_analyzed\P6a-20190722_GCaMP96uf\imaging\';
% wellsDir = dir(fullfile(plateFolder, '96Well*'));
% 
% stimFrames = [];
% for i = 1:length(wellsDir)
%     currentWellFolder = fullfile(plateFolder, wellsDir(i).name);
%     load(fullfile(currentWellFolder, 'ws_info_array.mat'), 'ws_info_array');
% 
%     % load stim times
%     % stimFrames = zeros(length(ws_info_array),1);
%     % for j = 1:length(ws_info_array)
%         [~, stimFrame] = min(abs(ws_info_array(1).ImageTime - ws_info_array(1).VPulseTime));
% %         stimFrames(j) = stimFrame;
%     % end
%     stimFrames = [stimFrames stimFrame];
%     %disp(stimFrame)
%     
% end
% 
% hist(stimFrames)