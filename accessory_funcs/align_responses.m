function [df_f_med_aligned] = align_responses(currentMutant, nStims, debugFlag)
%ALIGN_RESPONSES
%   align controlMutant.df_f_med. Only for new constructs with ws_info_array
%   shifts all traces to make sure they start at 1 sec mark
%   crawls through processed image directories
%   VARIABLES:
%     currentMutant: entry from piled data (e.g.
%                 pile_all_GCaMP96uf_upto_20190904_GCaMP96uf_analyzed.mat)
%                 nStims: number of AP stims presumed (e.g. 1,3,10,160)
%                 (e.g. 4 or 8)
%     debugFlag: set to 1 to plot a figure of 1AP traces before and after
%                 alignment

warning('FIX align_responses!')
trueStimIndex = 201; % make the plots centered on 1-second mark
imagingDataDir = 'Z:\GECIScreenData\GECI_Imaging_Data';
% APnumString = '001FP';
allPlateDir = dir(fullfile(imagingDataDir, '*_GCaMP96uf_*\P*')); % all plate directories
allPlateDir = allPlateDir([allPlateDir.isdir]);
allPlateNames = {allPlateDir.name}';

stimFrames = zeros(length(currentMutant.plate),nStims);
for i = 1:length(currentMutant.plate)
    
    currentPlate = currentMutant.plate{i};
    currentWell = currentMutant.well{i};
    currentPlateDir = allPlateDir(strcmp(allPlateNames, currentPlate));
    if(~isempty(currentPlateDir))
        
        currentWellDir = dir(fullfile(currentPlateDir.folder, currentPlateDir.name, 'imaging', ['96Well*' currentWell]));
        if (~isempty(currentWellDir))
            
            
            currentWellFolder = fullfile(currentWellDir.folder, currentWellDir.name);
            
            try
                load(fullfile(currentWellFolder, 'ws_info_array.mat'), 'ws_info_array');
                allStimFrames = zeros(1,4);
                for k = 1:length(ws_info_array)
                    [~,allStimFrames(k)] = min(abs(ws_info_array(k).ImageTime - ws_info_array(k).VPulseTime(1)));
                end
            catch e
                warning('ws_info_array not found! Must be an old plate')
                allStimFrames = ones(1,nStims) * trueStimIndex;
            end
            
        else
            warning(['Well ' currentWell ' of plate ' currentPlate ' not found!'])
            allStimFrames = ones(1,nStims) * trueStimIndex;
        end
        
    else
        warning(['Plate ' currentPlate ' not found!'])
        allStimFrames = ones(1,nStims) * trueStimIndex;
    end
    
    stimFrames(i,:) = allStimFrames;
    
end

df_f_med = currentMutant.df_f_med;
stimDelta = stimFrames - trueStimIndex;

df_f_med_aligned = zeros(size(df_f_med));


% figure,hist(stimFrames)
for i = 1:size(df_f_med,2) % cycle through stim number (1,3,10,160)
    for j = 1:size(df_f_med,3) % cycle through n replicates
        current_df_f_med = df_f_med(:,i,j);
        df_f_med_aligned(:,i,j) = circshift(current_df_f_med, -1*stimDelta(j,i), 1);
        % replace with nans
        if stimDelta(j,i) < 0
            df_f_med_aligned(1:abs(stimDelta(j,i)),i,j) = nan;
        elseif stimDelta(j,i) > 0
            df_f_med_aligned(end-stimDelta(j,i):end,i,j) = nan;
        end
    end
end

% for debugging
if debugFlag
    figure, subplot(1,2,1)
    title('not aligned')
    plot(squeeze(currentMutant.df_f_med(:,1,:)))
    subplot(1,2,2)
    plot(squeeze(df_f_med_aligned(:,1,:)))
    title('aligned')
end