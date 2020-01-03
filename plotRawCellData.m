%% find raw data associated with hit
% hitMutant: element from `mutants` struct array
% maxNumToShow: maximum number of figures to show to prevent slowdowns
% loading param files from Y:\GENIE_Pipeline\GECI Imaging Data\20190611_GCaMP96uf_analyzed\P1a-20190527_GCaMP96uf\imaging\96Well38-D02
% 20190913: returns stim frame numbers -- ws only?
function plotRawCellData(hitMutant, varargin)

%% user params
fs = 200; % imaging sampling rate
launchFiji = 0; % 1 to start Fiji window with movie (buggy...)
subplotsForImage = [1 2 3 6 7 8 11 12 13 16 17 18];
subplotsForAPs = [4 9 14 19];
subplotsForSingleTraces = [5 10 15 20];
imagingDataDir = 'Z:\GECIScreenData\GECI_Imaging_Data';
APnumString = '001FP'; % string to identify which AP num tiff stack to show
noRedSignal = 1; %set to 1 to NOT show red signal on image (e.g. if no mCherry)
%%
allPlateDir = dir(fullfile(imagingDataDir, '*_GCaMP96uf_analyzed*\P*')); % all plate directories
allPlateDir = allPlateDir([allPlateDir.isdir]);
allPlateNames = {allPlateDir.name}';

nWells = length(hitMutant.plate);

if nargin > 1
    maxNumToShow = varargin{1};
    if nWells > maxNumToShow
        warning(['Only showing first ' num2str(maxNumToShow) ' wells of ' num2str(nWells)])
        nWells = maxNumToShow;
    end
end

if nargin > 2
    launchFiji = varargin{2};
end

for i = 1:nWells
    currentPlate = hitMutant.plate{end-i};
    currentWell = hitMutant.well{end-i};
    
    currentPlateDir = allPlateDir(strcmp(allPlateNames, currentPlate));
    assert(~isempty(currentPlateDir), ['Plate ' currentPlate ' not found!'])
    
    currentWellDir = dir(fullfile(currentPlateDir.folder, currentPlateDir.name, 'imaging', ['96Well*' currentWell]));
    assert(~isempty(currentWellDir), ['Well ' currentWell ' of plate ' currentPlate ' not found!'])
    
    currentWellFolder = fullfile(currentWellDir.folder, currentWellDir.name);
    pArray = load(fullfile(currentWellFolder, 'para_array_cherry.mat'));
    
    % load stim times
    % load(fullfile(currentWellFolder, 'ws_info_array.mat'), 'ws_info_array');
    
    
%     stimFrames = zeros(length(ws_info_array),1);
%     for j = 1:length(ws_info_array)
%         [~, stimFrame] = min(abs(ws_info_array(1).ImageTime - ws_info_array(1).VPulseTime));
%         stimFrames(j) = stimFrame;
%     end
    segmentation = load(fullfile(currentWellFolder, 'segmentation_cherry.mat'));
    
    
    figure('Name', currentWellFolder);
    nTime = length(pArray.para_array(1).df_f);
    tArray = (1:nTime)/fs;
    nCells = size(pArray.para_array, 2);
    
    
    
    redSignal = zeros(size(segmentation.GCaMPbase2));
    if ~noRedSignal
        redSignal = segmentation.mCherry;
    end

    overlay = create_overlay(segmentation.GCaMPbase2, redSignal);
    imSize = size(overlay);
    
    % mark cells on overlay
    for j = 1:length(pArray.cell_list)
        BW = false(imSize(1:2));
        BW(pArray.cell_list(j).pixel_list) = true;
        B = bwperim(BW, 4);
        shadow = circshift(B, [1 1]);
        shadow(1, :) = false;
        shadow(:, 1) = false;
        % overlay(shadow, :) = 0;
        tempOverlay = overlay(:,:,1);
        tempOverlay(B) = 1;
        tempOverlay(shadow) = 0;
        overlay(:,:,1) = tempOverlay;
        
        tempOverlay = overlay(:,:,2);
        tempOverlay(B) = 1;
        tempOverlay(shadow) = 0;
        overlay(:,:,2) = tempOverlay;
        
        tempOverlay = overlay(:,:,3);
        tempOverlay(B) = 0;
        tempOverlay(shadow) = 0;
        overlay(:,:,3) = tempOverlay;
    end
    
    % draw well image
    subplot(4,5,subplotsForImage)
    imagesc(overlay)
    axis image
    set(gca, 'XTick', [], 'YTick', []);
    title([hitMutant.construct ' | ' num2str(nCells) ' cells '])
    % plot AP responses
    ylabels = {'1AP', '3AP', '10AP', '160AP'};
    for k = 1:4 % num APs
        df_f = [pArray.para_array(k,:).df_f];
        subplot(4,5,subplotsForAPs(k))
        hold on,plot(tArray, df_f, 'color', [.6 .6 .6])
        plot(tArray, mean(df_f,2), 'k-')
        ylabel(ylabels{k});
        
        % plot 1AP traces separately
        if k == 1
            subplot(4,5,subplotsForSingleTraces)
            plot(tArray,df_f+2*repmat(1:size(df_f,2), size(df_f,1),1), 'color', [.6 .6 .6])
            title('1AP raw')
            set(gca, 'YTick', [])
            axis tight; box off; 
        end
    
    end
    

    % show animation w/ Fiji
    if launchFiji
        FPfileDir = dir(fullfile(currentWellFolder, ['P*' APnumString '*.tif']));
        assert(~isempty(FPfileDir), [APnumString ' TIF file not found in ' currentWellFolder '!'])
        system(['C:\Users\kolbi\Downloads\Fiji.app\ImageJ-win64 "' ...
            fullfile(FPfileDir(1).folder, FPfileDir(1).name) '" -port1 &']);
    end
    
end

end