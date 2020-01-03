% 12/4/19
% quick script to get full plate directories from partial names
% used to rerun everything
% run from Z:\ilya\code\fastGCaMP_analysis\accessory_funcs
clearvars
load allPlates

allOldPlates = allPlates(~contains(allPlates, '2019'));

allListedPlates = dir('Y:\GENIE_Pipeline\GECI Imaging Data\*\*');
allListedPlateNames = {allListedPlates.name}';
imgPlate = {};
for i = 1:length(allOldPlates)
    foundIdx = find(contains(allListedPlateNames, allOldPlates{i}),1);
    plateFullDir = allListedPlates(foundIdx);
    topLvlFolder = plateFullDir.folder;
    topLvlFolder = strsplit(topLvlFolder, '\');
    
    imgPlate{end+1} = topLvlFolder{4};
end

imgPlate = imgPlate';
imgPlate = unique(imgPlate);
% transfectionDates = 