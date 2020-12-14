%% calculate F0 of 'best performer' weeks
% save CSVs for Prism export to D:\ufgcamp_paper_data\culture-APdata-csv\f0

clearvars -except mutants

% 10.641 must be first in list for normalization to work properly
control = '10.641';
hits = {'10.921','500.456', '500.686', '500.688', '500.712', '500.543', '500.707', '500.455', '10.1473', '10.1513', '10.1561', '538.1', '538.2', '538.3'};

base = 'Z:/';

if isempty(whos('mutants'))
    % all bests
    mutant1 = load(fullfile(base,'GECIScreenData\Analysis\pile_week_GCaMP96uf_upto_20201131_GCaMP96uf_analyzed.mat'), 'mutant');
    
    % all bests except 500.686
    mutant2 = load(fullfile(base,'GECIScreenData\Analysis\pile_week_GCaMP96uf_upto_20200310_GCaMP96uf_analyzed.mat'), 'mutant');
    
    mutants = [mutant1 mutant2]; % and mutant2
end

% get control F0s
f0_control = zeros(length(mutants),1);

for i = 1:length(f0_control)
    construct_list = {mutants(i).mutant.construct}';
    m = mutants(i).mutant(contains(construct_list, control));
    f0_control(i) = mean(m.f0);
end

for i = 1:length(hits)
    
    f0_array = [];
    % f0_control = 0;
    for j = 1:length(mutants) % cycle through weeks
        construct_list = {mutants(j).mutant.construct}';
        m = mutants(j).mutant(contains(construct_list, hits{i}));
        if ~isempty(m)
            f0 = m.f0';
            
            % normalize by in-week GCaMP6s fluorescence
            f0_array = [f0_array; f0 / f0_control(j)];
        end
        
    end

    writematrix(f0_array, fullfile('D:\ufgcamp_paper_data\culture-APdata-csv\f0\',[hits{i} '.csv']))
end

