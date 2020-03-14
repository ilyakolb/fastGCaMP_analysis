%% fast gcamp overall analysis
% RUN FROM DESK COMPUTER
%%
% 10.641 GCaMP6s   : pGP-SIV-(SalI)-syn-IRES-nls-mCherry-WPRE-GCaMP3 K78H T302L R303P D380Y T381R S383T R392G.10.641
% 10.693 GCaMP6f   : pGP-SIV-(SalI)-syn-IRES-nls-mCherry-WPRE-GCaMP3 T302L R303P A317E D380Y T381R S383T R392G.10.693
% 10.921 GCaMP7f  : pGP-SIV-(SalI)-syn-IRES-nls-mCherry-WPRE-GCaMP3 T302L R303P A317L D380Y.10.921
% close all
clc

plotHighlights = 1;

% construct/control parameters to plot
% xToPlot = 'x1_fp';
% yToPlot = 'rise_1_fp';
% sizeToPlot = 'decay_1_fp';
% colorToPlot = 'norm_f0';
xToPlot = 'x1_fp';
yToPlot = 'rise_1_fp'; %'rise_1_fp'; %'norm_f0'
sizeToPlot = 'decay_1_fp';
colorToPlot = 'norm_f0';

base = 'Z:/';
if ismac
    error('NO NOT RUN THIS CODE ON THE OLD MAC ANYMORE!')
	base = fullfile('/Volumes', 'genie');
end



% 6th round with fixed jgcamp7f control
% good = readtable(fullfile(base,'\GECIScreenData\Analysis\data_week_20200303_GCaMP96uf_analyzed_GCaMP96uf.xlsx'));

% ALL after week 2 of 6th round (updated ilastik params)
good = readtable(fullfile(base,'\GECIScreenData\Analysis\data_all_20200308_GCaMP96uf.xlsx'));

% ALL after week 2 of 6th round (old ilastik params)
% good = readtable(fullfile(base,'\GECIScreenData\Analysis\data_all_20200212_GCaMP96uf.xlsx'));

% 6th round up to 20200211 (week 2)
% good = readtable(fullfile(base,'\GECIScreenData\Analysis\data_week_20200211_GCaMP96uf_raw_GCaMP96uf.xlsx'));

% 6th round up to 20200205 (week 1)
% good = readtable(fullfile(base,'\GECIScreenData\Analysis\data_week_20200205_GCaMP96uf_analyzed_GCaMP96uf.xlsx'));

% All rounds w/ fixed pixel bug
% good = readtable(fullfile(base,'\GECIScreenData\Analysis\data_all_20200128_GCaMP96uf.xlsx'));

% 5th round with fixed bug + ilastik
% good = readtable(fullfile(base,'\GECIScreenData\Analysis\data_week_20191209_GCaMP96uf_raw_GCaMP96uf.xlsx'));

% all rounds w/ bug
% good = readtable(fullfile(base,'\GECIScreenData\Analysis\withBug\data_all_20191122_GCaMP96uf.xlsx'));

% Yan's combined 5th round (20191125)
% good = readtable(fullfile(base,'\GECIScreenData\Analysis\data_week_20191125_GCaMP96uf_analyzed_GCaMP96uf.xlsx'));

% Yan's combined 5th round (20191118)
% good = readtable(fullfile(base,'\GECIScreenData\Analysis\data_week_20191118_GCaMP96uf_analyzed_GCaMP96uf.xlsx'));

% Yan's combined 4th round with fixed pixel bug
% good = readtable(fullfile(base,'\GECIScreenData\Analysis\data_all_4thround.xlsx'));


% weekly XCaMP comparison (12/2/19)
% good = readtable(fullfile(base,'\GECIScreenData\Analysis\data_week_20191125_XCaMP_analyzed_GCaMP96uf.xlsx'));

% weekly XCaMP comparison
% good = readtable(fullfile(base,'\GECIScreenData\Analysis\data_week_20191005_GCaMP96uf_raw_GCaMP96uf.xlsx'));

% Yan's new round of screening (11/22/19)
% good = readtable(fullfile(base,'\GECIScreenData\Analysis\data_week_20191118_GCaMP96uf_analyzed_GCaMP96uf.xlsx'));

% Yan's new round of screening (11/27/19)
% good = readtable(fullfile(base,'\GECIScreenData\Analysis\data_week_20191125_GCaMP96uf_raw_GCaMP96uf.xlsx'));


% Abhi's sensor plate
% good = readtable(fullfile(base,'\GECIScreenData\Analysis\data_week_20191112_GCaMP96uf_Rig2_raw_mngGECO.xlsx')); % 8 AP stims
% good = readtable(fullfile(base,'\GECIScreenData\Analysis\data_week_20191006_GCaMP96uf_analyzed_GCaMP96uf.xlsx'));


% XCaMP data
% good = readtable(fullfile(base,'GENIE_Pipeline/Analysis/Results/data_week_20190904_GCaMP96uf_analyzed_GCaMP96uf.xlsx'));


% standardize column names when loading xlsx or txt files
switch width(good)
    case 52
        load varNames_4AP
    case 84
        load varNames_8AP
    otherwise
        error('Wrong number of columns')
end

good.Properties.VariableNames(~startsWith(good.Properties.VariableNames, 'Var'))= varNames;

% only include if at least X wells
good = good(good.replicate_number > 1, :);

controlNames = {'10.641', '10.693', '500.333', '10.921', '500.456'};
plotControls = 1;
nToSort = 5;

tau_on_variants = good.(sizeToPlot);
min_tau_on = min(tau_on_variants);
max_tau_on = max(tau_on_variants);
tau_on_variants_norm = 2*((tau_on_variants - min_tau_on) ./ max_tau_on);


% variants (from good.mat)
cn = good.construct; %.construct
x = good.(xToPlot);
y = good.(yToPlot);
cs = 10+200*tau_on_variants_norm;
cc = log(good.(colorToPlot)); % 2-log(good.(colorToPlot));

% save identifying x, y coords for interactive plotting
% save('plotted_xycoords.mat', 'x','y')

figure('position', [167          50        1194         908])
% plot variants
scatter(x,y,cs,cc,'filled')
% set(gca,'xscale','log')
% set(gca,'yscale','log')
text(x,y,strcat({'   '}, cn),'fontsize',7); %overlay construct names
grid on

% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% legend
title('fast GCaMP variants')
xlabel(xToPlot, 'Interpreter', 'none')
ylabel(yToPlot, 'Interpreter', 'none')
c = colorbar;
c.Label.String = colorToPlot;

% plot controls
if plotControls
    % plot highlights
    control_idx = cellfun(@(x) find(contains(cn, x)), controlNames, 'UniformOutput', false);
    found_controls_idx = vertcat(control_idx{:});
    if isempty(found_controls_idx)
        warning('Highlighted constructs not found!')
    end
    highlight_good = good(found_controls_idx,:);
    x_controls = highlight_good.(xToPlot);
    y_controls = highlight_good.(yToPlot);
    hold on, scatter(x_controls, y_controls, 20, [1 0 0])
end

%% sort good constructs
good_sorted = sortrows(good, xToPlot, 'descend'); % pos-going
head(good_sorted, nToSort)
% good_sorted = good(good.first_assay_date>2.019e7,:);
highlights = table2cell(good(contains(good.construct, {'500.666'}),'construct'));
% highlights = table2cell(good(startsWith(good.first_assay_date, '2020') ,'construct'));

% hits = table2cell(good_sorted(1:nToSort,1))'

%% highlight constructs w/ red circle
if plotHighlights
    % plot highlights
    highlight_idx = cellfun(@(x) find(contains(cn, x)), highlights, 'UniformOutput', false);
    found_highlight_idx = vertcat(highlight_idx{:});
    if isempty(found_highlight_idx)
        warning('Highlighted constructs not found!')
    end
    highlight_good = good(found_highlight_idx,:);
    x_highlight = highlight_good.(xToPlot);
    y_highlight = highlight_good.(yToPlot);
    scatter(x_highlight, y_highlight, 20, [1 0 0])
end