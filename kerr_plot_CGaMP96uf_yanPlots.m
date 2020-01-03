%% fast gcamp overall analysis
%% 20190928: generating figures for Yan
% RUN FROM DESK COMPUTER
% highlights plotted in red circles
% controls plotted in blue circles
% kerr plot taken from interactive_kerr_plot.m
% using CSV generated from Mac scripts
% Y:\GENIE_Pipeline\Analysis\Results\data_all_20190509_GCaMP96uf.txt

% table description
% x1_fp, x3_fp, x10_fp, x160fp: presumably df/f?

% 10.641 GCaMP6s   : pGP-SIV-(SalI)-syn-IRES-nls-mCherry-WPRE-GCaMP3 K78H T302L R303P D380Y T381R S383T R392G.10.641
% 10.693 GCaMP6f   : pGP-SIV-(SalI)-syn-IRES-nls-mCherry-WPRE-GCaMP3 T302L R303P A317E D380Y T381R S383T R392G.10.693
% 10.921 GCaMP7f??  : pGP-SIV-(SalI)-syn-IRES-nls-mCherry-WPRE-GCaMP3 T302L R303P A317L D380Y.10.921
% close all
clc

plotHighlights = 0;
plotControls = 1;
nKerrPlots = 3;

xToPlot_series = {'x3_fp', 'x3_fp', 'x3_fp'}; % df/f
yToPlot_series = {'rise_3_fp', 'decay_3_fp', 'x160_fp'}; % 1/time
xLabels = {'DF/F (3AP)', 'DF/F (3AP)', 'DF/F (3AP)'};
yLabels = {'1/rise (3AP)', '1/decay (3AP)', 'DF/F (160AP)'};

controlNames = {'10.641', '10.921'};

oneOverY = [1 1 0]; % 1 if y axis should be inverse

sizeToPlot = 'norm_f0'; % size: F0
colorToPlot = 'x3_fp_p'; % inv log df/f0 p-val

base = 'Z:/';
if ismac
    base = fullfile('/Volumes', 'genie');
end

% Yan's all rounds with fixed pixel bug
good = readtable(fullfile(base,'\GECIScreenData\Analysis\data_all_20191211_GCaMP96uf.xlsx'));

% sort to just 2nd round (500.2-500.303 + controls)
% good = good([1:3 67:328], :);

% sort to everything else
good = good([1:66 328:end], :);


% Yan's combined 4th round with fixed pixel bug
% good = readtable(fullfile(base,'\GECIScreenData\Analysis\data_all_4thround.xlsx'));

% all data
% good = readtable(fullfile(base,'GENIE_Pipeline/Analysis/Results/data_all_20190915_GCaMP96uf.xlsx'));

% XCaMP data
% good = readtable(fullfile(base,'GENIE_Pipeline/Analysis/Results/data_week_20190904_GCaMP96uf_analyzed_GCaMP96uf.xlsx'));


% standardize column names when loading xlsx or txt files
load varNames_4AP.mat
good.Properties.VariableNames(~startsWith(good.Properties.VariableNames, 'Var'))= varNames;

%% condition good
% only include if at least X wells
good = good(good.replicate_number > 1, :);

% #1 410dot1 to 410dot31
% good_410 = good(contains(good.construct, '410.'), :);
% good_410 = good_410(1:17, :);
% good_410 = [good_410; good(contains(good.construct, controlNames), :)];
% good = good_410;
%

% #2 up to 500.374 
% good_500374 = good(1:424, :);
% good_500374 = [good_500374; good(contains(good.construct, controlNames), :)];
% good = good_500374;

% #3 all data -- good = good

% #4 controls+best+XCaMPs
% good_bests = good(contains(good.construct, {'500.110', '500.543', '10.641', '10.921', '500.176', '500.333', '500.456', '500.640', '538.'}), :);
% good_bests = [good_bests; good(endsWith(good.construct, {'410.6', '410.80'}), :)]; % exact matches
% good = good_bests;

% these are bad:
good = good(~contains(good.construct, {'497.32', '497.440', '410.9', 'TE', 'none', '376.13'}), :);



for nPlot = 1:nKerrPlots
    
    xToPlot = xToPlot_series{nPlot};
    yToPlot = yToPlot_series{nPlot};
    
    size_variable = good.(sizeToPlot);
    min_size_variable = min(size_variable);
    max_size_variable = max(size_variable);
    size_variable_norm = (size_variable - min_size_variable) ./ max_size_variable;
    
    
    % variants (from good.mat)
    cn = good.construct; %.construct
    x = good.(xToPlot);
    y = good.(yToPlot);
    if oneOverY(nPlot)
        y = 1./y;
    end
    cs = 10+200*size_variable_norm;
    cc = good.(colorToPlot);
    % color for p values screws up for controls b/c there are so many
%     cc(contains(good.construct, {'TE', '10.921', '10.693', '500.333'})) = 1;
    cc = -1*log10(cc);
    % save identifying x, y coords for interactive plotting
    % save('plotted_xycoords.mat', 'x','y')
    
    figure('position', [2210         176         847         717])
    % plot variants
    scatter(x,y,cs,cc,'filled')
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    text(x,y,strcat({'   '}, cn),'fontsize',7); %overlay construct names
    grid on
    
    % set(gca, 'YScale', 'log')
    % set(gca, 'XScale', 'log')
    % legend
    % title('fast GCaMP variants (20190806)')
    xlabel(xLabels{nPlot}, 'Interpreter', 'none')
    ylabel(yLabels{nPlot}, 'Interpreter', 'none')
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
        % hold on, scatter(x_controls, y_controls, 20, [1 0 0])
    end
    
    %% sort good constructs
    good_sorted = sortrows(good, xToPlot, 'descend'); % pos-going
    head(good_sorted, 5)
    % good_sorted = good(good.first_assay_date>2.019e7,:);
    highlights = table2cell(good(contains(good.construct, {'500.456', '500.333'}),'construct'));
    
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
        hold on,scatter(x_highlight, y_highlight, 20, [1 0 0])
    end
    
end