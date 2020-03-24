function [cMap] = getColorMap(nLines)
%GETCOLORMAP Based on number of lines to plot, find optimal colormap
% INPUTS
%       nLines: number of lines to plot
% OUTPUTS
%       cMap: [nLines x 3] matrix of colors (R,G,B)
% To plot, use 
%              for i = 1:nLines, plot(x,y, 'color', cMap(i,:), end

% assume max distribution over colors

vec = .1:1/nLines:1;


cMap = [vec(randperm(nLines)); vec(randperm(nLines)); vec(randperm(nLines))]';
end

