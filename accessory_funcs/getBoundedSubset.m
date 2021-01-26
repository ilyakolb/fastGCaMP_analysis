function [subset] = getBoundedSubset(g, series, startIdx, endIdx)
% GETBOUNDEDSUBSET: return 'good' array matching a particular series,
% bounded by startIndx and endIdx
% g: good table
% series (str): series should end with . i.e. '500.'
% startIdx: (int) index with which to start/stop matching
% example: getBoundedSubset(g, '500.', 16, 27) returns table with
% constructs 500.16 - 500.27

% first match the series

subset = g(startsWith(g.construct, series), :);

idx_array = startIdx:endIdx;
idx_array_str = arrayfun(@(x) ['.' num2str(x)], idx_array, 'UniformOutput', 0);

% then match the startIdx:endIdx array
subset = subset(endsWith(subset.construct, idx_array_str), :);

end
