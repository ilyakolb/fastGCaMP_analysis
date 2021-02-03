function [normError] = normalized_error(A,B, dim)
%NORMALIZED_ERROR error propagation resulting from A/B
% calculates standard error
% see https://en.wikipedia.org/wiki/Propagation_of_uncertainty

normError = nanmean(A, dim)./nanmean(B, dim).*...
    sqrt( (nanstd(A,0,dim)./nanmean(A,dim)).^2 + (nanstd(B,0,dim)./nanmean(B,dim)).^2  );

normError = normError ./ sqrt(size(A, 2));
end

