function [normError] = normalized_error(A,B, dim)
%NORMALIZED_ERROR error propagation resulting from A/B
% calculates standard error
% see https://en.wikipedia.org/wiki/Propagation_of_uncertainty

normError = mean(A, dim)./mean(B, dim).*...
    sqrt( (std(A,0,dim)./mean(A,dim)).^2 + (std(B,0,dim)./mean(B,dim)).^2  );

normError = normError ./ sqrt(size(A, 2));
end

