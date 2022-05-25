function [p]=bernoull(x,eta);
% [p] = bernoull(x,eta)
%
% Computes Bernoulli distribution of x for "natural parameter" eta.
% The mean m of a Bernoulli distributions relates to eta as,
% m = exp(eta)/(1+exp(eta));

e = exp(eta);

p = ones(size(e));  

indx = find(~isinf(e));

p(indx) = exp(eta(indx).*x - log(1+e(indx)));
