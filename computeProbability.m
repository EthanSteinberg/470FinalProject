function [f,g] = computeProbability(x, network, tree)
% Compute the probability (and possibly the gradient) of a tree
% given a network and x parameter values

% f is the probability
% g is the gradient of said probability
% The gradient is only calculated when necessary

calllib('libnetworkprob', 'changeParams', network, x);

if nargout > 1 % gradient required
    [f, g] = calllib('libnetworkprob', 'computeProbability', network, tree, x);
else
    f = calllib('libnetworkprob', 'computeProbability', network, tree, []);
end