function weights = computeExpectedWeights(x, network, trees)
% Compute the expected weights for a network and a bunch of trees.
% The expected weights are simply the probability of that tree given the network.


for i = 1:length(trees)
    weights(i) = computeProbability(x, network, trees(i));
end

end