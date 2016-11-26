function [f,g] = computeNegativeTotalProbability(x, network, trees, weights)

% Compute the negative of thetotal probability of a bunch of trees, 
% weighting each tree according to weights.
% The total is the product of P(tree | network) ^ weight.

% Also computes the corresponding gradients when necessary.

treeProb = zeros(1, length(trees));
treeDeriv = zeros(length(trees), length(x));

for i = 1:length(trees)
    if nargout > 1
        [treeProb(i), treeDeriv(i,:)] = computeProbability(x, network, trees(i));
    else
        treeProb(i) = computeProbability(x, network, trees(i));
    end
end

f = -prod(treeProb .^ weights);

g = zeros(1, length(x));

if nargout > 1
    for i = 1:length(trees)
        term = ones(1, length(x));

        for j = 1:length(trees)
            if i == j
                term = (term .* treeDeriv(j,:)) * weights(j) * treeProb(j) ^ (weights(j) - 1);
            else
                term = term * treeProb(j) ^ weights(j);
            end
        end

        g = g + term;
    end
    
    g = -g';
end

