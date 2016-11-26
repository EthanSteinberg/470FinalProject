% Determine optimal branch lengths for the migration sample problem.

loadLibrary();

network = createMigrationNetwork();

treesAndWeights = tdfread('weights','\t');

weights = [];

for i=1:size(treesAndWeights.TREE, 1)
    trees(i) = generateTreeFromNewick(treesAndWeights.TREE(i,:));
    weights(i) = treesAndWeights.WEIGHT(i);
end

bestProb = 0;
bestLengths = [];

for i=1:10
    x0 = rand(10, 1);
    [negProb, lengths] = quasiNewton(x0, @(x) computeNegativeTotalProbability(x, network, trees, weights));

    if negProb < bestProb
        bestProb = negProb;
        bestLengths = lengths;
    end
end
    
calllib('libnetworkprob', 'freeNetworkBuffer', network.buffer);
 
for i=1:size(treesAndWeights.TREE, 1)
    calllib('libnetworkprob', 'freeTreeBuffer', trees(i).buffer);
end