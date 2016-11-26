function [bestProb, bestLength] = learningPlot(func, x, numIters)
% Create a learning plot for the function func for random starting points.
% Returns the best probability and the best lengths.

allProbs = zeros(numIters, 1);
allLengths = zeros(numIters, length(x));

for i=1:numIters
    x0 = rand(size(x));
    [negProb, lengths, values, funccount] = func(x0);
    plot(funccount, -values);
    allProbs(i) = negProb;
    allLengths(i,:) = lengths;
end

xlabel('Number of evaluations');
ylabel('P(geneTrees|speciesNetwork)');

bestProb = min(allProbs);
bestLength = allLengths(find(allProbs == bestProb), :);

end