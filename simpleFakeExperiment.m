% Perform experiments on a simple fake topology

%% Setup
loadLibrary();

network = createSimpleNetwork();
trees = createAllTrees('A':'C');

x = [0.5, 1, 0.5, 0.5, 1, 4, 1, .7];

low = zeros(size(x)) + 0.0001;
high = inf(size(x));
high(end) = 0.9999;

weights = computeExpectedWeights(x, network, trees);

func = @(x) computeNegativeTotalProbability(x, network, trees, weights);


%% Quasi newton
figure('Visible','off');
hold on;

newtonProbs = zeros(10, 1);
newtonLengths = zeros(10, length(x));

for i=1:10
    x0 = rand(size(x));
    [negProb, lengths, values, funccount] = quasiNewton(x0, func);
    plot(funccount, -values);
    newtonProbs(i) = negProb;
    newtonLengths(i,:) = lengths;
end

xlabel('Number of evaluations');
ylabel('P(geneTrees|speciesNetwork)');
title('10 Runs Of Quasi-Newton Method');

print('../finalProjectPaper/newtonSimple','-dpng')

bestNewtonLengths = newtonLengths(newtonProbs == min(newtonProbs),:);

%% Gradient
figure('Visible','off');
hold on;

gradProbs = zeros(10, 1);
gradLengths = zeros(10, length(x));

for i=1:10
    x0 = rand(size(x));
    [negProb, lengths, values] = gradientDescent(x0, func, low, high, 100);
    plot(1:100, -values);
    gradProbs(i) = negProb;
    gradLengths(i,:) = lengths;
end

xlabel('Number of evaluations');
ylabel('P(geneTrees|speciesNetwork)');
title('10 Runs Of Gradient Descent Method');

print('../finalProjectPaper/gradientSimple','-dpng')

bestGradLengths = gradLengths(gradProbs == min(gradProbs),:);

%% Hill climbing
figure('Visible','off');
hold on;

hillProbs = zeros(10, 1);
hillLengths = zeros(10, length(x));

for i=1:10
    x0 = rand(size(x));
    [negProb, lengths, values] = hillClimbing(x0, func, low, high, 900);
    plot(1:900, -values);
    hillProbs(i) = negProb;
    hillLengths(i,:) = lengths;
end

xlabel('Number of evaluations');
ylabel('P(geneTrees|speciesNetwork)');
title('10 Runs Of Hill Climbing Method');

print('../finalProjectPaper/hillClimbSimple','-dpng')

bestHillLengths = hillLengths(hillProbs == min(hillProbs),:);

calllib('libnetworkprob', 'freeNetworkBuffer', network.buffer);
calllib('libnetworkprob', 'freeTreeBuffer', trees(1).buffer);