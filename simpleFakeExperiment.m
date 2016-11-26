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

[newtonProb, newtonLength] = learningPlot(@(x0) quasiNewton(x0, func), x, 10);

title('10 Runs Of Quasi-Newton Method');

print('../finalProjectPaper/newtonSimple','-dpng')

%% Gradient
figure('Visible','off');
hold on;

[gradProb, gradLength] = learningPlot(@(x0) gradientDescent(x0, func, low, high, 100), x, 10);

title('10 Runs Of Gradient Descent Method');

print('../finalProjectPaper/gradientSimple','-dpng')

%% Hill climbing
figure('Visible','off');
hold on;

[hillProb, hillLength] = learningPlot(@(x0) hillClimbing(x0, func, low, high, 900), x, 10);

title('10 Runs Of Hill Climbing Method');

print('../finalProjectPaper/hillClimbSimple','-dpng')

%% Cleanup

calllib('libnetworkprob', 'freeNetworkBuffer', network.buffer);
calllib('libnetworkprob', 'freeTreeBuffer', trees(1).buffer);