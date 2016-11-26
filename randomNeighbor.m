function result = randomNeighbor(x, low, high)
% Find a random neighbor by randomly changing one of the parameters.
% Low is the lower bounds for the parameters.
% High is the upper bounds.

indexToChange = randi(length(x));

maxDelta = 0.1;

maxNewLength = min(high(indexToChange), x(indexToChange) + maxDelta / 2);

minNewLength = max(low(indexToChange), x(indexToChange) - maxDelta / 2);

result = x;
result(indexToChange) = ...
    minNewLength + rand() * (maxNewLength - minNewLength);
