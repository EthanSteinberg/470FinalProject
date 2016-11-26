function [fx, x, values, funccount] = hillClimbing(x0, f, low, high, maxIters)
% Performs the hill climbing algorithm to find then minimum of the function f
% x0 is the starting point
% f should provide both the value and gradient
% maxIters is the number of iterations to perform.

% Returns the final final, the final x, and the list of points

funccount = 1:maxIters;

x = x0;
iter = 1;

values(iter) = f(x);

while iter < maxIters
    iter = iter + 1;

    next = randomNeighbor(x, low, high);
    fnext = f(next);

    if (fnext < values(iter - 1))
        % We have found an improvement
        x = next;
        values(iter) = fnext;
    else
        values(iter) = values(iter - 1);
    end
end

fx = values(maxIters);