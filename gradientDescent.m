function [fx, x, values, funccount] = gradientDescent(x0, f, low, high, maxIters)
% Peform the gradient descent algorithm to find the minimum of f
% x0 is the starting point
% f should provide both the value and gradient
% maxIters is the number of iterations to perform.

% Returns the final final, the final x, and the list of points

funccount = 1:maxIters;

x = x0;
iter = 1;

[values(iter), gx] = f(x);

alpha = 0.1;

gamma = 0.9;

v = 0;

while iter < maxIters
    iter = iter + 1;

    v = gamma * v + alpha * gx';
    x = x - v;
    x = max(low, x);
    x = min(high, x);
    
    [values(iter), gx] = f(x);
end

fx = values(maxIters);