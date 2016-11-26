function [prob, result, values, funccount] = quasinewton(x0, f)
% Run the quasi newton algorithm from matlab

values = [];
funccount = [];

options = optimoptions(@fminunc,'Algorithm','quasi-newton', 'OutputFcn', @outfun);
options = optimoptions(options,'SpecifyObjectiveGradient',true);


[result,prob,~,output]  = fminunc(f,x0, options);

calls = output.funcCount;


function stop = outfun(x, optimValues, state)
stop = false;
values(length(values) + 1) = optimValues.fval;
funccount(length(funccount) + 1) = optimValues.funccount;
end

end
