%ZDT1
function [Out] = ZDT2(X)
    %Dimension
    n = numel(X);
    
    %f1(x)
    f1 = X(1);
    %g(x)
    g = 1 + 9 / (n-1) * sum(X(2:end));
    %h(f1,g)
    h = 1 - (f1 / g)^2;
    %f2(x) = g(x)h(f1(x),g(x))
    f2 = g * h;
Out = [f1
       f2];