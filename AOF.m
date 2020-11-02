%Aggregation Function
function [Out] = AOF(X,k)
    r = rand(k,1);
    w = zeros(k,1);
    for i=1:k
        w(i,1) = r(i,1)/sum(r);
    end
    
    P = ZDT1(X);
    F = w.*P;
Out = sum(F);