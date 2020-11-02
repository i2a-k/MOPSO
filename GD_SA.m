function [GD,SA] = GD_SA(di,S_PF,dm,PF)
    
    Di = [];
    Dm = [];
    K = size(di,1)/S_PF;
    Z = 1;
    
    for i=1:K:size(di,1)
        Di = cat(1,Di,min(di(Z:K,1)));
        dm_index = di(Z:K,1) == min(di(Z:K,1));
        Dm = cat(1,Dm,dm(dm_index,1));
        Z = Z + K;
    end
    
    GD = ( sqrt(sum(Di.^2)) ) / rank(PF);
    
    SA = ( sum(Dm) + sum(abs(Di-mean(Di))) ) / ( sum(Dm) + rank(PF)*mean(Di) );