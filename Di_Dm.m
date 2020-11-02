function [di,dm] = Di_Dm(PF,Z_PF)
    di = [];
    dm = [];
    
    for a=1:size(PF,1)
        for b=1:size(Z_PF,1)
            %Eu distance
            di = cat(1,di,sqrt(sum((PF(a,:)-Z_PF(b,:)).^2)));
            %dm
            dm = cat(1,dm,sum(abs((PF(a,:)-Z_PF(b,:)))));
        end
    end