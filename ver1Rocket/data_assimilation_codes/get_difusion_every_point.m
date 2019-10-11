function [Difus_HGT] = get_difusion_every_point(interp_difus_hght,height,difusion)

TMP = bsxfun(@(x,y) abs(x-y), interp_difus_hght(:), reshape(height,1,[]));

[~,index]=min(TMP,[],1); 
Difus_HGT = difusion(index);