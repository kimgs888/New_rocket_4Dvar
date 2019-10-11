function []=get_index_ofnear_pt_in_vect(A,B)
size( A(:))
size(reshape(B,1,[]))

TMP = bsxfun(@(x,y) abs(x-y), B(:), reshape(A,1,[]));
size(TMP)

[D, idxB] = min(TMP,[],2) ;
 Result = B(idxB); 
 size(Result)
 pause
% TFDiffLessThen3 = D < 3 