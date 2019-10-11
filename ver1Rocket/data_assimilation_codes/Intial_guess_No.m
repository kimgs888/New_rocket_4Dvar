function [RXo] = Intial_guess_No(H,b,U,No_EOF)
nmber_EOF =No_EOF;

U(:,1) = U(:,1).*-1;
Q = U;
M = H*Q;
a = pinv(M)*b;
RXo = Q*a;
while min(RXo)< 0
nmber_EOF=nmber_EOF-1;
Q = U(:,1:nmber_EOF);
HN = H'*H;
M = HN*Q;
BN = (H'*b);
a = pinv(M)*BN;
RXo = Q*a;
end 
disp(['number of EOFs utilized: ', num2str(nmber_EOF)])
end 
