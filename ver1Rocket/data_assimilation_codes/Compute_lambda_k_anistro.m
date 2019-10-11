function [lambda_k_1,gradient_langrange,R] = Compute_lambda_k_anistro(H,STEC,XK,bError,time,increase_sd,lambda_k,end_of_release,...
                                                                                start_of_release,diffusiontimex,everyreleasex,EoH2O,EoH2,temp_folder)

% By Nicholas Ssessanga while at Chungnam National University
% 2019-july-02

% Evaluate the adjoint
% evaluate  = 
%             d[Xk+1]^T
% lamda_k =   -------         
%             d[Xk]
[Adjnt] = Adjoint_nrt_anistro(end_of_release,start_of_release,diffusiontimex,everyreleasex,EoH2O,EoH2,XK,temp_folder)  ;                                                         
                                                             
% get the R 
[R]= R_covariance_Matrix2x(STEC,bError,time,increase_sd);

% evaluate  = 
%             d[Xk+1]^T
% lamda_k =   -------   lambda_k+1 + Hk^T * Rk^-1 (Hk*Xk-Yk ),k=N-1,N-2,0      
%             d[Xk]


gradient_langrange = Adjnt'*lambda_k;
HTR = (H'*R^-1);
HX_b = (H*XK-STEC);
lambda_k_1 = gradient_langrange - HTR *HX_b;
