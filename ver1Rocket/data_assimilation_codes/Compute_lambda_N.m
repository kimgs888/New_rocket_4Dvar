function [lambda_N] = Compute_lambda_N (H,STEC,Xn)

% formulate the coveraince matrix

[R]=R_covariance_Matrix(STEC);

% calculate lambda_N as in Adjoint model enhanced plume reconstruction from tomographic
% remote sensing measurements 

lambda_N = (H'*R^-1)*(STEC-H*Xn);
