function [lambda_N,R] = Compute_lambda_Nx (H,STEC,Xn,bError,Time,increase_sd)

% formulate the coveraince matrix

[R]=R_covariance_Matrix2x(STEC,bError,Time,increase_sd);


% calculate lambda_N as in Adjoint model enhanced plume reconstruction from tomographic
%remote sensing measurements 

lambda_N = (H'*R^-1)*(STEC-H*Xn);
