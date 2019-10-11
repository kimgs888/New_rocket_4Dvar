% This code aproximates the difussion 
% based on results from medillo1975
% Table 3: Spherical Diffusion Parameters at Various Ionospheric Heights
function [DFSion]=Diffusion_approx(heightpnt,hgtbgn,hgtend,molcule)
%molcule = 'H2O';
if strcmp(molcule,'H2O')
    Height = [0 250 350 450]; 
    Difusion = [1 2 12 67];

elseif strcmp(molcule,'H2')
    Height = [0 250 350 450]; 
    Difusion = [1 6 39 210];
end 

Heightd = hgtbgn:1:hgtend; 
Difusioni = interp1( Height,Difusion,Heightd,'linear','extrap') ;
RE = 6371e3;
[~,I]=min(abs(Heightd-(heightpnt-RE)/1000));
DFSion =  Difusioni(I);
