function [modelDensityPtr] = set_min_density(modelDensityPtr) 
dens_min = 0.01*modelDensityPtr;
den_min_indx = dens_min<1.0e6;
modelDensityPtr(den_min_indx) = 1.0e6;