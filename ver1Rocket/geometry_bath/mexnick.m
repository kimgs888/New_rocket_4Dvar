% This function creates the raytrace Mex function called in matlab 
% during raytrace
% you need to have all the header files in the same folder or specify there
% location. 
% The output is a raytrace which is then called in ratracer. 
% By Nicholas Ssessanga 
fclose all 
clear functions                             % Clear all MEX functions
WorkDir  = cd;                              % Copy current directory
MakePath = fileparts(which('make'));        % Locate MEX directory
cd(MakePath) 
mex -c matrices.c -output matrices
mex -c trace.c -output trace 
mex  -O mexraytrace.c -output raytrace trace.obj matrices.obj