function []= Create_raytrace_mex(working_dir) 
% This function creates the raytrace Mex function called in matlab 
% during raytrace
% you need to have all the header files in the same folder or specify there
% location. 
% The output is a raytrace which is then called in ratracer. 
% By Nicholas Ssessanga 

fclose all ;
clear functions                             % Clear all MEX functions

 if ispc 
             slashx = '\';
        elseif isunix 
            slashx = '/';
 end 
   
if exist([working_dir,slashx,'geometry_bath'], 'dir')==7
cd([working_dir,slashx,'geometry_bath'])
mex -c matrices.c -output matrices
mex -c trace.c -output trace 
mex  -O mexraytrace.c -output raytrace trace.o matrices.o
cd(working_dir)

else 
    error('Geometry folder is missing, please double check')
end 