function make()
% MatLAB function to compile MEX interfaces for MINDI
%    make()
%
% Notes:
%   Ensure that makefileunix or makefilewin have the correct paths and
% compiler setting for your system.  You may need to run 'mex -setup'
% from the matlab prompt prior to calling this function. Libraries for
% the fortran to C converter are included in bin. Unix libraries for
% 64 and 32 bit systems are located in f2c/linux32 and f2c/linux64, the
% appropriate libraries should be copied from there to bin.
%
% See also MEX
% --------------------------------------------------------------------------
clear functions                              % Clear all MEX functions
WorkDir  = cd;                               % Copy current directory
MakePath = fileparts(which('make'));         % Locate MEX directory
cd(MakePath);                                % Change to MEX directory


try
% --------------------------------------------------------------------------
% Windows code

   if ispc
      !nmake -f makefilewin clean
      !nmake -f makefilewin all
      Ext = ['.',mexext]; if strcmp(Ext,'.dll'), Ext = ''; else Ext = '.dll'; end
      fprintf('making raytrace.dll\n' );
      eval(['mex -O mexraytrace.c -output raytrace',Ext, 'trace.c matrices.c']);

%       fprintf('making iri95.dll\n' );
%       eval(['mex -O mexiri95.c -output iri95',Ext,' -IC:\korea work\temography\midas\nikiz',...
%          '../../bin/cira86.obj ../../bin/irif13.obj ../../bin/iris13.obj ',...
%          '../../bin/irit13.obj  ../../bin/bccf2c.lib']);
% 
% %      fprintf('making igrf.dll\n' );
% %      eval(['mex -O mexigrf.c -output igrf',Ext,' -I../../include ',...
% %         '../../bin/igrf.obj ../../bin/bccf2c.lib']);
% 
%       fprintf('making igrf.dll\n' );
%       eval(['mex -O mexigrf11.c -output igrf',Ext,' -I../../include ',...
%          '../../bin/igrf11.obj ../../bin/bccf2c.lib']);
% 
%       fprintf('making weimer.dll\n' );
%       eval(['mex -O mexweimer.c -output weimer',Ext,' -I../../include ',...
%          '../../bin/weimer.obj ../../bin/bccf2c.lib']);


% --------------------------------------------------------------------------
% Unix code
   else
      !make -fmakefileunix clean
      !make -fmakefileunix all

      fprintf('making raytrace.%s\n',mexext );
      mex -O mexraytrace.c -output raytrace -I../../include ...
         ../../bin/trace.o ../../bin/matrices.o

      fprintf('making iri95.%s\n',mexext );
      mex -O mexiri95.c -output iri95 -I../../include  ...
         ../../bin/cira86.o ../../bin/irif13.o ../../bin/iris13.o ...
         ../../bin/irit13.o ../../bin/libf2c.a

      fprintf('making igrf.%s\n',mexext );
      mex -O mexigrf11.c -output igrf -I../../include  ...
         ../../bin/igrf11.o ../../bin/libf2c.a

      fprintf('making weimer.%s\n',mexext );
      mex -O mexweimer.c -output weimer -I../../include  ...
         ../../bin/weimer.o ../../bin/libf2c.a

   end
end
% Change back to start directory
cd( WorkDir );

