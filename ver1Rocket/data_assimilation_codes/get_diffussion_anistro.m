function [D1,D2]=get_diffussion_anistro(temp_folder,savetime)
   if ispc 
         slashx = '\';
    elseif isunix 
        slashx = '/';
   end
    
tm = num2str(savetime);
tm(tm=='.')='_';


RE4 =  ['load ', temp_folder,slashx,'D1X',tm,'.mat ', 'D1X',tm,';'];

eval(RE4)

RE1 = ['D1=D1X',tm,'.H20;'];
eval(RE1)

RE2 = ['D2=D1X',tm,'.H2;']; 
eval(RE2)
D1 =  (D1*1e6)';  
D2 =  (D2*1e6)';

RE5 = ['clear D1X',tm];
eval(RE5)
