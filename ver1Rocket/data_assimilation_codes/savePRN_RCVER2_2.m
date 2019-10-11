function [] = savePRN_RCVER2_2(diffusiontimex,prnall,Rxltall,Rxlonall,obstimeall,simulation_analysis,time_to_collect_rays,time_sample,locationx)
if simulation_analysis ==0
prex ='O'; 
elseif simulation_analysis == 1
prex ='S'; 
end 

savetime=roundn(diffusiontimex,-2);
tm = num2str(savetime);

tm(tm=='.')='_';
colect=num2str(time_to_collect_rays);
colect(colect=='.')='_';

time_samplex=num2str(time_sample);
time_samplex(time_samplex=='.')='_';
RE1 = ['PRN_RCVER',prex,tm,'_',colect,'_',time_samplex,'.prnall=prnall;'];

eval(RE1);
RE1 = ['PRN_RCVER',prex,tm,'_',colect,'_',time_samplex,'.Rxltall=Rxltall;'];

eval(RE1);
RE1 = ['PRN_RCVER',prex,tm,'_',colect,'_',time_samplex,'.Rxlonall=Rxlonall;'];

eval(RE1);

RE1 = ['PRN_RCVER',prex,tm,'_',colect,'_',time_samplex,'.obstimeall=obstimeall;'];

eval(RE1);
cd (locationx)
RE3 =  ['save PRN_RCVER',prex,tm,'_',colect,'_',time_samplex,'.mat ', 'PRN_RCVER',prex,tm,'_',colect,'_',time_samplex,';'];


eval(RE3);
cd .. 
