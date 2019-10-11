function [] = saveHSTEC_temp(diffusiontimex,H,STEC,bError,simulation_analysis,time_to_collect_rays,time_sample,locationx)
if simulation_analysis ==0
prex ='O'; 
elseif simulation_analysis == 1
prex ='S'; 
end 
time_samplex=num2str(time_sample);
time_samplex(time_samplex=='.')='_';

colect=num2str(time_to_collect_rays);
colect(colect=='.')='_';

savetime=roundn(diffusiontimex,-2);
tm = num2str(savetime);
tm(tm=='.')='_';
RE1 = ['HSTEC',prex,tm,'_',colect,'_',time_samplex,'.H=H;'];

eval(RE1);
RE2 = ['HSTEC',prex,tm,'_',colect,'_',time_samplex,'.STEC=STEC;']; 
eval(RE2);

RE2 = ['HSTEC',prex,tm,'_',colect,'_',time_samplex,'.bError=bError;']; 
eval(RE2);
cd (locationx)
RE3 =  ['save HSTEC',prex,tm,'_',colect,'_',time_samplex,'.mat ', 'HSTEC',prex,tm,'_',colect,'_',time_samplex,';'];

eval(RE3);
cd .. 
