function [] = saveRxLat_RxLon_PRN_temp(diffusiontimex,RxLat,RxLon,PRN,Time,simulation_analysis,time_to_collect_rays,time_sample,locationx)
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
RE1 = ['RxLt_Ln_prn_tm',prex,tm,'_',colect,'_',time_samplex,'.RxLat=RxLat;'];

eval(RE1);
RE2 = ['RxLt_Ln_prn_tm',prex,tm,'_',colect,'_',time_samplex,'.RxLon=RxLon;']; 
eval(RE2);

RE2 = ['RxLt_Ln_prn_tm',prex,tm,'_',colect,'_',time_samplex,'.PRN=PRN;']; 
eval(RE2);
RE2 = ['RxLt_Ln_prn_tm',prex,tm,'_',colect,'_',time_samplex,'.Time=Time;']; 
eval(RE2);
cd (locationx)
RE3 =  ['save RxLt_Ln_prn_tm',prex,tm,'_',colect,'_',time_samplex,'.mat ', 'RxLt_Ln_prn_tm',prex,tm,'_',colect,'_',time_samplex,';'];

eval(RE3);
cd .. 

