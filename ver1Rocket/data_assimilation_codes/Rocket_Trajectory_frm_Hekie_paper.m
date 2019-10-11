function [Lat2,Long2,Heighti,Traje_time]=Rocket_Trajectory_frm_Hekie_paper(Long1,Lat1,bearig)
time = [0 1 2.5 4.5 6 6.5];
Traje_time = 0:.01:7;
Height = [0 5 90 200 290 320];
Horizontal =[0 0 90 300 590 720];
Heighti = interp1(time,Height,Traje_time,'PCHIP','extrap') ;
Horizontali = interp1(time,Horizontal,Traje_time,'PCHIP','extrap'); 
[Lat2,Long2]= trajectoryin_long_lat(Long1,Lat1,bearig,Horizontali);