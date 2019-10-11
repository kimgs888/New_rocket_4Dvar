function Tx = Extract_gps_position (navfile, time,PRN)
% This code reads in the GPS navigation files (.n) or igs files 
% and extracts the position of the required satellite in view. 
% Define the time when the satellite is need.
% by Nicholas 25 sept 2014 
% See also: read_rinexn.m, get_satpos and xyz2wgs




% call read_rinexn to extract the ephoc(eph) of each satellite
rinexn = navfile;
eph = read_rinexn(rinexn);
% pass the eph to  "get_satpos" fuction and extract the satellite position
% in ECEF (cartesian coordinates)
sv = PRN;

t= time;
flag =0;
 satpos = get_satpos(t,sv,eph,flag);  
  x = satpos(1);
  y = satpos(2);
  z = satpos(3);
  
 
  % change the satellite cartesian coordinates to the reciever frame (WGS84)
  % the out put is in the form (lon, lat, alt)
%  S = [t, x,y,z];
%  R = xyz2wgs(S);
%  t = R(1);
%  lon = R(2);
%  lat = R(3);
 %alt = R(4);
  Tx = [x,y,z]; % this is the required output to be used in the Ray generation coord

 
  