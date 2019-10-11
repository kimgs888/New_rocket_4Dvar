
function [sat]= calcsatposnk (Elmat, Azmat, stationLong,stationLat,  satalitude)
% This function computes the Trasimitter (Tx) location at a specific altitude 
% Here we are assume the Tx is located at "satalitude" set by the user. 
% We use the location of the of the receiver together with the elevation and azimuth 
% of the satellite to trace where Tx would be located. 
% By Nicholas Ssessanga while at Chungnam National University
% 2019-july-02

Re = 6371e3; 
hipp = satalitude; 
ru = stationLat;
rln = stationLong; 
chng2rad = pi/180;
chng2deg = 180/pi;
Elev = Elmat*chng2rad;
Azim = Azmat*chng2rad;

rcu = ru*chng2rad;
rlon =rln*chng2rad;

phipp =((pi/2)-Elev- asin((Re/(Re+hipp)).*cos(Elev)));

ipplat = asin(sin(rcu).*cos(phipp) + cos(rcu).*sin(phipp).*cos(Azim));

ipplong = rlon + asin((sin(phipp).*sin(Azim))./cos(ipplat));

sat.lat = ipplat*chng2deg;
sat.long= ipplong*chng2deg;







