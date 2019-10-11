function [Lat2,Long2 ] = trajectoryin_long_lat(Long1,Lat1,bearig,d)
% This functions computes the longitude and latitude of the trajectory of
% the Rocket on ground  using the distance measured on ground (d in km). 
% The user must provide the starting Log1 Lat1,brearing and d(Km)
% all data should be provided in degrees. 
% Reference: http://www.movable-type.co.uk/scripts/latlong.html
% Destination point given distance and bearing from start point

bearig = bearig*pi/180;
R = 6371; % raduis of the Earth. 
Long1 = Long1*pi/180;
Lat1 =Lat1*pi/180;
% Now find the longitude and latitude  
Lat2 = asin(sin(Lat1)*cos(d/R) + cos(Lat1)*sin(d/R)*cos(bearig));
Long2 = Long1 + atan2(sin(bearig)*sin(d/R)*cos(Lat1), cos(d/R)-sin(Lat1)*sin(Lat2));
% include the starting point as part of the data
Lat2 =Lat2.*(180/pi);
Long2 = Long2.*(180/pi);