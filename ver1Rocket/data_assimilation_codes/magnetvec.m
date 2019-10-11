function [magcord] = magnetvec(All_lat,All_HGT,All_lon)

[latggv,altgv,longgv]=meshgrid(All_lat,All_HGT,All_lon);
latstep = .1;
longstep=.1;
np = length(longgv(:));
magcord =zeros(np,3);
latggv = latggv(:);
longgv = longgv(:);
altgv = altgv (:);
 for i=1:np
in_lat = latggv(i);%*GetRadToDeg() % input latitude
in_lon = longgv(i);%*GetRadToDeg() % input longitude
in_alt = altgv(i);

time=12;
year=2012;
month=07;
dom=23;
 [lat,lon, alt] = plotgeomagnick(in_lat,in_lat,latstep,in_lon ,in_lon ,longstep,in_alt,in_alt,time,year,month,dom);
 
lat = mean(lat);
lon = mean(lon);
alt = mean(alt);

 % convert lon to 0 to 360
 if lon<0
 lon = lon + 360;
 end
 % convert to radians
 latrad = lat*pi/180; 
 lonrad = lon*pi/180; 

 magcord (i,1:3)=  [latrad ,lonrad,alt];
end 