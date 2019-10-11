function [All_lat,All_lon,All_HGT]=gridvectors(Glon,Glat,Gheight)


All_lat=[];
All_lon=[];
All_HGT = [];


for I = 1:length(Glat)
All_lat = [All_lat;Glat{I}'];
end 
for I = 1:length(Glon)
All_lon= [All_lon;Glon{I}'];
end 
for I = 1:length(Gheight)
All_HGT= [All_HGT;Gheight{I}'];
end 

