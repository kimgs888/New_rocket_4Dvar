
function [combined_elecX,combined_elec,All_lat,All_lon,All_HGT]=density (Glon,Glat,Gheight,year,dn,hry)
% change to IRI directory  to compute densities 

cd('IRI')
jm=0;
htec_max = 2000;
ivar = 1;
icontinue = 0; 

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

LZ = length( All_lat)*(length(All_lon))*(length(All_HGT));

LT= length(hry)*length(dn);

combined_elecX = zeros(LZ, LT);
cntx=1;
initial_pt=0;
if initial_pt ==0   
combined_elec=[];
All_lat=[];
All_lon=[];
All_HGT = [];


for KI = 1:length(dn)
    
for hour =hry 

 lnstpos = 0;
 lnendpos = 0;   
for Glonn=1:length(Glon)  
    
LNH = Glon{Glonn};
longbegin =LNH (1);
longstep = LNH(2)-LNH(1);
longend =  LNH(end); 
lnstpos= lnendpos+1;
lnendpos =  lnendpos + length(LNH);
ltstpos = 0;
ltendpos = 0;

for Glatn = 1:length(Glat)
    
LTH = Glat{Glatn};
latbegin = LTH(1);
latstep = LTH(2)-LTH(1);
latend =  LTH(end);
ltstpos=ltendpos+1;
ltendpos = ltendpos + length(LTH);
hghtstpos=0;
hghtendpos=0;

for GI = 1:length(Gheight)
    
HGT = Gheight{GI};
hgtbgn = HGT(1);
hgtstep = HGT(2)-HGT(1);
hgtend =HGT(end);
hghtstpos = hghtendpos+1;
hghtendpos=hghtendpos+length(HGT);
    


years =year(KI);
doys =  -1*dn(KI); % in IRI ,mmdd(or -ddd) multiply by negative -1 for doy
hours = hour;
fid=fopen('IRI2016input.txt','w+'); % open the input file and write to it
fprintf(fid,'%1.0f\n',jm); %'jmag(=0/1,geog/geom),lati/deg,long/deg'
fprintf(fid, '%5.4f\n',latbegin);
fprintf(fid, '%5.4f\n',latstep);
fprintf(fid, '%5.4f\n',latend);
fprintf(fid, '%5.4f\n',longbegin);
fprintf(fid, '%5.4f\n',longstep);
fprintf(fid, '%5.4f\n',longend);
fprintf(fid,'%1.0f\n', 1); %iut(=0/1,LT/UT)
fprintf(fid, '%5.4f\n', -1); % 'height/km'
fprintf(fid, '%5.0f\n',0); % 'output_option'
fprintf(fid, '%5.4f\n', htec_max); % 'upper height [km] for TEC integration (0 for no TEC)'
fprintf(fid, '%5.0f\n',ivar);%  %'variable? (1/2/../8 for height/lat/long/year/month/', 'day/day of year/hour)'
fprintf(fid, '%5.4f\t %5.4f\t %5.4f\n',hgtbgn,hgtend,hgtstep); %'begin, end, and stepsize for the selected variable'
fprintf(fid, '%5.0f\n', 0); % 'Enter 0 to use standard or 1 to enter your own'   
%fprintf(fid, '%5.0f\n', icontinue);  % flag to read to terminate the IRI.
%hours = [11 12 13];
fprintf(fid, '%1.0f\n', length(years));  % flag to read to terminate the IRI.
fprintf(fid, '%1.0f\n', length(doys));  % flag to read to terminate the IRI.
fprintf(fid, '%1.0f\n', length(hours));  % flag to read to terminate the IRI.
f1 ='%5.4f\t ';
flst = '%5.4f\n';
Fy = repmat (f1,1,length(years));
Fy(end-1)='n';
Sy = ['fprintf(fid,','''',Fy,''',[',num2str(years),']);'];
eval(Sy)
 
Fd = repmat (f1,1,length(doys));
Fd(end-1)='n';
Sd = ['fprintf(fid,','''',Fd,''',[',num2str(doys),']);'];
eval(Sd)
 
Fh = repmat (f1,1,length(hours));
Fh(end-1)='n';
Sh = ['fprintf(fid,','''',Fh,''',[',num2str(hours),']);'];
eval(Sh)
 
fclose(fid); 
delete('fort.20')
delete('fort.15')

if ispc 
 y = dos ('IRI2016_nick.exe  &');
elseif isunix 
system('./IRI2016_nick')
end 
while 1
if exist('fort.20','file')
%TASKKILL [/S system [/U username [/P [password]]]]
%         { [/FI filter] [/PID processid | /IM imagename] } [/T] [/F]

%Description:
%This tool is used to terminate tasks by process id (PID) or image name.
%/f forces the the kill process 
% nul sets the screen output as nul i.e not shown 
if ispc 
[R]=dos('taskkill /im cmd.exe /f  > nul');  
end 

fid1 = fopen('fort.15');


out = textscan (fid1, '%*f %f');
fclose all;

lt = latbegin:latstep:latend; 
ln = longbegin:longstep:longend;
HGT =  hgtbgn:hgtstep:hgtend;

[YY1,ZZ1,XX1]=meshgrid(lt,HGT,ln);

elec = out{1}*1e6;
Sz=size(XX1); 

elecnew = reshape(elec, size(XX1));
% put all the variable latitudes longs and heights into one vector

if initial_pt == 0
    
    if (GI==1)&&(Glonn==1)
    All_lat=[All_lat,lt];
    end 

    if (GI==1)&&(Glatn==1)
    All_lon=[All_lon,ln];
    end 

    if (Glatn==1)&&(Glonn==1)
    All_HGT = [All_HGT,HGT];
    end 
end

% also store the combined electron density in pages 
combined_elec(hghtstpos:hghtendpos, ltstpos:ltendpos,lnstpos:lnendpos) = elecnew;

break
end 
end 
end 
end
end

combined_elecX(:,cntx)= combined_elec(:);
cntx=cntx+1;
initial_pt = 1;
end
end
end
% change to IRI directory  to compute densities 
cd('..')
