 function [] = save_difused_nrt_mendillo_anisttrotopic(latbegin,latend,longbegin,longend,hgtbgn,hgtend,longstep,latstep,...
                                            hgtstep,Long1,Lat1,bearig,diffusiontimex,everyreleasex)

[xi yi zi] = meshgrid(latbegin:latstep:latend,hgtbgn:hgtstep:hgtend,longbegin:longstep:longend);

time = [0 1 2.5 4.5 6 6.5];
xii = 0:.01:7;
Height = [0 5 90 200 290 320];
Horizontal =[0 0 90 300 590 720];
Heighti = interp1(time,Height,xii,'cubic','extrap') ;
Horizontali = interp1(time,Horizontal,xii,'cubic','extrap'); 
[Lat2,Long2 ]= trajectoryin_long_lat(Long1,Lat1,bearig,Horizontali);
everyrelease = everyreleasex;

diffusiontime = diffusiontimex;

maxrealse = round((6.5-2.5)*60/10);

if  (diffusiontime*60)/10 >= maxrealse 
number_of_realse = maxrealse;
else 
number_of_realse = round(diffusiontime*60/10);
end

tp = 2.5;
for tpx = 1:number_of_realse
 
  [val,indxh] = min(abs(xii-tp)); 
  
 point_of_relz = sphcart([64E5+Heighti(indxh)*1e3,[Lat2(indxh),Long2(indxh)].*(pi/180)]);
 Pts_from_source = sphcart([64E5+yi(:).*1e3,[xi(:),zi(:)].*(pi/180)]);
 
 dist_from_source = (((Pts_from_source(:,1)-point_of_relz(:,1)).^2 + (Pts_from_source(:,2)-point_of_relz(:,2)).^2 ...
                       + (Pts_from_source(:,3)-point_of_relz(:,3)).^2));
                   
                  
                  pointofrelzhgt = 64E5+Heighti(indxh)*1e3; % altitude or rocket (hr)
                  altpointfrmsource = 64E5+yi(:).*1e3;       % altitude of points from source (hn)
                   
                   T = (pointofrelzhgt.^2 + altpointfrmsource.^2); 
          
                   newdiffusion_height = (((2*T - dist_from_source)./4).^.5);                 
      
 
 
 
% generate new diffusion mid point between rocket and n(r,t) point
savetime=roundn(tp,-4)

 for XI = 1:length(newdiffusion_height)
 molcule = 'H2O'; 
 [D1]=Diffusion_approx(newdiffusion_height(XI),hgtbgn,hgtend,molcule);
  molcule = 'H2';  
 [D2]=Diffusion_approx(newdiffusion_height(XI),hgtbgn,hgtend,molcule);
 tm = num2str(savetime);
 tm(tm=='.')='_';
 RE1 = ['D1X',tm,'.H20(XI)=D1;'];

 eval(RE1)
 RE2 = ['D1X',tm,'.H2(XI)=D2;']; 
 eval(RE2)
 
 end
 cd ('F:\Data assimilati\codes\Reviewers\Reveiwcode_withchanged_diffusion_coffient\diffusion_data_anistro2')
 RE3 =  ['save D1X',tm,'.mat ', 'D1X',tm];
 eval(RE3)
%  RE4 =  ['load D1X',tm,'.mat ', 'D1X',tm]
%  eval(RE4)
%  RE1 = ['D1X',tm,'.H20']
%  eval(RE1)
tp = tp + everyrelease;
 end 






