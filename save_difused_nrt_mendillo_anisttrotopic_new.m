 function [] = save_difused_nrt_mendillo_anisttrotopic_new(xi, yi, zi,Traje_time,Lat2,Long2,Heighti,hgtbgn,hgtend, Total_release_time,everyreleasex, end_of_release, start_of_release, temp_folder)

 % Find out whether the platform is pc or unix 
    if ispc 
         slashx = '\';
    elseif isunix 
        slashx = '/';
    end
    
 
RE = 6371e3; 
Heightd = 0:1:hgtend; 

%if strcmp(molcule,'H2O')
    HeightH20 = [0 250 350 450]; 
    DifusionH20 = [1 2 12 67];
    DifusionH20i = interp1( HeightH20,DifusionH20,Heightd,'linear','extrap') ;
%elseif strcmp(molcule,'H2')
    HeightH2 = [0 250 350 450]; 
    DifusionH2 = [1 6 39 210];
    DifusionH2i = interp1( HeightH2,DifusionH2,Heightd,'linear','extrap') ;
    
 

everyrelease = everyreleasex;
 
%diffusiontime = Total_release_time;

%maxrealse = round((end_of_release-start_of_release)*60/10);

% if  (diffusiontime*60)/10 >= maxrealse 
% number_of_realse = maxrealse;
% else 
% number_of_realse = round(diffusiontime*60/10);
% end

%tp = start_of_release;

for  tp = start_of_release:everyrelease:end_of_release

  [~,indxh] = min(abs(Traje_time-tp));
  
  

  
 point_of_relz = sphcart([RE+Heighti(indxh)*1e3,[Lat2(indxh),Long2(indxh)].*(pi/180)]);

 Pts_from_source = sphcart([RE+(yi(:).*1e3),[xi(:),zi(:)].*(pi/180)]);
 
 dist_from_source = (((Pts_from_source(:,1)-point_of_relz(:,1)).^2 + (Pts_from_source(:,2)-point_of_relz(:,2)).^2 ...
                       + (Pts_from_source(:,3)-point_of_relz(:,3)).^2));
                   
                  
                  pointofrelzhgt = RE+(Heighti(indxh)*1e3); % altitude or rocket (hr)
                  altpointfrmsource = RE+(yi(:).*1e3);       % altitude of points from source (hn)
                   
                   T = (pointofrelzhgt.^2 + altpointfrmsource.^2); 
                  % look at the derived equation for altitude of  mid point between
                  % rocket and n(r,t), this is the new defusion height ....
                  % We did this to mimic anistrotopic diffusion 
                   newdiffusion_height = (((2*T - dist_from_source)./4).^.5);                 
      

% generate new diffusion mid point between rocket and n(r,t) point
savetime=roundn(tp,-4); 


%molcule = 'H2O';


HGTd = ((newdiffusion_height-RE)/1000);



[D1] = get_difusion_every_point(Heightd,HGTd, DifusionH20i);

[D2] = get_difusion_every_point(Heightd,HGTd, DifusionH2i);

disp(['Difusion time in min: ', num2str(savetime)])

 
%tic 
%  for XI = 1:length(newdiffusion_height)
%  molcule = 'H2O'; 
%  [D1]=Diffusion_approx(newdiffusion_height(XI),hgtbgn,hgtend,molcule);
%   molcule = 'H2';  
%  [D2]=Diffusion_approx(newdiffusion_height(XI),hgtbgn,hgtend,molcule);
 tm = num2str(savetime);
 tm(tm=='.')='_';
 RE1 = ['D1X',tm,'.H20=D1;'];

 eval(RE1)
 RE2 = ['D1X',tm,'.H2=D2;']; 
 eval(RE2)
  
  
 RE3x =  ['save ',temp_folder,slashx, 'HGX',tm,' ', 'newdiffusion_height;'];
 eval(RE3x)

  RE3d =  ['save ',temp_folder,slashx , 'dst_frm_src_sqd',tm,' ', 'dist_from_source;'];
 eval(RE3d)
 
 RE3 =  ['save ',temp_folder,slashx , 'D1X',tm,'.mat ', 'D1X',tm,';'];
 eval(RE3)
 RE3 =  ['save ',temp_folder,slashx , 'D1X',tm,'.mat ', 'D1X',tm,';'];
 eval(RE3)

%  RE4 =  ['load D1X',tm,'.mat ', 'D1X',tm]
%  eval(RE4)
%  RE1 = ['D1X',tm,'.H20']
%  eval(RE1)
%tp = tp + everyrelease;
RE2 = ['clear D1X',tm, ' ','HGX',tm,' ' ];
eval(RE2)

end 







