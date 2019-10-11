function [lt,lg,Rr,fx0,fX,df,X]=creat_H_matrix3 (S,Time,year,doy)

% NeedH = [];
% H = 0;
% v = 0;
% n_file = dir(S);
CMNfile= dir(S);

% %figure 
% %hold on 
% for file = 1:length(n_file)
%   
%  CMNfile = n_file(file).name;  % now get the file name 
%  
%  %Call this file and extract the Rx, and latitude and logitude to
% %generate the Ray in constracting the the H matrix 


[info]= read_CMN_extract_Rx_pos (CMNfile,Time);
% if flag ==1  %check to see if the time was available, if not,  do not do any thing
     
stlat = info.satlat;   % get the satellite latitudes 
stlon = info.satlon;   % get the satellite longitudes 
stalt = info.Satalt;   % get the satellite altitudes
Rxlat = info.RxLat;    % get the Receiver  latitude
Rxlon = info.RxLon;    % get the Receiver  Longitude
Rxalt = info.RxAlt;    % get the Receiver  altitudes
%rxlt = repmat(64E5, size(stlat));



Rx = sphcart([64E5+Rxalt',[Rxlat',Rxlon']*pi/180]);
Tx = sphcart([(64E5+stalt)',[stlat',stlon']*pi/180]);

Grdlatst = 20;
Grdlatend = 50;
Grdlongst = 120;
Grdlongend = 150;
Grdhgtst = 64E5+90E3;
Grdhgtend = 64E5+1500E3;

hstep = 50E3;
latstep =2; 
longstep =2; 
% Creat the  H matrix
   Grid  = setgrid('Rad',Grdhgtst:hstep:Grdhgtend,'Lat',Grdlatst:latstep:Grdlatend,'Lon',Grdlongst:longstep:Grdlongend);
   Ray   = gpsray(Rx,Tx,Grid);
   [H] = raytracer(Grid, Ray);
   size(H)
   pause
  
 
  
  

% extract the cartesian cordinates and use them to obtain xo intial guess.
%   gx = Grid.X(:); % : was j
%   gy = Grid.Y(:); % : was j
%   gz =Grid.Z(:); % : was j
%   longbegin = Grdlatst ;
%   longend = Grdlatend ;
%   longstep = longstep;
  
%   cf = 0;
%   
%     for longpt = Grdlongst:longstep:Grdlongend
%       
%       for latpt = Grdlatst:latstep:Grdlatend
%          hour = Time;
%          dn = doy;
%        hgtbegin = (Grdhgtst-64E5)/1E3;
%         hgtend= (Grdhgtend-64E5)/1E3;
%         hgtstep=hstep/1E3;  
%        [outcolne] = profile_for_H_matrix2(latpt,longpt,year, hour,dn,hgtbegin,hgtend,hgtstep);
%        if cf ==0
%         xintial.xo = outcolne;
%         cf = 1;
%         else
%          xintial.xo(length(xintial.xo)+1:length(xintial.xo)+length(outcolne)) = outcolne;
%         end 
%         
%       end 
%   
%     end
dn1 = doy;
hour = Time;
month = 07;
[pp] = calibratefof2(Grdlatst, Grdlatend, Grdlongst, Grdlongend, Grdhgtst, Grdhgtend,hstep,latstep,longstep,dn1,hour,month,year);
  
%   for h =Grdhgtst:hstep:Grdhgtend
%       hgt=(h-64E5)/1E3;
%       for lt = Grdlongst:latstep:Grdlongend
%           latpt = lt;
%           hour = Time;
%           dn = doy;
%            
%         [outcolne] = profile_for_H_matrix2(latpt, year, hour,dn,hgt,longbegin,longend,longstep);
%        if cf ==0
%         xintial.xo = outcolne;
%         cf = 1;
%         else
%          xintial.xo(length(xintial.xo)+1:length(xintial.xo)+length(outcolne)) = outcolne;
%         end 
%         
%       end 
   
 %end 
% length(xintial.xo);
% length(Grid.X(:));



%for K = 1:29:30
K =100;  % Now set up the number iterations 
b = info.stec';
A = H; 
%pp = xintial.xo;
gg = abs(sqrt(pp*1e6/1.24e10));
x0 = 1.24e10*(gg).^2;
size(x0);
size(A)
clear pp 
%   [X] = kaczmarz(A,b,K,x0);
% figure
% plot(X, 'g')
% title('kr')
% hold on 
% figure
% plot(x0, 'b')
% title('kr')
 %[X] = MART2(A,b,K,x0);
%  plot(b*10E16,'*'); 
%  figure
%  plot(A*x0,'*r')

 b = b*10E16;
 %[X] = kaczmarz(A,b,K,x0);
  [X] = MART1(A,b,K,x0);
% size(X) 
% size(A)
% size(x0)

%  figure 
%  plot(A*X,'*k')
%  figure
%  plot((A*X-b), '*g')
%  figure
%  plot((A*X-b)./b, '*b')
 
% [X] = kaczmarz(A,b,K,x0);
 
%  figure 
%  plot(A*X,'*k')
%  figure
%  plot((A*X-b), '*g')
%  figure
%  plot((A*X-b)./b, '*b')
%  
 
%  [X] = kaczmarz(A,b,K,x0);
%  X2 = X;
%  figure 
%  plot(A*X2,'ok')
%  figure
%  plot((A*X2-b), 'og')
%   figure
%  plot((A*X2-b)./b, '+b')
%  
 

% figure
% plot(X, 'r')
% title('mart')
% hold on 
% figure
% plot(x0, 'b')
% title('mart')
%  figure
% scatter3(Grid.X(:),Grid.Y(:),Grid.Z(:),20,x0,'filled');
% colorbar

[lat,long, r]= cart2sph(Grid.X(:),Grid.Y(:),Grid.Z(:));
% figure
% scatter3(Grid.X(:),Grid.Y(:),Grid.Z(:),20,X,'filled');
% colorbar
% figure
% scatter3(long*180/pi,lat*180/pi,(r-64E5)/1E3,10,X,'filled');
% title('Electron Density (MART)')
%  ylabel ('Longitude (deg)', 'Fontsize', 12, 'FontWeight', 'bold','Fontname', 'Times New Roman')
%  xlabel ('Latitude (deg)', 'Fontsize', 12, 'FontWeight', 'bold','Fontname', 'Times New Roman')
%  zlabel ('Altitude (km)', 'Fontsize', 12, 'FontWeight', 'bold','Fontname', 'Times New Roman')
%  set (gca,'Fontsize', 12, 'FontWeight', 'bold','Fontname', 'Times New Roman')
%  zlim([90 1800])
% h = colorbar;
% view([-82 26])
% title(h,'Ne','Fontsize', 12, 'FontWeight', 'bold','Fontname', 'Times New Roman')
% figure 
% scatter3(long*180/pi,lat*180/pi,r-64E5,10,x0,'filled');
% title('orign')
% title('Electron Density (IRI-2012)')
%  ylabel ('Longitude (deg)', 'Fontsize', 12, 'FontWeight', 'bold','Fontname', 'Times New Roman')
%  xlabel ('Latitude (deg)', 'Fontsize', 12, 'FontWeight', 'bold','Fontname', 'Times New Roman')
%  zlabel ('Altitude (km)', 'Fontsize', 12, 'FontWeight', 'bold','Fontname', 'Times New Roman')
%  set (gca,'Fontsize', 12, 'FontWeight', 'bold','Fontname', 'Times New Roman')
% view([-82 26])
% colorbar

% figure 
% scatter3(long*180/pi,lat*180/pi,r-64E5,10,X-x0,'filled');
% title('orign')
% title('diff (IRI-2012)')
%  ylabel ('Longitude (deg)', 'Fontsize', 12, 'FontWeight', 'bold','Fontname', 'Times New Roman')
%  xlabel ('Latitude (deg)', 'Fontsize', 12, 'FontWeight', 'bold','Fontname', 'Times New Roman')
%  zlabel ('Altitude (km)', 'Fontsize', 12, 'FontWeight', 'bold','Fontname', 'Times New Roman')
%  set (gca,'Fontsize', 12, 'FontWeight', 'bold','Fontname', 'Times New Roman')
% view([-82 26])
% colorbar


% figure
fX = ((X/1.24e10).^0.5);
%  plot (fX(100:300),'.r')
%  hold on 
fx0 = ((x0/1.24e10).^0.5);
% plot(fx0(100:300),'.k')
df = fX-fx0;
% figure
% plot (abs(df),'.g')
% figure
% plot (df(100:300),'.b')


lt = long*180/pi;
lg = lat*180/pi;
Rr = (r-64E5)/1E3; 
% xind = find (abs(lg-round(140))<=eps(round(140))& abs(lt-round(25))<=eps(round(25)));
% 
% 
% 
% % lg = lg(indx1);
% % lt = lg(indx1);
% 
% %xind = find ((lg == 140)& (lt==30));
% figure
% plot(fx0(xind),Rr(xind), '.r')
% hold on 
% plot(fX(xind),Rr(xind),'.k')
% hold off
% figure 
% scatter3(long*180/pi,lat*180/pi,r,20,X-x0,'filled');
% title('diff')
% colorbar
% figure
%   x =long*180/pi;
%   y = lat*180/pi;
%   z =r;
%   c=X;
%   plot4(x,y,z,c,'.');
%   title('plt4')
%   colorbar
% figure
% plot(Rxlon,Rxlat, '*k' )
% size(Rxlon)
% hold on 
% cd('C:\korea work\TEC MAP')
% lat = load('worldLat.txt');
% long = load('worldLong.txt');
% plot (long,lat,'b')
% grid on 





%end

