function [NeedH,H,v]=creat_H_matrix2 (S,Time,year,doy)

 NeedH = [];
 H = 0;
 v = 0;
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


[info,flag]= read_CMN_extract_Rx_pos (CMNfile,Time);
% if flag ==1  %check to see if the time was available, if not,  do not do any thing
     
stlat = info.satlat;
stlon = info.satlon;
stalt = info.Satalt;
Rxlat = info.RxLat;
Rxlon = info.RxLon;
Rxalt = info.RxAlt;
rxlt = repmat(64E5, size(stlat));



Rx = sphcart([64E5+Rxalt',[Rxlat',Rxlon']*pi/180]);
Tx = sphcart([(64E5+stalt)',[stlat',stlon']*pi/180]);
hstep = 750E3;
latstep = 10; 
longstep = 10; 
% Creat the  H matrix
   Grid  = setgrid('Rad',[64E5+90E3:hstep:64E5+1500E3],'Lat',20:latstep:50,'Lon',120:longstep:150);
   Ray   = gpsray(Rx,Tx,Grid);
   [H] = raytracer(Grid, Ray);
 
  
  

% extract the cartesian cordinates and use them to obtain xo intial guess.
  gx = Grid.X(:); % : was j
  gy = Grid.Y(:); % : was j
  gz =Grid.Z(:); % : was j
  

 for ln = 1:length(Grid.X(:))   
 Ru = cartsph([gx(ln),gy(ln),gz(ln)]);
 height = (Ru(1)-64E5)/1E3;
 if height<0
     height = 0;
 end
%setup input to iri for that specific grid point
latpt = Ru(2)*180/pi;
longpt = Ru(3)*180/pi;
heightbegin =height;
heightend = height;
ele.lat(ln) = latpt;
ele.lon(ln) = longpt;
ele.hght(ln) = heightend;
heightstep = 1;
foF2value = 10; %use this as dummy input for this iri model but its not used since jchoice is set to 0(stardard)
dn = doy;
hour = Time;

 % cd('C:\korea work\TEC MAP')
  % call this function to generate the electron density at that grid point 
   [selectedfof2,outcolne, heightrange] = profile_for_H_matrix(longpt,latpt, year, hour, dn,heightbegin,heightend, heightstep, foF2value);
  xintial.xo(ln) = outcolne(1); % store the values as intial guess from IRI code 
%  cd('C:\korea work\temography\nk\nikiz')
 end 
%for K = 1:3:15
K = 3 ;  % Now set up the number iterations 
b = info.stec';
A = H; 
size(A)
x0 = xintial.xo';
%   [X,infm] = kaczmarz(A,b,K,x0);
%    figure
%    plot(X, 'g')
%    title('kr')
%    hold on 
%    plot(x0, 'b')
%     title('kr')
 
[X,infm] = MART1(A,b,K,x0);
%figure

%  plot(abs(sqrt(x0*1e6/1.24e10)), 'b')
%  hold on
%  title('orgn')
%  figure
%  plot(abs(sqrt(X*1e6/1.24e10)), 'r')
%  title('mart')
%  hold on 
 
% [X,infm] = MART1(A,b,K,x0);
% figure
% plot(abs(sqrt(x0*1e6/1.24e10)), 'b')
%  hold on 
%  plot(abs(sqrt(X*1e6/1.24e10)), 'r')
%  title('mart')
%  hold on 
 
%   figure
%  scatter3( ele.lat,ele.lon,ele.hght,99,X,'filled');
%   colorbar
%   title('scter X')
%  figure
%  scatter3(ele.lat,ele.lon,ele.hght,99,x0,'filled');
%  title('scter xo')
%  colorbar
%  figure 
  x = ele.lat;
  y = ele.lon;
  z =ele.hght;
  c=X;
  plot4(x,y,z,c,'.');
   title('plt4')
  
figure 
[x,y,z,f] = ndgrid(Grid.X(:), Grid.Y(:),Grid.Z(:),c);
scatter3(x(:),y(:),z(:),40,f(:),'filled')
colorbar
  
% end 

% 
% %  length(xo)
% %  length(H)
% 
%  
% dlmwrite('coordinates.txt', coord, 'delimiter', '\t','precision',6) 
%  
% %   NBB- check slice box for plotting 
% 
 %  plotrays(S,'Color','k','LineWidth',2); axis vis3d;
% % %   size(S);
% % %   hold on 
% % %   grid on 
% %    plot3(Rx(2),Rx(3),64E5,'*r')
% %      set(gca,'xtick',[20:30:50])
% %      set(gca,'ytick',[120:30:150])
% %      set(gca,'ztick',[64E5:1500E3:64E5+1500E3])
% % %  % tx = ['station[',num2str(Rx(2)), ' ', num2str(Rx(3)), ']'];
% % %  % text(Rx(2),Rx(3),64E5,tx) 
% %    xlim([33 37])
% %    ylim([125 127])
% %    zlim([64E5 80E5])
% %    view(99, 12)
% 
% 
%  
%   
%  % plot (needLat(i),needLong(i),'*g','Linewidth',2)
%  
% %   for k =0:10E3:1500E3
% %        m=0;
% %   for  i = 1:7
% %  
% %  line([20+m 20+m],[120 150],[65E5 65E5+k],  'color', 'k','LineStyle',':' ); % Horizontal lines in increment of 10 steps 
% %  line([20 50],[120+m 120+m],[65E5 65E5+k],  'color', 'k','LineStyle',':' ); %  vertical lines  in increment of 10 steps
% %  m = m+5;
% %   end
% %   
% %  end 
% % zlim ([65E5 65E5+1500E3])
% %  xlim([35 37])
% %   ylim([125 126]) 
%   
% %gx = Grid.X;
% 
% %  whos
% % %  % Ru = [];
% % %   gx = Grid.X(j);
% % %   gy = Grid.Y(j);
% % %   gz =Grid.Z(j);
% % %  for ln = 1:length(Grid.X(j))
% % %  Ru = cartsph([gx(ln),gy(ln),gz(ln)]);
% % %  Rud = Ru(1);
% % %  ltu= Ru(2)*180/pi;
% % %  lnu = Ru(3)*180/pi;
% % %  hold on
% % %  plot3(ltu,lnu,Rud,'.','LineWidth',3)
% % %  end 
% %  v =  size(v);
% %  s =  size(S);
% %  xt = [];
% %  yt = [];
% %  zt = [];
% %  for mn = 1:length(S)
% % c = S(:,mn) ;
% % Rx = cartsph([c(2),c(3),c(4)]);
% % Rd = Rx(1);
% % lt= Rx(2)*180/pi;
% % ln = Rx(3)*180/pi;
% % %plot3(c(2),c(3),c(4)); axis vis3d; view(45,45);
% % plot3(lt,ln,Rd,'.','LineWidth',3)
% % xt = [xt;lt];
% % yt = [yt;ln];
% % zt = [zt,Rd];
% % 
% % 
% % 
% %  end 
% % [v] = find(H);
% % ph = H(v);
% % sum(ph);
% % 
% % size(v);
% % size(zt);
% %set(gca,'YTick',yt) 
% %    [ii,j,v] = find(H);
% %   plotrays(S,'Color','k','LineWidth',3);% axis vis3d;
% %   hold on 
% % % scatter3(Grid.X(j),Grid.Y(j),Grid.Z(j),99,v,'filled');
% % for k =1:length(Grid.X(j))
% %     
% %    plot3(Grid.X(j),Grid.Y(j),Grid.Z(j), '*')
% % end 
%  
% %   xlabel('x')
% %   xlabel('y')
% %   zlabel('z')
% %     
% %   Y = [H needSTEC(i)];
% %   NeedH = [NeedH; Y];
% %   x = Grid.X;
% %   y = Grid.Y;
% %   z = Grid.Z;
% %   h = plotrays(S);
% %   grid on 
% 
%   end 
% end 
% % end 
% % size(  NeedH )
% %    [i,j,v] = find(H);
% % 
% % 
% % Grid  = setgrid('Rad',[50:20E3:20E3],'Lat',20:1:50,'Lon',120:1:150)
% % Assume that the reciever on surface of the earth with constant radius RE
% % Ray   = gpsray(sphcart([64E5,[Rx(2),Rx(3)]*pi/180]),sphcart([2.6611E7,[latn(i),lonn(i)]*pi/180]),Grid);
% % [H,S] = raytracer( Grid, Ray);
% 
% % for  i = 1:7
% % line([120 150],[20+m 20+m] ,  'color', 'k','LineStyle',':' ); % Horizontal lines in increment of 10 steps 
% % line([120+m 120+m],[20 50] ,  'color', 'k','LineStyle',':' ); %  vertical lines  in increment of 10 steps
% % m = m+5;
% %end
% %  plot(Rx(2),Rx(3),'*r')
%  for i = 1:1%length(needLong)
% %   line([Rx(2),needLat(i)],[Rx(3),needLong(i)], [0 20000])
% % %  line([Rx(2),needLat(i)],[Rx(3),needLong(i)])
% %  
% %   xlim([20 50])
% %   ylim ([120 150])
% %   xlabel('Latitude (deg)')
% %   ylabel('Longitude (deg)')
% %   zlabel('Altitude (km)')
% %   zlim([0 20000])
% %  view(-5,7); % set the view, elevation, azimuth 
% % %  grid on 
% %  % hold on; 
%  end
%%end
% data = [];
% for i = 1:length(ele.lat)
%     
% end 



