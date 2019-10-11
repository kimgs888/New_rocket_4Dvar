function [NeedH,H,v]=creat_H_matrix (S, Time)

NeedH = [];
n_file = dir(S);
%figure 
%hold on 
for file = 1:length(n_file)
  
 CMNfile = n_file(file).name;  % now get the file name 
 
 %Call this file and extract the Rx, and latitude and logitude to
%generate the Ray in constracting the the H matrix 
[Rx,needPRN,needLong,needLat,needSTEC,flag]= read_CMN_extract_Rx_pos (CMNfile,Time);
% if flag ==1  %check to see if the time was available, if not,  do not do any thing
%   for i = 1:1%length(needLong)
%     needLat = 21.6271; %75 25.4237;  % 80 29.1919;
%     needLong =123.3198; %75  123.9566; %80 124.6300; %3198;
%     hstep = 50E3;
%     latstep = 10; 
%     longstep = 10; 
%     Lg = [20; 30];
%     Lh = [120; 130];
%     Lr = [64E5; 64E5];
%     ndLat = [22; 25];
%     ndLong = [125; 130];
%     ndht = [64E5+20000E3 64E5+20000E3]; 
%  Lv =   sphcart([Lr,[Lg(:),Lh(:)]*pi/180]);
%  Lf =    sphcart([ndht(:),[ndLat(:),ndLong(:)]*pi/180]);
%  
%   Grid  = setgrid('Rad',[64E5+90E3:hstep:64E5+1500E3],'Lat',20:latstep :40,'Lon',120:longstep:130);
%   Ray   = gpsray(sphcart([64E5,[Rx(2),Rx(3)]*pi/180]),sphcart([(64E5+20000E3),[needLat(i),needLong(i)]*pi/180]),Grid);
%   Ray1 = gpsray(Lv,Lf,Grid)  ;
%  [H,S]  = raytracer(Grid,Ray1);
%   
%  % [H,S] = raytracer(Grid, Ray);
%   [v] = find(H);
%   V = (full(H(v)))';
%   sizev = size(V);
% Length_MIDAS = sum(V);
%  dlmwrite('H.txt', V/1E3)
%  [ri,j, m] = find(H);
%    gx = Grid.X(:); % : was j
%    gy = Grid.Y(:); % : was j
%    gz =Grid.Z(:); % : was j
%   girdcord =  [gx gy gz];
%   
%  
%  Rud =[];
%  ltu=[];
%  lnu =[]; 
%  xo = [];
%  
%  for ln = 1:length(Grid.X(:))
%    
%  Ru = cartsph([gx(ln),gy(ln),gz(ln)]);
%  cm = (Ru(1)-64E5)/1E3;
%  if cm<0
%      cm = 0;
%  end 
%  Rud =[Rud; cm];
%  ltu=[ltu; Ru(2)*180/pi];
%  lnu =[lnu; Ru(3)*180/pi];
%  longpt = Ru(3);
%  latpt = Ru(2);
% 
%   heightbegin = cm;
%   heightend = cm;
%   heightstep = 1;
%   foF2value = 10; 
%  dn = 10;
%  hour = 10;
%  year = 2010;
%   cd('C:\korea work\TEC MAP')
%    
%   [selectedfof2,outcolne, heightrange] = profile_for_H_matrix(longpt,latpt, year, hour, dn,heightbegin,heightend, heightstep, foF2value);
%   xo = [xo; outcolne(1)];
%   cd('C:\korea work\temography\nk\nikiz')
%  end 
%  K = 10;
%  b = [30 32]';
%  A = H; 
%  x0 = xo;
%  [X,info] = kaczmarz(A,b,K,x0);
%  plot(X, 'g')
%  figure
%  plot(xo, 'r')
% 
% %  length(xo)
% %  length(H)
% 
%  
% % dlmwrite('coordinates.txt', coord, 'delimiter', '\t','precision',6) 
%  
% %   NBB- check slice box for plotting 
% 
% %    plotrays(S,'Color','k','LineWidth',2);% axis vis3d;
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
end




