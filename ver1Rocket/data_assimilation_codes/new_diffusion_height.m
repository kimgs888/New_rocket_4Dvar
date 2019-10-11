function [newdiffusion_height, dist_from_source ] = new_diffusion_height(Heighti,Lat2,Long2,yi,xi,zi,indxh )


point_of_relz = sphcart([64E5+Heighti(indxh)*1e3,[Lat2(indxh),Long2(indxh)].*(pi/180)]);
Pts_from_source = sphcart([64E5+yi(:).*1e3,[xi(:),zi(:)].*(pi/180)]);

dist_from_source = (((Pts_from_source(:,1)-point_of_relz(:,1)).^2 + (Pts_from_source(:,2)-point_of_relz(:,2)).^2 ...
                      + (Pts_from_source(:,3)-point_of_relz(:,3)).^2)); % this is sqaured distance 
                  
                  
                 pointofrelzhgt = 64E5+Heighti(indxh)*1e3; % altitude or rocket (hr)
                  altpointfrmsource = 64E5+yi(:).*1e3;       % altitude of points from source (hn)
                  
                  T = (pointofrelzhgt.^2 + altpointfrmsource.^2); 
                  % look at the derived equation for altitude of  mid point between
                  % rocket and n(r,t), this is the new defusion height ....
                  % We did this to mimic anistrotopic diffusion 
                  newdiffusion_height = (((2*T - dist_from_source)./4).^.5); 