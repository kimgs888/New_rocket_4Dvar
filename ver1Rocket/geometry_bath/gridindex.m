function Ix = gridindex( Grid, Crd )
% Compute indices into grid of cartesian point(s)
%    Ix = gridindex( Grid, Crd )
% Input Arguments:
%    Grid  = Grid definition structure, see setgrid
%    Crd   = Cartesian coordinates of point(s), [m,3]
% Output Arguments:
%    Ix   = Fractional index into grid, [m,3]
%
% Notes:
%    Computes the vertex indices of points located within an optionally
% rotated grid volume. Spherical grids must be uniform in latitude and 
% longitude, but may be non-uniform in the radial. Cartesian grids must be
% monotonic in each dimension.  All indices are clamped to be within 
% the grid.
%
% Example: Get vertex indices of a point within the grid volume
%    Grid = setgrid('Rad',64E5:1E5:74E5,'Lat',-80:10:80,'Lon',-80:10:80);
%    Ix   = gridindex(Grid,sphcart([64E5,1*pi/180,1*pi/180]))
% Example: Plot interpolated radial density distribution
%    Grid = setgrid('Rad',64E5:5E4:74E5,'Lat',-80:10:80,'Lon',-80:10:80);
%    Grid = setgrid(Grid,'F',iri(Grid,datenum(2003,1,1,12,0,0)));
%    Rad  = [64E5:1E4:74E5].';
%    Ix   = gridindex(Grid,[Rad,repmat([2,2]*pi/180,length(Rad),1)]);
%    Ne   = interpn(Grid.F,Ix(:,1),Ix(:,2),Ix(:,3),'spline');
%    plot(Ne,Rad/1000-6400,'.-'); 
%
% See also RAYTRACE SETGRID CARTSPH
% -------------------------------------------------------------------------
if isfield(Grid,'Rot'), R = Grid.Rot.'; else R = [1,0,0;0,1,0;0,0,1]; end
if isempty(Crd), Ix = []; return; end

Crd = [sum(repmat(R(1,:),size(Crd,1),1).*Crd,2),...
       sum(repmat(R(2,:),size(Crd,1),1).*Crd,2),...
       sum(repmat(R(3,:),size(Crd,1),1).*Crd,2)];

if isfield(Grid,'Lat')
   p       = cartsph(Crd);
   dLat    = (Grid.Lat(2)-Grid.Lat(1))*pi/180;
   dLon    = (Grid.Lon(2)-Grid.Lon(1))*pi/180;
   Lat     = round((p(:,2)-Grid.Lat(1)*pi/180)/dLat)*dLat + Grid.Lat(1)*pi/180;
   Lon     = round((p(:,3)-Grid.Lon(1)*pi/180)/dLon)*dLon + Grid.Lon(1)*pi/180;
   r       = sphcart([ones(size(Lat,1),1),Lat,Lon]);
   cosa    = dot(Crd,r,2)./(p(:,1)+1E-9);
   p(:,1)  = p(:,1) + p(:,1).*(1-cosa);

   Ix(:,1) = interp1([0,Grid.Rad,1E9],[0,1:length(Grid.Rad)+1],p(:,1),'linear');
   Ix(:,2) = (p(:,2)-Grid.Lat(1)*pi/180)./dLat + 1;
   Ix(:,3) = (p(:,3)-Grid.Lon(1)*pi/180)./dLon + 1;
   Size    = [length(Grid.Rad),length(Grid.Lat),length(Grid.Lon)];
elseif isfield(Grid,'Xc')
   Ix(:,1) = interp1([-1E9,Grid.Xc(2:end-1),1E9],1:length(Grid.Xc),Crd(:,1));
   Ix(:,2) = interp1([-1E9,Grid.Yc(2:end-1),1E9],1:length(Grid.Yc),Crd(:,2));
   Ix(:,3) = interp1([-1E9,Grid.Zc(2:end-1),1E9],1:length(Grid.Zc),Crd(:,3));
   Size    = [length(Grid.Xc),length(Grid.Yc),length(Grid.Zc)];
else
   error('Grids must be spherical or cartesian');
end

Ix(Ix<1) = 1;
i = find(Ix(:,1)>Size(1)); Ix(i,1)=Size(1);
i = find(Ix(:,2)>Size(2)); Ix(i,1)=Size(2);
i = find(Ix(:,3)>Size(3)); Ix(i,1)=Size(3);
