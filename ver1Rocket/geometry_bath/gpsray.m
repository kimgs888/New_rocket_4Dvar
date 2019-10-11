function Ray = gpsray( Rx, Tx, Grid, Ext )
% Compute ray structure for GPS propagation
%    Ray = gpsray( Rx, Tx, Grid )
% Input Arguments:
%    Rx   = Cartesian coordinates of receiver, metres, [m,3]
%    Tx   = Cartesian coordinates of transmitter, metres, [m,3]
% Optional Input Arguments:
%    Ext  = Grid extension in radius {[lo,hi]}
%    Grid = Grid structure from setgrid
% Output Arguments:
%    Ray  = Ray structure used by raytrace
% Notes:
%   This function generates a Ray structure as used by the raytracer for
% LINEAR propagation.  Rays start and end points are clipped in radius 
% when they are outside the Grid such that start and end points are 2% 
% from the Grid radial limits.  
%
% Example: Plot per vertex line integral contributions for a ray
%   Grid  = setgrid('Rad',[65:73]*1E5,'Lat',-20:20:20,'Lon',-20:20:20);
%   Ray   = gpsray(sphcart([64E5,[9,3]*pi/180]),sphcart([74E5,[4,3]*pi/180]),Grid);
%   [H,S] = raytracer( Grid, Ray  ); [i,j,v] = find(H);
%   scatter3(Grid.X(j),Grid.Y(j),Grid.Z(j),99,v,'filled'); colourbar; hold on; 
%   geo('Line',Grid); plotrays(S,'Color','k','LineWidth',3); axis vis3d;
%
% See also LINEARTRACE RAYTRACE SETGRID CARTSPH SHELLSECT
% -------------------------------------------------------------------------
if nargin == 4 & ~isempty(Ext)
   if Ext{1}(1)~=0, Grid.Rad = [Grid.Rad(1)+Ext{1}(1),Grid.Rad]; end
   if Ext{1}(2)~=0, Grid.Rad = [Grid.Rad,Grid.Rad(end)+Ext{1}(2)]; end
end

Ray = struct('Prop','LINEAR','A',1,'Rx',[],'Tx',[],'Rv',[],'Tv',[]);
if isempty(Rx) || isempty(Tx), return; end

% Allow for single point specification
if size(Rx,1)==1 && size(Tx,1)>1, Rx = repmat(Rx,size(Tx,1),1); end
if size(Tx,1)==1 && size(Rx,1)>1, Tx = repmat(Tx,size(Rx,1),1); end

% Prevent propagation exactly along grid boundary
Tx = Tx + .99;
Rx = Rx + .99;

rx = cartsph(Rx);
tx = cartsph(Tx);

% Clip rays in radius
i  = find(tx(:,1) > Grid.Rad(end)*1.02); 
Tx(i,:) = shellsect(Rx(i,:),Tx(i,:),Grid.Rad(end)*1.02);
i  = find(tx(:,1) < Grid.Rad(1)*0.98); 
Tx(i,:) = shellsect(Rx(i,:),Tx(i,:),Grid.Rad(1)*0.98);

i  = find(rx(:,1) > Grid.Rad(end)*1.02); 
Rx(i,:) = shellsect(Rx(i,:),Tx(i,:),Grid.Rad(end)*1.02);
i  = find(rx(:,1) < Grid.Rad(1)*0.98); 
Rx(i,:) = shellsect(Rx(i,:),Tx(i,:),Grid.Rad(1)*0.98);

Ray.Prop = 'LINEAR';
Ray.A    = 1;
Ray.Rx   = Rx;
Ray.Rv   = Tx-Rx;
Ray.Tx   = Tx;
Ray.Tv   = Rx-Tx;
