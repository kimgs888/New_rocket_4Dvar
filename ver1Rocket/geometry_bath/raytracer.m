function [H,S] = raytracer( Grid, Ray, varargin )
% Wrapper for raytrace with extended grid volume
%     [H,[S]] = raytracer( Grid, Ray, varargin )
% Input Arguments:
%    Grid = Grid structure
%    Ray  = Ray structure
% Optional Input Arguments:
%    'Ext',Ext     = Grid extensions, {[lo,hi],[lo,hi],[lo,hi]}
%    'Iter',Iter   = Number of iterations for non-linear trace
%    'Angle',Angle = Initial angle, degree, for non-linear trace
%    'Relax,Relax  = Angular relaxation factor for non-linear trace  
% Output Arguments:
%    H    = Decomposed path within each voxel, metres [m,n]
% Optional Output Arguments:
%    S    = Ray-path/grid intersection coordinates from raytrace
% Notes:
%    This function provides a wrapper for the function raytrace that
% allows for the fact that transmission points may be outside of a grid.
% The input grid is extended along each dimension by Ext{dim}(lo,hi), 
% where lo and hi represent offsets to be added to the first and last
% terms on the grid axis for each dimension. For spherical grids only
% radial extensions are supported.  By default spherical grid extensions
% are defined as 5% the radial limits.  
%   The parameter Iter enables an iterative approach to estimating the
% transmission vector that minimises the departure of the ray from the
% receiver.  
%
% Example: Acoustic trace through volume from external rays
%   Grid   = setgrid('Xc',-5:5,'Yc',-5:5,'Zc',-5:5);
%   Grid   = setgrid(Grid,'F',1+6*exp( -(Grid.X.^2+Grid.Y.^2+Grid.Z.^2)/10 ) );
%   Ray    = struct('Prop','CIRCULAR','A',1); N = 100;
%   Ray.Tx = sphcart([ones(N,1)*10,rand(N,1)*pi-pi/2,rand(N,1)*2*pi-pi]);
%   Ray.Tv = -Ray.Tx;
%   [H,S]  = raytracer(Grid,Ray,'Ext',{[-1,1],[-1,1],[-1,1]}.*10);
%   h      = plotrays(S,'CVal',sum(H,2)); axis vis3d; grid on; view(45,45);
%
% See also RAYTRACE LINEARTRACE SETGRID GRIDINDEX GPSRAY
% -------------------------------------------------------------------------
if isempty(Ray.Tx), H = zeros(0,length(Grid.X(:))); S = []; return; end

Param.Ext   = {};
Param.Iter  = 0;
Param.Angle = 1;
Param.Relax = 0.5;
Param = structcat(Param,varargin);

if isempty(Param.Ext)
   if isfield(Grid,'Xc'),      Param.Ext = {[-1,1],[-1,1],[-1,1]}; 
   elseif isfield(Grid,'Rad'), Param.Ext = {[-Grid.Rad(1)*.05,Grid.Rad(end)*.05]}; 
   else error('Not a standard Grid structure'); end
end

% Construct extended grid
if isfield(Grid,'Xc') % Cartesian
   Del = {min(abs(diff(Grid.Xc))),min(abs(diff(Grid.Yc))),min(abs(diff(Grid.Zc)))}./1E7;
   Xc  = Grid.Xc; Vx  = ones(1,length(Grid.Xc));
   Yc  = Grid.Yc; Vy  = ones(1,length(Grid.Yc));
   Zc  = Grid.Zc; Vz  = ones(1,length(Grid.Zc));

   if Param.Ext{1}(1)~=0, Xc = [Xc(1)+[Param.Ext{1}(1),-Del{1}],Xc];  Vx=[0,0,Vx]; end
   if Param.Ext{1}(2)~=0, Xc = [Xc,Xc(end)+[Del{1},Param.Ext{1}(2)]]; Vx=[Vx,0,0]; end

   if Param.Ext{2}(1)~=0, Yc = [Yc(1)+[Param.Ext{2}(1),-Del{2}],Yc];  Vy=[0,0,Vy]; end
   if Param.Ext{2}(2)~=0, Yc = [Yc,Yc(end)+[Del{2},Param.Ext{2}(2)]]; Vy=[Vy,0,0]; end

   if Param.Ext{3}(1)~=0, Zc = [Zc(1)+[Param.Ext{3}(1),-Del{3}],Zc];  Vz=[0,0,Vz]; end
   if Param.Ext{3}(2)~=0, Zc = [Zc,Zc(end)+[Del{3},Param.Ext{3}(2)]]; Vz=[Vz,0,0]; end

   G = setgrid('Xc',Xc,'Yc',Yc,'Zc',Zc);
   pause 
elseif isfield(Grid,'Rad') % Spherical, radial extension only
   Del = min(abs(diff(Grid.Rad)))./1E7;   
   Rad = Grid.Rad; Vx  = ones(1,length(Grid.Rad));
   Lat = Grid.Lat; Vy  = ones(1,length(Grid.Lat));
   Lon = Grid.Lon; Vz  = ones(1,length(Grid.Lon));

   if Param.Ext{1}(1)~=0, Rad = [Rad(1)+[Param.Ext{1}(1),-Del],Rad];  Vx=[0,0,Vx]; end
   if Param.Ext{1}(2)~=0, Rad = [Rad,Rad(end)+[Del,Param.Ext{1}(2)]]; Vx=[Vx,0,0]; end

   G = setgrid('Rad',Rad,'Lat',Lat,'Lon',Lon);
   if isfield(Grid,'Rot'), G = setgrid(G,'Rot',Grid.Rot); end
else error('Unsupported geometry'); end

% Indices of vertices in origional grid
[Vx,Vy,Vz] = ndgrid(Vx,Vy,Vz); V = find(Vx(:) & Vy(:) & Vz(:)); Vo = find(Vx(:)==0 | Vy(:)==0 | Vz(:)==0);


if isfield(Grid,'F') % Refractive index outside grid is set to unity
   if strcmp(Ray.Prop,'CIRCULAR'), G.F = ones(size(G.X))*Ray.A;
   else G.F = zeros(size(G.X)); end
   G.F(V) = Grid.F(1:length(V(:))); 
end

% Allow for shortcuts
if size(Ray.Tx,1)==1 && size(Ray.Tv,1)>1, Ray.Tx = repmat(Ray.Tx,size(Ray.Tv,1),1); end

% Compute vertex indices of transmission points within the grid
Ray.Ix = fix(gridindex(G,Ray.Tx))-1;

if Param.Iter >= 1
   R     = Ray;
   R.Rx  = repmat(Ray.Rx,3,1);
   R.Rv  = repmat(Ray.Rv,3,1);
   R.Tx  = repmat(Ray.Tx,3,1);
   R.Ix  = repmat(Ray.Ix,3,1);
end

for i = 1:Param.Iter
% Orthogonal vectors to transmission vector
   Ray.Tv = normal(Ray.Tv,2,2);
   a = cross(Ray.Tv,circshift(Ray.Tv,[0,1]),2);
   b = cross(a,Ray.Tv,2);
p
% Scale to angular range
   a = a * Param.Angle /45;
   b = b * Param.Angle /45;

% Construct triad of rays about central ray at 120 degrees sepn
   R.Tv  = repmat(Ray.Tv,3,1) + [a ; (-a-b)/sqrt(2) ; (-a+b)/sqrt(2)];

   if isfield(G,'Rad') && ~isempty(Param.Ext{1}) % Tx outside spherical grid
      [R.Tx,R.Ix] = sect(G,repmat(Ray.Tx,3,1),R.Ix,R.Tv,G.Rad(end-1)+Del);
   end

 
   [H,S] = raytrace( G, R );

% Update, Require uniform distance scaling
   x = rayexit(S); x(sum(H,2)==0,1:3) = NaN;
   n = log(1/3)/log(2/3);
   s = size(Ray.Tx,1);
   d = reshape(sqrt( sum((R.Rx-x).^2, 2 ) ),size(Ray.Tx,1),3);
   D = (d(:,1)+d(:,2)+d(:,3)).^n;
   Ray.Tv = R.Tv(1:s,:)      .*repmat( (d(:,2)+d(:,3)).^n./D,1,3) + ...
            R.Tv(s+1:2*s,:)  .*repmat( (d(:,1)+d(:,3)).^n./D,1,3) + ...
            R.Tv(2*s+1:end,:).*repmat( (d(:,1)+d(:,2)).^n./D,1,3);
   Ray.Tv = real(Ray.Tv);
   Param.Angle = Param.Angle*Param.Relax;
end

R = Ray;
if isfield(G,'Rad') && ~isempty(Param.Ext{1}) % Tx outside spherical grid
   [R.Tx,R.Ix] = sect(G,Ray.Tx,Ray.Ix,Ray.Tv,G.Rad(end-1)+Del);
else Dx = []; end
 
if nargout == 1
   H = raytrace(G,R); 
else 
 
   [H,S] = raytrace(G,R);
   
   if isfield(G,'Rad') && ~isempty(Param.Ext{1}) && isfield(Ray,'Rx') % Tx outside spherical grid
      L = sum(H(:,Vo),2);
      i = find(diff([-1,S(1,:)])); 
      j = S(1,i)+1; 
%       Ray.Tx(j,:)
%       repmat(L(j),[1,3])
%      normal(Ray.Tx(j,:)-Ray.Rx(j,:),2,2)
%      repmat(L(j),[1,3])

      S(2:4,i) = (Ray.Tx(j,:)+ normal(Ray.Tx(j,:)-Ray.Rx(j,:),2,2).*repmat(L(j),[1,3])).';
      
   end
end

H = H(:,V);


% Correct for Tx outside spherical grid
function [Tx,Ix] = sect(G,Tx,Ix,Tv,Shell)
   Rad  = sqrt(sum(Tx.^2,2));
   i    = find( (Rad > G.Rad(end)) & ~isnan(Tv(:,1)) );
   if ~isempty(i)
      To(i,:) = Tx(i,:);   
      Tx(i,:) = shellsect(Tx(i,:),Tx(i,:)+ normal(Tv(i,:),2,2).*repmat(Rad(i),1,3),Shell);
      Ix(i,:) = fix(gridindex(G,Tx(i,:)))-1;
   end


