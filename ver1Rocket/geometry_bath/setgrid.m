function Grid = setgrid( varargin )
% Create grid structure for representing cartesian vertex grids
%    Grid = setgrid( varargin )
% Required Input Params (structures or param,value pairs):
%   'Rad','Lat','Lon' - Spherical grid axes (angles in degrees)
%   'Xc','Yc','Zc'    - OR cartesian grid axes
% Optional Input Params:
%   'F'      - Scalar field or cell array vector field
%   'Rot'    - Cartesian rotation matrix [3,3]
%   'Interp' - Interpolations in each dimension [1,3]
%   'Method' - Inperpolation method for F, default 'spline'
%   'Sample' - Sampling cell array of strings, default {':',':',':'}
%   'Wrap'   - Wrapping in each dimension, boolean [1,3] 
%    boolean fxn in matlab coverts the vector X into boolean. by nick
% Output Arguments:
%    Grid = Grid structure with appended cartesian vertex coords;
%    'X'     - Cartesian vertex x coordinates
%    'Y'     - Cartesian vertex y coordinates
%    'Z'     - Cartesian vertex z coordinates
%
% Notes:
%   This function computes the three-dimensional matrices of vertex
% coordinates defined by spherical or cartesian axis vectors using
% MatLAB's function ndgrid.  Wrapping grids should not have duplicate
% points, ie set longitudes from -180:10:170 only wrap when plotting.
%
% The function rayentry requires that grids be defined with uniform 
% monatonic changes in the Lat and Lon fields, Rad may be variable.
%
% Example: Construct a grid, set interpolation and wrapping and plot
%    Grid = setgrid('Rad',6400E3+[100,1000]*1E3,'Lat',[20:5:50],...
%                   'Lon',[-180:10:170]);
%    Grid = setgrid(Grid,'Interp',[0,1,1],'Wrap',[0,0,1]);
%    geo('OutLine'); geo('Line',Grid,'Color',[0,.7,0]);
%
% See also NDGRID STRUCTCAT STRUCTSECT CARTSPH SPHCART CARTROT
% ----------------------------------------------------------------------
Ref  = struct('Xc',[],'Yc', [],'Zc',  [],'Rad',[],'Lat',[],'Lon',[],'F',[],...
              'Rot',[],'Wrap',[],'Interp',[],'Sample',[]);

Grid = structcat('Method','spline',varargin);

Grid = structcat(Grid,structsect(struct(Ref),Grid));


if isfield(Grid,'Rad') && isfield(Grid,'Lat') && isfield(Grid,'Lon') % Spherical grid
   x = Grid.Rad(:).';
   y = Grid.Lat(:).';
   z = Grid.Lon(:).';
elseif isfield(Grid,'Xc') && isfield(Grid,'Yc') && isfield(Grid,'Zc') % Cartesian grid
   x = Grid.Xc(:).';
   y = Grid.Yc(:).';
   z = Grid.Zc(:).';
else
   warning('Grid requires fields Rad,Lat,Lon or Xc,Yc,Zc'); return
end

if isfield(Grid,'F') && length(Grid.F(:))==length(x)*length(y)*length(z)
   Si = [length(x),length(y),length(z)];
   if iscell(Grid.F)
      Grid.F = {reshape(Grid.F{1},Si),reshape(Grid.F{2},Si),reshape(Grid.F{3},Si)};
   else Grid.F = reshape(Grid.F,Si); end
end


if isfield(Grid,'Sample') && ~isempty(Grid.Sample)
   if ~isempty(Grid.Sample{1}), Sx = Grid.Sample{1}; else Sx = ':'; end
   if ~isempty(Grid.Sample{2}), Sy = Grid.Sample{2}; else Sy = ':'; end
   if ~isempty(Grid.Sample{3}), Sz = Grid.Sample{3}; else Sz = ':'; end

   x = eval(['x(',Sx,')']); y = eval(['y(',Sy,')']); z = eval(['z(',Sz,')']);
   x = reshape(x,1,length(x));
   y = reshape(y,1,length(y));
   z = reshape(z,1,length(z));

   if isfield(Grid,'F')
      if ~iscell(Grid.F), Grid.F = eval(['Grid.F(',Sx,',',Sy,',',Sz,')']);
      else Grid.F = {eval(['Grid.F{1}(',Sx,',',Sy,',',Sz,')']),...
                     eval(['Grid.F{2}(',Sx,',',Sy,',',Sz,')']),...
                     eval(['Grid.F{3}(',Sx,',',Sy,',',Sz,')'])}; end
   end
   Grid = rmfield(Grid,'Sample');
end

if isfield(Grid,'Wrap') && ~isempty(Grid.Wrap)
   if Grid.Wrap(1) & ~isempty(x)
      x = [x , x(end)+x(2)-x(1)];
      if isfield(Grid,'F')
         if ~iscell(Grid.F), Grid.F = cat(1,Grid.F,Grid.F(1,:,:));
         else Grid.F = {cat(1,Grid.F{1},Grid.F{1}(1,:,:)),...
                 cat(1,Grid.F{2},Grid.F{2}(1,:,:)),cat(1,Grid.F{3},Grid.F{3}(1,:,:))};
         end
      end
   end
   if Grid.Wrap(2) & ~isempty(y)
      y = [y , y(end)+y(2)-y(1)];
      if isfield(Grid,'F')
         if ~iscell(Grid.F), Grid.F = cat(2,Grid.F,Grid.F(:,1,:));
         else Grid.F = {cat(2,Grid.F{1},Grid.F{1}(:,1,:)),...
                 cat(2,Grid.F{2},Grid.F{2}(:,1,:)),cat(2,Grid.F{3},Grid.F{3}(:,1,:))};
         end
      end
   end
   if Grid.Wrap(3) & ~isempty(z)
      z = [z , z(end)+z(2)-z(1)];
      if isfield(Grid,'F')
         if ~iscell(Grid.F), Grid.F = cat(3,Grid.F,Grid.F(:,:,1));
         else Grid.F = {cat(3,Grid.F{1},Grid.F{1}(:,:,1)),...
                 cat(3,Grid.F{2},Grid.F{2}(:,:,1)),cat(3,Grid.F{3},Grid.F{3}(:,:,1))};
         end
      end
   end
   Grid = rmfield(Grid,'Wrap');
end


if isfield(Grid,'Interp') && ~isempty(Grid.Interp)
   if length(Grid.Interp) == 1, Grid.Interp = repmat(Grid.Interp,1,3); end
   xi = [1:1/(Grid.Interp(1)+1):length(x)];
   if length(x)>1, x = interp1(x,xi,'linear'); end
   yi = [1:1/(Grid.Interp(2)+1):length(y)];
   if length(y)>1, y = interp1(y,yi,'linear'); end
   zi = [1:1/(Grid.Interp(3)+1):length(z)];
   if length(z)>1, z = interp1(z,zi,'linear'); end

   if isfield(Grid,'F') && ~isempty(Grid.F)
      if ~iscell(Grid.F) 
         if     length(xi)==1, Grid.F = interpn(squeeze(Grid.F),yi,zi.',Grid.Method);
         elseif length(yi)==1, Grid.F = interpn(squeeze(Grid.F),xi,zi.',Grid.Method);
         elseif length(zi)==1, Grid.F = interpn(squeeze(Grid.F),xi,yi.',Grid.Method);
         else Grid.F = interpn(Grid.F,xi.',yi,zi,Grid.Method); end
      else Grid.F = {interpn(Grid.F{1},xi.',yi,zi,Grid.Method),...
                     interpn(Grid.F{2},xi.',yi,zi,Grid.Method),...
                     interpn(Grid.F{3},xi.',yi,zi,Grid.Method)}; end
      if ~iscell(Grid.F), Grid.F = reshape(Grid.F,[length(xi),length(yi),length(zi)]);
      else Grid.F = {reshape(Grid.F{1},[length(xi),length(yi),length(zi)]),...
                     reshape(Grid.F{2},[length(xi),length(yi),length(zi)]),...
                     reshape(Grid.F{3},[length(xi),length(yi),length(zi)])}; end
   end
   Grid = rmfield(Grid,'Interp');
end

[Grid.X,Grid.Y,Grid.Z] = ndgrid(x,y,z);

if isfield(Grid,'Rad') && isfield(Grid,'Lat') && isfield(Grid,'Lon') % Spherical grid
   Crd = sphcart([Grid.X(:),Grid.Y(:)*pi/180,Grid.Z(:)*pi/180]);
   Grid.Rad = x;
   Grid.Lat = y;
   Grid.Lon = z;
else
   Crd = [Grid.X(:),Grid.Y(:),Grid.Z(:)];
   Grid.Xc = x;
   Grid.Yc = y;
   Grid.Zc = z;
end

if isfield(Grid,'Rot') && ~isempty(Grid.Rot), Crd = cartrot(Grid.Rot,Crd); end

Grid.X = reshape(Crd(:,1),size(Grid.X));
Grid.Y = reshape(Crd(:,2),size(Grid.Y));
Grid.Z = reshape(Crd(:,3),size(Grid.Z));

Grid = rmfield(Grid,'Method');
