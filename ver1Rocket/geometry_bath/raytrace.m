% MEX function for ray tracing through a 3D cartesian vertex grid
%   [H,[S]] = raytrace( Grid, Ray );
% Input Arguments:
%   Grid = Grid definition structure with the fields;
%          X[m,n,o](double)   = Cartesian vertex coordinates in x
%          Y[m,n,o](double)   = Cartesian vertex coordinates in y
%          Z[m,n,o](double)   = Cartesian vertex coordinates in z
%          F[m,n,o](double)   = Per vertex field values, optional
%   Ray  = Ray definition structure with the fields;
%          Prop[1,:](char)    = Ray propagator for refractive index RI;
%             'LINEAR'          : RI   = 1 (Grid.F not required)
%             'CIRCULAR'        : RI   = A / F
%             'PARABOLIC'       : RI^2 = 1 - A F, H contains phase path
%             'PARABOLIC_GROUP' : RI^2 = 1 - A F, H contains group path
%          A[1],[r,1](double) = Propagation constant, global or per ray(r)
%          Tx[r,3](double)    = Transmitter cartesian coordinates
%          Tv[r,3](double)    = Transmitter cartesian propagation vector
%          Ix[r,3](double)    = Transmitter vertex indices, 0..m-1,...
%          Rx[r,3](double)    = Optional receiver cartesian coordinates
%          Rv[r,3](double)    = Optional receiver cartesian reception vector
% Ouput Arguments:
%   H[r,v](sparse double)     = Optical path for rays(r), vertices(v)
%   S[4,:](double)            = Ray grid intersections, [r;x;y;z,...]
%
% Notes:
%   This function computes the ray paths for rays propagating through a
% linear gradient medium defined by the field variable, F, and the
% relationship this variable has to the refractive index, Prop.
%   Rays are traced from points Tx with starting vectors Tv. With no Rx,Rv
% defined the ray terminates on exit from the grid.  With Rx,Rv defined
% termination occurs when the ray intersects Rx to within |Rv| in the
% plane perpendicular to Rv. Tx must lie within the grid and the indices
% of the voxel containing Tx must be provided in Ix.
%   The matrix H contains a per vertex decomposition of the optical path
% where each row corresponds to a ray, r, and each column to a sequentially
% ordered vertex, v = sub2ind(size(Grid.X),m,n,o). The optical path matrix
% is defined such that sum(H,2) is the optical path integral of 'RI ds' for
% path element ds along the ray (ds/RI for group path). Failure for a ray
% results in sum(H(r,:),2)=0 for ray index r.
%   The S matrix contains the ray indices and x,y,z coordinates of the
% intersections of the ray with the grid.  The first point is the ray
% starting point.
%   Grids are wrapped in the third dimension if the seperation of the first
% and last vertices is less than twice that of adjacent vertices.
%   Shortcuts for the fields Tx,Tv,Rx,Rv,A can be applied such that a single
% value or vector applies to all rays.
%
%
% See also RAYTRACER LINEARTRACE SETGRID GPSRAY GRIDINDEX PLOTRAYS
% ----------------------------------------------------------------------------
