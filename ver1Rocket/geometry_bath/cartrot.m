function Y = cartrot( Matrix, X  )
% Function to rotate cartesian coordinates given rotation matrix
%    Y = cartrot( Matrix, X  )
% Input arguments:
%    Matrix = Rotation matrix [3,3]
%    X      = Cartesian coordinates, [n,3] or cell array {x,y,z}
% Output arguments:
%    Y      = Rotated coordinates
%
% Notes:
%
% See also SPHCART CARTSPH
% --------------------------------------------------------------------------
if ~iscell(X)
   Y = [sum(repmat(Matrix(1,:),size(X,1),1).*X,2),...
        sum(repmat(Matrix(2,:),size(X,1),1).*X,2),...
        sum(repmat(Matrix(3,:),size(X,1),1).*X,2)];
else
   Y = {reshape(sum(repmat(Matrix(1,:),length(X{1}(:)),1).*[X{1}(:),X{2}(:),X{3}(:)],2),size(X{1})),...
        reshape(sum(repmat(Matrix(2,:),length(X{1}(:)),1).*[X{1}(:),X{2}(:),X{3}(:)],2),size(X{1})),...
        reshape(sum(repmat(Matrix(3,:),length(X{1}(:)),1).*[X{1}(:),X{2}(:),X{3}(:)],2),size(X{1}))};
end
