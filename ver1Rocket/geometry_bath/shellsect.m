function [SectA,SectB,Points] = shellsect( A, B, Radius )
% Function to compute intersection point of a ray on a shell
%    [CrdA,CrdB,Points] = shellsect( A, B, r )
% Input Arguments:
%    A      = Cartesian coordinates of point A (m,3)
%    B      = Cartesian coordinates of point B (m,3)
%    Radius = Radius of shell
% Output arguments:
%    SectA  = Cartesian coordinates of first intersection point
%    SectB  = Cartesian coordinates of second intersection point
%    Points = Number of intersections (0,1,2)
% Notes:
%    If only one intersection exists it is placed in SectA
%
% See also
% -----------------------------------------------------------------------
if size(A,1)==1, A = repmat(A,size(B,1),1); end
if size(B,1)==1, B = repmat(B,size(A,1),1); end

a = (B(:,1)-A(:,1)).^2+(B(:,2)-A(:,2)).^2+(B(:,3)-A(:,3)).^2;
b = A(:,1).*(B(:,1)-A(:,1))+A(:,2).*(B(:,2)-A(:,2))+A(:,3).*(B(:,3)-A(:,3));
c = A(:,1).^2+A(:,2).^2+A(:,3).^2-Radius^2;

Root     =  b.^2 - a.*c;
Lambda_a = -b./a + sqrt(Root)./a;
Lambda_b = -b./a - sqrt(Root)./a;

SectA(:,1) = A(:,1) + Lambda_a .* (B(:,1)-A(:,1));
SectA(:,2) = A(:,2) + Lambda_a .* (B(:,2)-A(:,2));
SectA(:,3) = A(:,3) + Lambda_a .* (B(:,3)-A(:,3));

SectB(:,1) = A(:,1) + Lambda_b .* (B(:,1)-A(:,1));
SectB(:,2) = A(:,2) + Lambda_b .* (B(:,2)-A(:,2));
SectB(:,3) = A(:,3) + Lambda_b .* (B(:,3)-A(:,3));

Points = (Lambda_a>0).*(Lambda_a<1) + (Lambda_b>0).*(Lambda_b<1);
Points = Points.*(Root > 0);

Root = repmat( (Root > 0), 1, 3);
La   = repmat( (Lambda_a>0).*(Lambda_a<1), 1, 3);
Lb   = repmat( (Lambda_b>0).*(Lambda_b<1), 1, 3);

SectA = SectA.*Root.*La;
SectB = SectB.*Root.*Lb;

i = find( (sum(SectA.')==0).*(sum(SectB.')~=0) );
SectA(i,:) = SectB(i,:); % Single intersection in SectA


