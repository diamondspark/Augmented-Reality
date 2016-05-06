function [ R,t ] = intrinsicParams2D( A,H )
%INTRINSICPARAMS2D Summary of this function goes here
%   Detailed explanation goes here
A_inv=inv(A);
scale=1/norm(A\H(:,1));
r1= scale*A_inv*H(:,1);
r2= scale*A_inv*H(:,2);
 r3 = cross(r1,r2);
 t  = scale*inv(A)*H(:,3);
 R=[r1,r2,r3];

end

