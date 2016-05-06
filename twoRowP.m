%% 1.2 
function [P]= twoRowP(M,m)
    % M= Corner Point in 3d
     %m= Image of corner in 2D
     P= [M,0,0,0,0,m(1).*M;0,0,0,0,M,m(2).*M];
end