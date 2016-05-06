% Generate 16x12 P Matrix
clear all
M = [2,2,2,1;-2,2,2,1;-2,2,-2,1;2,2,-2,1;2,-2,2,1;-2,-2,2,1;-2,-2,-2,1;2,-2,-2,1];
m= [422,323,1;178,323,1;118,483,1;482,483,1;438,73,1;162,73,1;78,117,1;522,117,1];
%% Q1.1
figure(1), plot(m(:,1),m(:,2),'o'), title('Image Points')
%% Q1.3
P=[];
for i = 1:8
    P=vertcat(P,twoRowP(M(i,:),m(i,:)))
end
clear M
%% Q1.4
[U,S,V]=svd(P);
M(1,1:4) = V(1:4,end);
M(2,1:4) = V(5:8,end);
M(3,1:4) = V(9:12,end);
disp('Projection Matrix M')
disp(M)
%% Q1.5
[U,S,V]=svd(M);
temp=V(:,4)/V(4,4);
Euclidiean_Camera_Coord=[temp(1),temp(2),temp(3)]
%% Q1.6
 M_prime= M(1:3,1:3);
 M_prime= M_prime/M_prime(3,3);
 disp('Matrix M''')
 disp(M_prime)
 %% Q1.7
 cos_theta_x= M_prime(3,3)/sqrt(M_prime(3,3)^2+M_prime(3,2)^2);
 sin_theta_x= -M_prime(3,2)/sqrt(M_prime(3,3)^2 + M_prime(3,2)^2);
 theta_x_deg=rad2deg(acos(cos_theta_x))
 R_x=[1,0,0;0,cos_theta_x,-sin_theta_x;0,sin_theta_x,cos_theta_x]
 N=M_prime*R_x
%  N=N/N(3,3);
 %% Q1.8
 cos_theta_z= (N(2,2)/sqrt(N(2,1)^2+N(2,2)^2));
 sin_theta_z= -N(2,1)/sqrt(N(2,1)^2+N(2,2)^2);
 theta_z_deg= rad2deg(acos(cos_theta_z))
 R_z=[cos_theta_z,-sin_theta_z,0;sin_theta_z,cos_theta_z,0;0,0,1]
%  N= N*R_z;
 %% Q1.9
 % R= R1*R2*R3
 R=R_x*eye(3)*R_z
 K= M_prime*(R);
 K = K/K(3,3)
fl= (K(1,1)+K(2,2))/2;
disp('Focal length = ')
disp(fl)
disp('Image center = ')
disp([K(1,3),K(2,3)])