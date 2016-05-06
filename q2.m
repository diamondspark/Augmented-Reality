%% q2.1
im= imread('images2.png');
imshow (im)
h1= ginput(4);
im= imread('images9.png');
imshow (im)
h2= ginput(4);
im= imread('images12.png');
imshow (im)
h3= ginput(4);
im= imread('images20.png');
imshow (im)
h4= ginput(4);
h1(:,3)=1;
x1=h1';
h2(:,3)=1;
x2=h2';
h3(:,3)=1;
x3=h3';
h4(:,3)=1;
x4=h4';

 x2_i= [0,0,1;270,0,1;270,210,1;0,210,1]';
%  x2_i=x2_i/1000;
 H1= homography2d(x2_i,x1); 
 H2= homography2d(x2_i,x2); 
 H3=homography2d(x2_i,x3); 
 H4=homography2d(x2_i,x4);
 H1 = H1/H1(3,3);
 H2 = H2/H2(3,3);
 H3 = H3/H3(3,3);
 H4 = H4/H4(3,3);
% Display results
    disp('Homography for images2.png')
    disp(H1)
    disp('Homography for images9.png')
    disp(H2)
    disp('Homography for images12.png')
    disp(H3)
    disp('Homography for images20.png')
    disp(H4)
 %% Q2.2
 H=H1; % selecting homography for which to work
 v12=calcV(H,1,2);
 v11=calcV(H,1,1);
 v22=calcV(H,2,2);
 V1=[v12';(v11-v22)'];
 
  H=H2; % selecting homography for which to work
 v12=calcV(H,1,2);
 v11=calcV(H,1,1);
 v22=calcV(H,2,2);
 V2=[v12';(v11-v22)'];
 
  H=H3; % selecting homography for which to work
 v12=calcV(H,1,2);
 v11=calcV(H,1,1);
 v22=calcV(H,2,2);
 V3=[v12';(v11-v22)'];
 
  H=H4; % selecting homography for which to work
 v12=calcV(H,1,2);
 v11=calcV(H,1,1);
 v22=calcV(H,2,2);
 V4=[v12';(v11-v22)'];
 
 V=[V1;V2;V3;V4];
%   V=V1;
[s,u,v]= svd(V);
b= v(:,end);
B=[b(1),b(2),b(4);b(2),b(3),b(5);b(4),b(5),b(6)];

disp('Intrinsic Params')
v0= (b(2)*b(4)-b(1)*b(5))/(b(1)*b(3)-b(2)*b(2))
lambda = b(6)-(b(4)*b(4)+v0*(b(2)*b(4)-b(1)*b(5)))/b(1)
alpha= sqrt(lambda/b(1))
beta=sqrt(lambda*b(1)/(b(1)*b(3)-b(2)*b(2)))
gamma= -1*b(2)*alpha*alpha*beta/lambda
u0=gamma*v0/alpha-(b(4)*alpha*alpha/lambda)
A=[alpha,gamma,u0;0,beta,v0;0,0,1];
A_inv= inv(A);

%Intrinsic features for all 4 images
disp('Rotation matrix for image2')
[R,t]= intrinsicParams2D(A,H1);
 disp(R)
disp('Translation vector for image2') 
 disp(t)
 
        [U,S,Vprime] = svd(R);
        Rnew = U*Vprime;
       disp('New Rotation matrix for image2 after enforcing constraints')
       disp(Rnew)
       disp('New Rotation matrix * Rotation Matrix_transpose for image2 after enforcing constraints')
       disp(Rnew'*Rnew)
 
 disp('Rotation matrix for image9')
[R,t]= intrinsicParams2D(A,H2);
 disp(R)
disp('Translation vector for image9') 
 disp(t)
 
 [U,S,Vprime] = svd(R);
        Rnew = U*Vprime;
       disp('New Rotation matrix for image9 after enforcing constraints')
       disp(Rnew)
       disp('New Rotation matrix * Rotation Matrix_transpose for image9 after enforcing constraints')
       disp(Rnew'*Rnew)
 
 
 disp('Rotation matrix for image12')
[R,t]= intrinsicParams2D(A,H3);
 disp(R)
disp('Translation vector for image12') 
 disp(t)
 
 [U,S,Vprime] = svd(R);
        Rnew = U*Vprime;
       disp('New Rotation matrix for image12 after enforcing constraints')
       disp(Rnew)
       disp('New Rotation matrix * Rotation Matrix_transpose for image12 after enforcing constraints')
       disp(Rnew'*Rnew)
 
 
 
 disp('Rotation matrix for image20')
[R,t]= intrinsicParams2D(A,H4);
 disp(R)
disp('Translation vector for image20') 
 disp(t)
 
 [U,S,Vprime] = svd(R);
        Rnew = U*Vprime;
       disp('New Rotation matrix for image20 after enforcing constraints')
       disp(Rnew)
       disp('New Rotation matrix * Rotation Matrix_transpose for image20 after enforcing constraints')
       disp(Rnew'*Rnew)
 
 
 %% Q2.3.1
 %All grid corners in 3D world cordinate
 corner_wc=[];
 for i = 1:10
     for j =1:8
        corner_wc=vertcat(corner_wc,[(i-1)*30,(j-1)*30,1]);
     
     end
 end
 
 figure, imshow(imread('images2.png')),title('Projected grid corners for images2')
  hold on
  p_approx1 = H1*corner_wc';
  % Converting into homogenous coordinates by normalizing
  for i=1:size(p_approx1,2)
    p_approx1(:,i)=  p_approx1(:,i)/p_approx1(3,i);
  end
  plot(p_approx1(1,:),p_approx1(2,:),'ro');
  hold off
  
 figure, imshow(imread('images9.png')),title('Projected grid corners for images9')
  hold on
  p_approx2 = H2*corner_wc';
  % Converting into homogenous coordinates by normalizing
  for i=1:size(p_approx2,2)
    p_approx2(:,i)=  p_approx2(:,i)/p_approx2(3,i);
  end
  plot(p_approx2(1,:),p_approx2(2,:),'ro');
  hold off
  
  figure, imshow(imread('images12.png')),title('Projected grid corners for images12')
  hold on
  p_approx3 = H3*corner_wc';
  % Converting into homogenous coordinates by normalizing
  for i=1:size(p_approx3,2)
    p_approx3(:,i)=  p_approx3(:,i)/p_approx3(3,i);
  end
  plot(p_approx3(1,:),p_approx3(2,:),'ro');
  hold off
  
  figure, imshow(imread('images20.png')),title('Projected grid corners for images20')
  hold on
  p_approx4 = H4*corner_wc';
  % Converting into homogenous coordinates by normalizing
  for i=1:size(p_approx4,2)
    p_approx4(:,i)=  p_approx4(:,i)/p_approx4(3,i);
  end
  plot(p_approx4(1,:),p_approx4(2,:),'ro');
  hold off
  %% Q2.3.2
  im2 = imread('images2.png');
  sigma = 2; thresh = 500; radius = 2;
  [cim, r, c, rsubp1, csubp1] = harris(rgb2gray(im2), sigma, thresh, radius, 1);
  title('Figure 2 Harris Corners in image2')
  
   im9 = imread('images9.png');
  sigma = 2; thresh = 500; radius = 2;
  [cim, r, c, rsubp2, csubp2] = harris(rgb2gray(im9), sigma, thresh, radius, 1);
  title('Figure 2 Harris Corners in images9')
  
  im12 = imread('images12.png');
  sigma = 2; thresh = 500; radius = 2;
  [cim, r, c, rsubp3, csubp3] = harris(rgb2gray(im12), sigma, thresh, radius, 1);
  title('Figure 2 Harris Corners in images12')
  
  im20 = imread('images20.png');
  sigma = 2; thresh = 500; radius = 2;
  [cim, r, c, rsubp4, csubp4] = harris(rgb2gray(im20), sigma, thresh, radius, 1);
  title('Figure 2 Harris Corners in images20')
  
  %% Q2.3.3
  D = dist2(p_approx1(1:2,:)',[csubp1, rsubp1]);
  [D_sorted, D_index] = sort(D, 2);
  p_correct(:,:)=[csubp1(D_index(:,1)),rsubp1(D_index(:,1))];
  p_correct(:,3)=1;
  p_correct=p_correct';
  figure,imshow(im2), title('Figure 3 : grid points Images2')
  hold on
  plot(p_correct(1,:),p_correct(2,:),'g+')
  hold off
  p_correct1=p_correct;
  clear p_correct
  
  D = dist2(p_approx2(1:2,:)',[csubp2, rsubp2]);
  [D_sorted, D_index] = sort(D, 2);
  p_correct(:,:)=[csubp2(D_index(:,1)),rsubp2(D_index(:,1))];
  p_correct(:,3)=1;
  p_correct=p_correct';
  figure,imshow(im9), title('Figure 3 : grid points Images9')
  hold on
  plot(p_correct(1,:),p_correct(2,:),'g+')
  hold off
    p_correct2=p_correct;
  clear p_correct
  
    D = dist2(p_approx3(1:2,:)',[csubp3, rsubp3]);
  [D_sorted, D_index] = sort(D, 2);
  p_correct(:,:)=[csubp3(D_index(:,1)),rsubp3(D_index(:,1))];
  p_correct(:,3)=1;
  p_correct=p_correct';
  figure,imshow(im12), title('Figure 3 : grid points Images12')
  hold on
  plot(p_correct(1,:),p_correct(2,:),'g+')
  hold off
    p_correct3=p_correct;
  clear p_correct
  
  D = dist2(p_approx4(1:2,:)',[csubp4, rsubp4]);
  [D_sorted, D_index] = sort(D, 2);
  p_correct(:,:)=[csubp4(D_index(:,1)),rsubp4(D_index(:,1))];
  p_correct(:,3)=1;
  p_correct=p_correct';
 figure, imshow(im20), title('Figure 3 : grid points Images20')
  hold on
  plot(p_correct(1,:),p_correct(2,:),'g+')
  hold off
    p_correct4=p_correct;
  clear p_correct
  
  %% Q2.3.4 new homography
  H1_new= homography2d(corner_wc',p_correct1);
  H2_new= homography2d(corner_wc',p_correct2);
  H3_new= homography2d(corner_wc',p_correct3);
  H4_new= homography2d(corner_wc',p_correct4);
  H1_new= H1_new/H1_new(3,3);
  H2_new= H2_new/H2_new(3,3);
  H3_new= H3_new/H3_new(3,3);
  H4_new= H4_new/H4_new(3,3);
  disp('New Homography for images2')
  disp(H1_new);
  disp('New Homography for images9')
  disp(H2_new);
  disp('New Homography for images12')
  disp(H3_new);
  disp('New Homography for images20')
  disp(H4_new);
  %% Q2.3.5 Calculate intrinsic and extrinsic parameters using Hnew
 v12=calcV(H1_new,1,2);
 v11=calcV(H1_new,1,1);
 v22=calcV(H1_new,2,2);
 V1=[v12';(v11-v22)'];
 
 v12=calcV(H2_new,1,2);
 v11=calcV(H2_new,1,1);
 v22=calcV(H2_new,2,2);
 V2=[v12';(v11-v22)'];
 
 v12=calcV(H3_new,1,2);
 v11=calcV(H3_new,1,1);
 v22=calcV(H3_new,2,2);
 V3=[v12';(v11-v22)'];
 
 v12=calcV(H4_new,1,2);
 v11=calcV(H4_new,1,1);
 v22=calcV(H4_new,2,2);
 V4=[v12';(v11-v22)'];
 
  V=[V1;V2;V3;V4];

[s,u,v]= svd(V);
b= v(:,end);
B=[b(1),b(2),b(4);b(2),b(3),b(5);b(4),b(5),b(6)];

disp('Intrinsic Params')
v0= (b(2)*b(4)-b(1)*b(5))/(b(1)*b(3)-b(2)*b(2))
lambda = b(6)-(b(4)*b(4)+v0*(b(2)*b(4)-b(1)*b(5)))/b(1)
alpha= sqrt(lambda/b(1))
beta=sqrt(lambda*b(1)/(b(1)*b(3)-b(2)*b(2)))
gamma= -1*b(2)*alpha*alpha*beta/lambda
u0=gamma*v0/alpha-(b(4)*alpha*alpha/lambda)
A=[alpha,gamma,u0;0,beta,v0;0,0,1];
A_inv= inv(A);

save('A.mat','A')
%Intrinsic features for all 4 images
disp('Rotation matrix for image2')
[R,t]= intrinsicParams2D(A,H1_new);
 disp(R)
disp('Translation vector for image2') 
 disp(t)
 R1=R;
 t1=t;
 save('k2.mat','R1','t1');
 
        [U,S,Vprime] = svd(R);
        Rnew = U*Vprime;
       disp('New Rotation matrix for image2 after enforcing constraints')
       disp(Rnew)
       disp('New Rotation matrix * Rotation Matrix_transpose for image2 after enforcing constraints')
       disp(Rnew'*Rnew)
 
 disp('Rotation matrix for image9')
[R,t]= intrinsicParams2D(A,H2_new);
 disp(R)
disp('Translation vector for image9') 
 disp(t)
 R2=R;
 t2=t;
 save('k9.mat','R2','t2');
 [U,S,Vprime] = svd(R);
       Rnew = U*Vprime;
       disp('New Rotation matrix for image9 after enforcing constraints')
       disp(Rnew)
       disp('New Rotation matrix * Rotation Matrix_transpose for image9 after enforcing constraints')
       disp(Rnew'*Rnew)
 
 
 disp('Rotation matrix for image12')
[R,t]= intrinsicParams2D(A,H3_new);
 disp(R)
disp('Translation vector for image12') 
 disp(t)
  R3=R;
 t3=t;
 save('k12.mat','R3','t3');
 
 [U,S,Vprime] = svd(R);
        Rnew = U*Vprime;
       disp('New Rotation matrix for image12 after enforcing constraints')
       disp(Rnew)
       disp('New Rotation matrix * Rotation Matrix_transpose for image12 after enforcing constraints')
       disp(Rnew'*Rnew)
 
 
 
 disp('Rotation matrix for image20')
[R,t]= intrinsicParams2D(A,H4_new);
 disp(R)
disp('Translation vector for image20') 
 disp(t)
 R4=R;
 t4=t;
 save('k20.mat','R4','t4');
 
 
 [U,S,Vprime] = svd(R);
        Rnew = U*Vprime;
       disp('New Rotation matrix for image20 after enforcing constraints')
       disp(Rnew)
       disp('New Rotation matrix * Rotation Matrix_transpose for image20 after enforcing constraints')
       disp(Rnew'*Rnew)
 
H1_new= homography2d(corner_wc',p_correct1);
H1_new= H1_new/H1_new(3,3);
p_projection = H1_new*corner_wc';
  for i=1:size(p_projection,2)
    p_projection(:,i)=  p_projection(:,i)/p_projection(3,i);
  end       
              
figure, imshow(im2), title('err_reprojection images2')
hold on
plot(p_correct1(1,:),  p_correct1(2,:),   'g+')
plot(p_projection(1,:), p_projection(2,:),  'bo')
hold off
% 
% err_xline = [p_correct1(1,:) p_projection(1,:)];
% err_yline = [p_correct1(2,:) p_projection(2,:)];


H2_new= homography2d(corner_wc',p_correct2);
H2_new= H2_new/H2_new(3,3);
p_projection = H2_new*corner_wc';
  for i=1:size(p_projection,2)
    p_projection(:,i)=  p_projection(:,i)/p_projection(3,i);
  end       
              
figure,imshow(im9), title('err_reprojection images9')
hold on
plot(p_correct2(1,:),  p_correct2(2,:),   'g+')
plot(p_projection(1,:), p_projection(2,:),  'bo')
hold off
% 
% err_xline = [p_correct1(1,:) p_projection(1,:)];
% err_yline = [p_correct1(2,:) p_projection(2,:)];


H3_new= homography2d(corner_wc',p_correct3);
H3_new= H3_new/H3_new(3,3);
p_projection = H3_new*corner_wc';
  for i=1:size(p_projection,2)
    p_projection(:,i)=  p_projection(:,i)/p_projection(3,i);
  end       
              
figure,imshow(im12), title('err_reprojection images12')
hold on
plot(p_correct3(1,:),  p_correct3(2,:),   'g+')
plot(p_projection(1,:), p_projection(2,:),  'bo')
hold off
% 
% err_xline = [p_correct1(1,:) p_projection(1,:)];
% err_yline = [p_correct1(2,:) p_projection(2,:)];


H4_new= homography2d(corner_wc',p_correct4);
H4_new= H4_new/H4_new(3,3);
p_projection = H4_new*corner_wc';
  for i=1:size(p_projection,2)
    p_projection(:,i)=  p_projection(:,i)/p_projection(3,i);
  end       
              
figure,imshow(im20), title('err_reprojection images20')
hold on
plot(p_correct4(1,:),  p_correct4(2,:),   'g+')
plot(p_projection(1,:), p_projection(2,:),  'bo')
hold off


% err_reprojection1 = mean(sqrt(sum((-p_correct1(:,1:2)'+projection1').^2)))
% err_reprojection2 = mean(sqrt(sum((-p_correct1(:,1:2)'+projection2').^2)))
% err_reprojection3 = mean(sqrt(sum((-p_correct1(:,1:2)'+projection3').^2)))
% err_reprojection4 = mean(sqrt(sum((-p_correct1(:,1:2)'+projection4').^2)))