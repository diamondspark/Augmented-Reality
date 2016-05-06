mouseim = imread('4.png');
mouseIn = imresize(mouseim,[210,NaN]); % reshape mouse clipart while preserving its aspect ratio
[rowsM, colsM,~] = size(mouseIn);

  mouseIn = mouseIn + 1;
  mouseIn = imrotate(mouseIn,-90);
% We need to find for every pixel position in mouse image, its
% corresponding co-ordinates in the world scene. This is possible because
% we have the 2D pixel coordinates from mouse image and also the homography
% for each of the 4 grid image on which we want to overlay this mouse
% image.
mouse_pixel_coord=[];
mouse_rgb=[];
for i = 1:colsM
    for j = 1:rowsM
        mouse_pixel_coord= [mouse_pixel_coord;i,j,1];
         mouse_rgb= [mouse_rgb;mouseIn(i,j,1),mouseIn(i,j,2),mouseIn(i,j,3)]; % Correspondng RGB values of each pixel from mouse_pixel_coord
    end
end
mouse_wc= double(mouse_pixel_coord(:,1:3))*H1_new';
%Normalize world coordinates to make it homogeneous co-ordinates i.e. make
%3rd coordinate 1
mouse_wc(:,1)=mouse_wc(:,1)./mouse_wc(:,3);
mouse_wc(:,2)=mouse_wc(:,2)./mouse_wc(:,3);
mouse_wc(:,3)=mouse_wc(:,3)./mouse_wc(:,3);
mouse_wc=abs(round(mouse_wc+1)); % Round off to get integral world coordinates

mouseNew = zeros(max(mouse_wc(:,1)),max(mouse_wc(:,2)),3); % Create a blank RGB image

for i = 1: size(mouse_wc,1)
    r = mouse_wc(i,2);
    c = mouse_wc(i,1);
    mouseNew(r,c,1)= mouse_rgb(i,1);
    mouseNew(r,c,2)= mouse_rgb(i,2);
    mouseNew(r,c,3)= mouse_rgb(i,3);
end

 mouseNew = uint8(mouseNew);
%  figure, imshow(mouseNew),title('aligned clip-art')
    
 bw = (mouseNew > 0) & (mouseNew < 255); % Make all non white (255) pixels in mouse image transparent
% Pad with zeros
dif = abs(size(im2)-size(mouseNew));
mask_pad = padarray(bw, [dif(1),dif(2)], 'post');
img_masked_pad = mouseNew .* uint8(bw);

new_image= zeros(size(im2,1),size(im2,2),3);
for i = 1:size(im2,1)
 for j = 1:size(im2,2)
    if mask_pad(i,j,:)
        new_img(i,j,:) = img_masked_pad(i,j,:);
    else
        new_img(i,j,:) = im2(i,j,:);
    end
 end
end
figure(), imshow(uint8(new_img)),title('Overlayed on Images2')

%% Images9
mouse_wc= double(mouse_pixel_coord(:,1:3))*H2_new';
%Normalize world coordinates to make it homogeneous co-ordinates i.e. make
%3rd coordinate 1
mouse_wc(:,1)=mouse_wc(:,1)./mouse_wc(:,3);
mouse_wc(:,2)=mouse_wc(:,2)./mouse_wc(:,3);
mouse_wc(:,3)=mouse_wc(:,3)./mouse_wc(:,3);
mouse_wc=abs(round(mouse_wc+1)); % Round off to get integral world coordinates

mouseNew = zeros(max(mouse_wc(:,1)),max(mouse_wc(:,2)),3); % Create a blank RGB image

for i = 1: size(mouse_wc,1)
    r = mouse_wc(i,2);
    c = mouse_wc(i,1);
    mouseNew(r,c,1)= mouse_rgb(i,1);
    mouseNew(r,c,2)= mouse_rgb(i,2);
    mouseNew(r,c,3)= mouse_rgb(i,3);
end

 mouseNew = uint8(mouseNew);
%  figure, imshow(mouseNew),title('aligned clip-art')
    
 bw = (mouseNew > 0) & (mouseNew < 255); % Make all non white (255) pixels in mouse image transparent
% Pad with zeros
dif = abs(size(im9)-size(mouseNew));
mask_pad = padarray(bw, [dif(1),dif(2)], 'post');
img_masked_pad = mouseNew .* uint8(bw);

new_image= zeros(size(im9,1),size(im9,2),3);
for i = 1:size(im9,1)
 for j = 1:size(im9,2)
    if mask_pad(i,j,:)
        new_img(i,j,:) = img_masked_pad(i,j,:);
    else
        new_img(i,j,:) = im9(i,j,:);
    end
 end
end
figure(), imshow(uint8(new_img)),title('Overlayed on Images9')

%% Images12
mouse_wc= double(mouse_pixel_coord(:,1:3))*H3_new';
%Normalize world coordinates to make it homogeneous co-ordinates i.e. make
%3rd coordinate 1
mouse_wc(:,1)=mouse_wc(:,1)./mouse_wc(:,3);
mouse_wc(:,2)=mouse_wc(:,2)./mouse_wc(:,3);
mouse_wc(:,3)=mouse_wc(:,3)./mouse_wc(:,3);
mouse_wc=abs(round(mouse_wc+1)); % Round off to get integral world coordinates

mouseNew = zeros(max(mouse_wc(:,1)),max(mouse_wc(:,2)),3); % Create a blank RGB image

for i = 1: size(mouse_wc,1)
    r = mouse_wc(i,2);
    c = mouse_wc(i,1);
    mouseNew(r,c,1)= mouse_rgb(i,1);
    mouseNew(r,c,2)= mouse_rgb(i,2);
    mouseNew(r,c,3)= mouse_rgb(i,3);
end

 mouseNew = uint8(mouseNew);
%  figure, imshow(mouseNew),title('aligned clip-art')
    
 bw = (mouseNew > 0) & (mouseNew < 255); % Make all non white (255) pixels in mouse image transparent
% Pad with zeros
dif = abs(size(im12)-size(mouseNew));
mask_pad = padarray(bw, [dif(1),dif(2)], 'post');
img_masked_pad = mouseNew .* uint8(bw);

new_image= zeros(size(im12,1),size(im12,2),3);
for i = 1:size(im12,1)
 for j = 1:size(im12,2)
    if mask_pad(i,j,:)
        new_img(i,j,:) = img_masked_pad(i,j,:);
    else
        new_img(i,j,:) = im12(i,j,:);
    end
 end
end
figure(), imshow(uint8(new_img)),title('Overlayed on Images12')

%% Images 20
mouse_wc= double(mouse_pixel_coord(:,1:3))*H4_new';
%Normalize world coordinates to make it homogeneous co-ordinates i.e. make
%3rd coordinate 1
mouse_wc(:,1)=mouse_wc(:,1)./mouse_wc(:,3);
mouse_wc(:,2)=mouse_wc(:,2)./mouse_wc(:,3);
mouse_wc(:,3)=mouse_wc(:,3)./mouse_wc(:,3);
mouse_wc=abs(round(mouse_wc+1)); % Round off to get integral world coordinates

mouseNew = zeros(max(mouse_wc(:,1)),max(mouse_wc(:,2)),3); % Create a blank RGB image

for i = 1: size(mouse_wc,1)
    r = mouse_wc(i,2);
    c = mouse_wc(i,1);
    mouseNew(r,c,1)= mouse_rgb(i,1);
    mouseNew(r,c,2)= mouse_rgb(i,2);
    mouseNew(r,c,3)= mouse_rgb(i,3);
end

 mouseNew = uint8(mouseNew);
%  figure, imshow(mouseNew),title('aligned clip-art')
    
 bw = (mouseNew > 0) & (mouseNew < 255); % Make all non white (255) pixels in mouse image transparent
% Pad with zeros
dif = abs(size(im20)-size(mouseNew));
mask_pad = padarray(bw, [dif(1),dif(2)], 'post');
img_masked_pad = mouseNew .* uint8(bw);

new_image= zeros(size(im20,1),size(im20,2),3);
for i = 1:size(im20,1)
 for j = 1:size(im20,2)
    if mask_pad(i,j,:)
        new_img(i,j,:) = img_masked_pad(i,j,:);
    else
        new_img(i,j,:) = im20(i,j,:);
    end
 end
end
figure(), imshow(uint8(new_img)),title('Overlayed on Images20')