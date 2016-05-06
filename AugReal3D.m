cube=[0,          0
        90,    0
        0,          90
        0,          0
        0,          90
        90,    0
        90,    90
        90,    90];
    cube=horzcat(cube,[0,          1;
0,          1;
0,          1;
90,    1;
90,    1;
90,    1;
0,          1;
90,    1;]);
newCube(:,1) = cube(:,3);
newCube(:,2) = cube(:,2);
newCube(:,3) = cube(:,1);
newCube(:,4) = cube(:,4);
cube=newCube;
cubeEdges = [cube(1,1:3) cube(2,1:3);
             cube(1,1:3) cube(3,1:3);
             cube(1,1:3) cube(4,1:3);
             cube(2,1:3) cube(6,1:3);
             cube(2,1:3) cube(7,1:3);
             cube(3,1:3) cube(5,1:3);
             cube(3,1:3) cube(7,1:3);
             cube(4,1:3) cube(5,1:3);
             cube(4,1:3) cube(6,1:3);
             cube(8,1:3) cube(5,1:3);
             cube(8,1:3) cube(6,1:3);
             cube(8,1:3) cube(7,1:3);];
cubeEdgesX = [cubeEdges(:,1),cubeEdges(:,4)]';
cubeEdgesY = [cubeEdges(:,2),cubeEdges(:,5)]';
cubeEdgesZ = [cubeEdges(:,3),cubeEdges(:,6)]';

figure(), plot3(cube(:,1), cube(:,2), cube(:,3), 'ro')
hold on
plot3(cubeEdgesX, cubeEdgesY, cubeEdgesZ)
hold off
extrinsic = zeros(3,4);
Projection = zeros(3,4);


im = imread('images2.png');
extrinsic(:,:)= [R1(:,:),t1];
Projection(:,:) = A*extrinsic(:,:);
pixels = Projection(:,:)*cube';
pixels = pixels';
for i=1:length(pixels)
pixels(i,:) = pixels(i,:)/pixels(i,3);
end
pixelEdges = [pixels(1,1:3) pixels(2,1:3);
pixels(1,1:3) pixels(3,1:3);
pixels(1,1:3) pixels(4,1:3);
pixels(2,1:3) pixels(6,1:3);                                                
pixels(2,1:3) pixels(7,1:3);
pixels(3,1:3) pixels(5,1:3);
pixels(3,1:3) pixels(7,1:3);
pixels(4,1:3) pixels(5,1:3);
pixels(4,1:3) pixels(6,1:3);
pixels(8,1:3) pixels(5,1:3);
pixels(8,1:3) pixels(6,1:3);
pixels(8,1:3) pixels(7,1:3);];
pixelEdgesX = [pixelEdges(:,1),pixelEdges(:,4)]';
pixelEdgesY = [pixelEdges(:,2),pixelEdges(:,5)]';
pixelEdgesZ = [pixelEdges(:,3),pixelEdges(:,6)]';
disp('3D coordinates of cube for images2.png')
disp(pixelEdges)
figure, imshow(im)
hold on
plot(pixelEdgesX, pixelEdgesY)
plot(pixels(:,1), pixels(:,2), 'ro')
hold off



im = imread('images9.png');
extrinsic(:,:)= [R2(:,:),t2];
Projection(:,:) = A*extrinsic(:,:);
pixels = Projection(:,:)*cube';
pixels = pixels';
for i=1:length(pixels)
pixels(i,:) = pixels(i,:)/pixels(i,3);
end
pixelEdges = [pixels(1,1:3) pixels(2,1:3);
pixels(1,1:3) pixels(3,1:3);
pixels(1,1:3) pixels(4,1:3);
pixels(2,1:3) pixels(6,1:3);                                                
pixels(2,1:3) pixels(7,1:3);
pixels(3,1:3) pixels(5,1:3);
pixels(3,1:3) pixels(7,1:3);
pixels(4,1:3) pixels(5,1:3);
pixels(4,1:3) pixels(6,1:3);
pixels(8,1:3) pixels(5,1:3);
pixels(8,1:3) pixels(6,1:3);
pixels(8,1:3) pixels(7,1:3);];
pixelEdgesX = [pixelEdges(:,1),pixelEdges(:,4)]';
pixelEdgesY = [pixelEdges(:,2),pixelEdges(:,5)]';
pixelEdgesZ = [pixelEdges(:,3),pixelEdges(:,6)]';
disp('3D coordinates of cube for images9.png')
disp(pixelEdges)
figure, imshow(im)
hold on
plot(pixelEdgesX, pixelEdgesY)
plot(pixels(:,1), pixels(:,2), 'ro')
hold off

im = imread('images12.png');
extrinsic(:,:)= [R3(:,:),t3];
Projection(:,:) = A*extrinsic(:,:);
pixels = Projection(:,:)*cube';
pixels = pixels';
for i=1:length(pixels)
pixels(i,:) = pixels(i,:)/pixels(i,3);
end
pixelEdges = [pixels(1,1:3) pixels(2,1:3);
pixels(1,1:3) pixels(3,1:3);
pixels(1,1:3) pixels(4,1:3);
pixels(2,1:3) pixels(6,1:3);                                                
pixels(2,1:3) pixels(7,1:3);
pixels(3,1:3) pixels(5,1:3);
pixels(3,1:3) pixels(7,1:3);
pixels(4,1:3) pixels(5,1:3);
pixels(4,1:3) pixels(6,1:3);
pixels(8,1:3) pixels(5,1:3);
pixels(8,1:3) pixels(6,1:3);
pixels(8,1:3) pixels(7,1:3);];
pixelEdgesX = [pixelEdges(:,1),pixelEdges(:,4)]';
pixelEdgesY = [pixelEdges(:,2),pixelEdges(:,5)]';
pixelEdgesZ = [pixelEdges(:,3),pixelEdges(:,6)]';
disp('3D coordinates of cube for images12.png')
disp(pixelEdges)
figure, imshow(im)
hold on
plot(pixelEdgesX, pixelEdgesY)
plot(pixels(:,1), pixels(:,2), 'ro')
hold off

im = imread('images20.png');
extrinsic(:,:)= [R4(:,:),t4];
Projection(:,:) = A*extrinsic(:,:);
pixels = Projection(:,:)*cube';
pixels = pixels';
for i=1:length(pixels)
pixels(i,:) = pixels(i,:)/pixels(i,3);
end
pixelEdges = [pixels(1,1:3) pixels(2,1:3);
pixels(1,1:3) pixels(3,1:3);
pixels(1,1:3) pixels(4,1:3);
pixels(2,1:3) pixels(6,1:3);                                                
pixels(2,1:3) pixels(7,1:3);
pixels(3,1:3) pixels(5,1:3);
pixels(3,1:3) pixels(7,1:3);
pixels(4,1:3) pixels(5,1:3);
pixels(4,1:3) pixels(6,1:3);
pixels(8,1:3) pixels(5,1:3);
pixels(8,1:3) pixels(6,1:3);
pixels(8,1:3) pixels(7,1:3);];
pixelEdgesX = [pixelEdges(:,1),pixelEdges(:,4)]';
pixelEdgesY = [pixelEdges(:,2),pixelEdges(:,5)]';
pixelEdgesZ = [pixelEdges(:,3),pixelEdges(:,6)]';
disp('3D coordinates of cube for images20.png')
disp(pixelEdges)
figure, imshow(im)
hold on
plot(pixelEdgesX, pixelEdgesY)
plot(pixels(:,1), pixels(:,2), 'ro')
hold off


