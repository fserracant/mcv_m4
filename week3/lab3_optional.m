%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. OPTIONAL: Photo-sequencing with your own images

% 4.1 Take a set of images of a moving scene from different viewpoints at 
%     different time instants. At least two images have to be taken from
%     roughly the same location by the same camera.
%
% 4.2 Implement the first part (until line 16) of the Algorithm 1 of the 
%     Photo-sequencing paper with a selection of the detected dynamic
%     features. You may reuse the code generated for the previous question.
%
clear all;
addpath('../week2');
addpath('../week2/sift'); % ToDo: change 'sift' to the correct path where you have the sift functions

% Read images
im1rgb = imread('Data/sequence3/view0_time0.jpg');
im2rgb = imread('Data/sequence3/view0_time1.jpg');
im3rgb = imread('Data/sequence3/view1_time2.jpg');
im4rgb = imread('Data/sequence3/view2_time3.jpg');
im5rgb = imread('Data/sequence3/view2_time4.jpg');

im1 = sum(double(im1rgb), 3) / 3 / 255;
im2 = sum(double(im2rgb), 3) / 3 / 255;
im3 = sum(double(im3rgb), 3) / 3 / 255;
im4 = sum(double(im4rgb), 3) / 3 / 255;
im5 = sum(double(im5rgb), 3) / 3 / 255;

% show images
figure;
subplot(2,3,1); imshow(im1rgb); axis image; title('Image 1');
subplot(2,3,2); imshow(im2rgb); axis image; title('Image 2');
subplot(2,3,3); imshow(im3rgb); axis image; title('Image 3');
subplot(2,3,4); imshow(im4rgb); axis image; title('Image 4');
subplot(2,3,5); imshow(im4rgb); axis image; title('Image 4');

%% Compute SIFT keypoints
[points_1, desc_1] = sift(im1, 'Threshold', 0.015); % Do not change this threshold!
[points_2, desc_2] = sift(im2, 'Threshold', 0.015);
[points_3, desc_3] = sift(im3, 'Threshold', 0.015);
[points_4, desc_4] = sift(im4, 'Threshold', 0.015);
[points_5, desc_5] = sift(im5, 'Threshold', 0.015);


%% compute the fundamental matrix between 1 and 2
matches12 = siftmatch(desc_1, desc_2);
p1 = [points_1(1:2, matches12(1,:)); ones(1, length(matches12))];
p2 = [points_2(1:2, matches12(2,:)); ones(1, length(matches12))];
[F12, inliers2] = ransac_fundamental_matrix(p1, p2, 2); 

%% compute the fundamental matrix between 1 and 3
matches13 = siftmatch(desc_1, desc_3);
p1 = [points_1(1:2, matches13(1,:)); ones(1, length(matches13))];
p3 = [points_3(1:2, matches13(2,:)); ones(1, length(matches13))];
[F13, inliers3] = ransac_fundamental_matrix(p1, p3, 2); 

%% compute the fundamental matrix between 1 and 4
matches14 = siftmatch(desc_1, desc_4);
p1 = [points_1(1:2, matches14(1,:)); ones(1, length(matches14))];
p4 = [points_4(1:2, matches14(2,:)); ones(1, length(matches14))];
[F14, inliers4] = ransac_fundamental_matrix(p1, p4, 2); 

%% compute the fundamental matrix between 1 and 5
matches15 = siftmatch(desc_1, desc_5);
p1 = [points_1(1:2, matches15(1,:)); ones(1, length(matches15))];
p5 = [points_5(1:2, matches15(2,:)); ones(1, length(matches15))];
[F15, inliers5] = ransac_fundamental_matrix(p1, p5, 2); 

%% Match SIFT keypoints between a and b
figure;
plotmatches(im1, im2, points_1(1:2,:), points_2(1:2,:), matches12, 'Stacking', 'v');

% show inliers
figure;
plotmatches(im1, im2, points_1(1:2,:), points_2(1:2,:), matches12(:,inliers2), 'Stacking', 'v');
title('Inliers');

vgg_gui_F(im1rgb, im2rgb, F12');

%% Plot the object trajectory (keypoint idx_obj_I1 in image 1)
idx_obj_I1 = 307;
idx_obj_I2 = matches12(2, matches12(1,:)==idx_obj_I1);% ToDo: identify the corresponding point of idx_obj_I1 in image 2
idx_obj_I3 = matches13(2, matches13(1,:)==idx_obj_I1); % ToDo: identify the corresponding point of idx_obj_I1 in image 3
idx_obj_I4 = matches14(2, matches14(1,:)==idx_obj_I1); % ToDo: identify the corresponding point of idx_obj_I1 in image 4
idx_obj_I5 = matches15(2, matches15(1,:)==idx_obj_I1); % ToDo: identify the corresponding point of idx_obj_I1 in image 4

%% Show descriptor's anchor points into images
figure();
subplot(2,3,1); imshow(im1); hold on;
show_anchor_points(points_1(:,idx_obj_I1));
title('frame0: view0 @ time0', 'FontSize', 24);
subplot(2,3,2); imshow(im2); hold on;
show_anchor_points(points_2(:,idx_obj_I2));
title('frame1: view0 @ time1', 'FontSize', 24);
subplot(2,3,3); imshow(im3); hold on;
show_anchor_points(points_3(:,idx_obj_I3));
title('frame2: view1 @ time2', 'FontSize', 24);
subplot(2,3,4); imshow(im4); hold on;
show_anchor_points(points_4(:,idx_obj_I4));
title('frame3: view2 @ time3', 'FontSize', 24);
subplot(2,3,5); imshow(im5); hold on;
show_anchor_points(points_5(:,idx_obj_I5));
title('frame4: view2 @ time4', 'FontSize', 24);

%%
% coordinates (in image 1) of the keypoint idx_obj_I1 (point in a obj). 
% point1_1 is the projection of a 3D point in the 3D trajectory of the obj
point1_1 = [points_1(1:2,idx_obj_I1); 1];
% coordinates (in image 1) of another 3D point in the same 3D trajectory of
% the obj
point1_2 = [1070; 210; 1]; % (this is a given data)

% l1 is the projection of the 3D trajectory of keypoint idx_obj_I1
% (it is the line that joins point1_1 and point1_2)
%l1 = F12 * point1_1; % ToDo: compute the line
l1 = cross( point1_1, point1_2); % ToDo: compute the line
% plot the line
figure;imshow(im1);
hold on;
t=1:0.1:2000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(points_1(1,idx_obj_I1), points_1(2,idx_obj_I1), 'y*');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_obj_I1 in image 2
point2 = [points_2(1:2,idx_obj_I2); 1];%
% ToDo: compute the epipolar line of point2 in the reference image
l2 = F12'*point2 ;%
% plot the epipolar line
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'c');
% ToDo: compute the projection of point idx_obj_I2 in the reference image 
pi2 = cross(l1,l2); %
% plot this point
plot(pi2(1)/pi2(3), pi2(2)/pi2(3), 'c*');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_obj_I1 in image 3
point3 = [points_3(1:2,idx_obj_I3); 1];%
% ToDo: compute the epipolar line of point3 in the reference image
l3 = F13'*point3 ;%
% plot the epipolar line
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'b');
% ToDo: compute the projection of point idx_obj_I3 in the reference image
pi3 = cross(l1,l3);%
plot(pi3(1)/pi3(3), pi3(2)/pi3(3), 'b*');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_obj_I1 in image 4
point4 = [points_4(1:2,idx_obj_I4); 1];%
% ToDo: compute the epipolar line of point4 in the reference image
l4 = F14'*point4 ;%
% plot the epipolar line
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'g');
% ToDo: compute the projection of point idx_obj_I4 in the reference image
pi4 = cross(l1,l4);%
plot(pi4(1)/pi4(3), pi4(2)/pi4(3), 'g*');
% 
% ToDo: write the homogeneous coordinates of the corresponding point of idx_obj_I1 in image 4
point5 = [points_5(1:2,idx_obj_I5); 1];%
% ToDo: compute the epipolar line of point5 in the reference image
l5 = F15'*point5 ;%
% plot the epipolar line
plot(t, -(l5(1)*t + l5(3)) / l5(2), 'r');
% ToDo: compute the projection of point idx_obj_I5 in the reference image
pi5 = cross(l1,l5);%
plot(pi5(1)/pi5(3), pi5(2)/pi5(3), 'r*');

legend("trajectory on t0", "obj on t0", ...
       "t1 epipolar line", "object on t1", ...
       "t2 epipolar line", "object on t2", ...
       "t3 epipolar line", "object on t3", ...
       "t4 epipolar line", "object on t4", ...
       "Location", "SouthEast")
