%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 3: The geometry of two views 
% (application: photo-sequencing)

addpath('../week2');
addpath('../sift'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Compute the fundamental matrix

clear all;
% Two camera matrices for testing purposes
P1 = eye(3,4);
c = cosd(15); s = sind(15);
R = [c -s 0; s c 0; 0 0 1];
t = [.3 0.1 0.2]';
P2 = [R t];
n = 8;
X = [rand(3,n); ones(1,n)] + [zeros(2,n); 3 * ones(1,n); zeros(1,n)];
x1_test = P1 * X;
x2_test = P2 * X;

% Estimated fundamental matrix
% ToDo: create the following function that estimates F using the normalised 8 point algorithm
F_es = fundamental_matrix(x1_test, x2_test);

% Real fundamental matrix
% translation vector must express translation from camera 2 to camera 1
% so t2 need to be as the opposite as t
t2 = t .* -1;
% skew symmetric matrix
S = [0 -t2(3) t2(2); t2(3) 0 -t2(1); -t2(2) t2(1) 0];
% Essential matrix
E = S * R;
% Intrinsic parameters of cameras (let's suppose it is eye(3))
K = eye(3);
F_gt = inv(K)' * E * inv(K); % ToDo: write the expression of the real fundamental matrix for P1 and P2

% Evaluation: these two matrices should be very similar (up-to-scale)
gt_f = F_gt / norm(F_gt);
predicted_F = F_es / norm(F_es) * sign(F_es(1,1));
error_on_F = sum(abs(gt_f(:) - predicted_F(:)));
assert(error_on_F < 1e-9, 'Predicted F is not similar to ground truth')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Robustly fit fundamental matrix

close all;

% Read images
im1rgb = imread('Data/0000_s.png');
im2rgb = imread('Data/0001_s.png');
im1 = sum(double(im1rgb), 3) / 3 / 255;
im2 = sum(double(im2rgb), 3) / 3 / 255;

% show images
figure;
subplot(1,2,1); imshow(im1rgb); axis image; title('Image 1');
subplot(1,2,2); imshow(im2rgb); axis image; title('Image 2');


%% Compute SIFT keypoints

% (make sure that the sift folder provided in lab2 is on the path)

[points_1, desc_1] = sift(im1, 'Threshold', 0.01);
[points_2, desc_2] = sift(im2, 'Threshold', 0.01);

%% Match SIFT keypoints between a and b
matches = siftmatch(desc_1, desc_2);
figure;
plotmatches(im1, im2, points_1(1:2,:), points_2(1:2,:), matches, 'Stacking', 'v');

% p1 and p2 contain the homogeneous coordinates of the matches
p1 = [points_1(1:2, matches(1,:)); ones(1, length(matches))];
p2 = [points_2(1:2, matches(2,:)); ones(1, length(matches))];

% ToDo: create this function (you can use as a basis 'ransac_homography_adaptive_loop.m')
[F, inliers] = ransac_fundamental_matrix(p1, p2, 2, 'sampson'); 

% show inliers
figure;
plotmatches(im1, im2, points_1(1:2,:), points_2(1:2,:), matches(:,inliers), 'Stacking', 'v');
title('Inliers');

vgg_gui_F(im1rgb, im2rgb, F');

%% Plot some epipolar lines

l2 = F * p1; % epipolar lines in image 2 % ToDo
l1 = F' * p2; % epipolar lines in image 1 % ToDo

% choose some random indices
n_inliers = 3;
inlier = randsample(inliers,n_inliers);

% image 1 (plot the three points and their corresponding epipolar lines)
figure;
imshow(im1rgb);

hold on;
for n = 1:length(inlier)
  plot(p1(1, inlier(n)), p1(2, inlier(n)), '+g');
  plot_homog_line(l1(:, inlier(n)));
end
hold off;

% image 2 (plot the three points and their corresponding epipolar lines)
figure;
imshow(im2rgb);

hold on;
for n = 1:length(inlier)
  plot(p2(1, inlier(n)), p2(2, inlier(n)), '+g');
  plot_homog_line(l2(:, inlier(n)));
end
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Photo-sequencing with aerial images

% In this part we will compute a simplified version of the algorithm
% explained in the Photo-sequencing paper. 
% Since we do not have two images
% taken from roughly the same viewpoint at two different time instants we
% will manually pick a dynamic point corresponding to a point in a van 
% (identified by index 'idx_car_I1') and the projection of its 3D trajectory 
% in the reference image. Then we will compute the projection (to the reference image) 
% of three points on this 3D trajectory at three different time instants 
% (corresponding to the time when the three other provided images where taken). 

clear all;

% Read images
im1rgb = imread('Data/frame_00000.tif');
im2rgb = imread('Data/frame_00001.tif');
im3rgb = imread('Data/frame_00002.tif');
im4rgb = imread('Data/frame_00003.tif');

im1 = sum(double(im1rgb), 3) / 3 / 255;
im2 = sum(double(im2rgb), 3) / 3 / 255;
im3 = sum(double(im3rgb), 3) / 3 / 255;
im4 = sum(double(im4rgb), 3) / 3 / 255;

% show images
figure;
subplot(2,2,1); imshow(im1rgb); axis image; title('Image 1');
subplot(2,2,2); imshow(im2rgb); axis image; title('Image 2');
subplot(2,2,3); imshow(im3rgb); axis image; title('Image 3');
subplot(2,2,4); imshow(im4rgb); axis image; title('Image 4');

% Compute SIFT keypoints
[points_1, desc_1] = sift(im1, 'Threshold', 0.015); % Do not change this threshold!
[points_2, desc_2] = sift(im2, 'Threshold', 0.015);
[points_3, desc_3] = sift(im3, 'Threshold', 0.015);
[points_4, desc_4] = sift(im4, 'Threshold', 0.015);

%% ToDo:

% Take image im1 as reference image (image 1) and compute the fundamental 
% matrices needed for computing the trajectory of point idx_car_I1
% (use the SIFT keypoints previously computed)


%% compute the fundamental matrix between 1 and 2
matches12 = siftmatch(desc_1, desc_2);
p1 = [points_1(1:2, matches12(1,:)); ones(1, length(matches12))];
p2 = [points_2(1:2, matches12(2,:)); ones(1, length(matches12))];

[F12, inliers] = ransac_fundamental_matrix(p1, p2, 2); 


%% compute the fundamental matrix between 1 and 3
matches13 = siftmatch(desc_1, desc_3);
p1 = [points_1(1:2, matches13(1,:)); ones(1, length(matches13))];
p3 = [points_3(1:2, matches13(2,:)); ones(1, length(matches13))];

[F13, inliers] = ransac_fundamental_matrix(p1, p3, 2); 

%% compute the fundamental matrix between 1 and 4
matches14 = siftmatch(desc_1, desc_4);
p1 = [points_1(1:2, matches14(1,:)); ones(1, length(matches14))];
p4 = [points_4(1:2, matches14(2,:)); ones(1, length(matches14))];

[F14, inliers] = ransac_fundamental_matrix(p1, p4, 2); 


%% Plot the car trajectory (keypoint idx_car_I1 in image 1)

% ToDo: complete the code

idx_car_I1 = 1197;
idx_car_I2 = matches12(2, matches12(1,:)==idx_car_I1);% ToDo: identify the corresponding point of idx_car_I1 in image 2
idx_car_I3 = matches13(2, matches13(1,:)==idx_car_I1); % ToDo: identify the corresponding point of idx_car_I1 in image 3
idx_car_I4 = matches14(2, matches14(1,:)==idx_car_I1); % ToDo: identify the corresponding point of idx_car_I1 in image 4

% coordinates (in image 1) of the keypoint idx_car_I1 (point in a van). 
% point1_1 is the projection of a 3D point in the 3D trajectory of the van
point1_1 = [points_1(1:2,idx_car_I1); 1];
% coordinates (in image 1) of another 3D point in the same 3D trajectory of
% the van
point1_2 = [334; 697; 1]; % (this is a given data)

% l1 is the projection of the 3D trajectory of keypoint idx_car_I1
% (it is the line that joins point1_1 and point1_2)
%l1 = F12 * point1_1; % ToDo: compute the line
l1 = cross( point1_1, point1_2); % ToDo: compute the line
% plot the line
figure;imshow(im1);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(points_1(1,1197), points_1(2,1197), 'y*');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_car_I1 in image 2
point2 = [points_2(1:2,idx_car_I2); 1];%
% ToDo: compute the epipolar line of point2 in the reference image
l2 = F12'*point2 ;%
% plot the epipolar line
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'c');
% ToDo: compute the projection of point idx_car_I2 in the reference image 
pi2 = cross(l1,l2); %
% plot this point
plot(pi2(1)/pi2(3), pi2(2)/pi2(3), 'c*');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_car_I1 in image 3
point3 = [points_3(1:2,idx_car_I3); 1];%
% ToDo: compute the epipolar line of point3 in the reference image
l3 = F13'*point3 ;%
% plot the epipolar line
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'b');
% ToDo: compute the projection of point idx_car_I3 in the reference image
pi3 = cross(l1,l3);%
plot(pi3(1)/pi3(3), pi3(2)/pi3(3), 'b*');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_car_I1 in image 4
point4 = [points_4(1:2,idx_car_I4); 1];%
% ToDo: compute the epipolar line of point4 in the reference image
l4 = F14'*point4 ;%
% plot the epipolar line
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'g');
% ToDo: compute the projection of point idx_car_I4 in the reference image
pi4 = cross(l1,l4);%
plot(pi4(1)/pi4(3), pi4(2)/pi4(3), 'g*');
legend("trajectory on t0", "van on t0", ...
       "t1 epipolar line", "van on t1", ...
       "t2 epipolar line", "van on t2", ...
       "t3 epipolar line", "van on t3")

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
