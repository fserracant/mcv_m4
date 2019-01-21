%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 4: Reconstruction from two views (knowing internal camera parameters) 

close all;
clear all;
addpath('../sift'); % ToDo: change 'sift' to the correct path where you have the sift functions
addpath('../week3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Triangulation

% ToDo: create the function triangulate.m that performs a triangulation
%       with the homogeneous algebraic method (DLT)
%
%       The entries are (x1, x2, P1, P2, imsize), where:
%           - x1, and x2 are the Euclidean coordinates of two matching 
%             points in two different images.
%           - P1 and P2 are the two camera matrices
%           - imsize is a two-dimensional vector with the image size

% Test the triangulate function
% Use this code to validate that the function triangulate works properly

P1 = eye(3,4);
c = cosd(15); s = sind(15);
R = [c -s 0; s c 0; 0 0 1];
t = [.3 0.1 0.2]';
P2 = [R t];
n = 8;
X_test = [rand(3,n); ones(1,n)] + [zeros(2,n); 3 * ones(1,n); zeros(1,n)];
x1_test = euclid(P1 * X_test);
x2_test = euclid(P2 * X_test);

N_test = size(x1_test,2);
X_trian = zeros(4,N_test);
for i = 1:N_test
    X_trian(:,i) = triangulate(x1_test(:,i), x2_test(:,i), P1, P2, [2 2]);
end

% Check error (euclidean distance between vectors)
error = sqrt(sum((euclid(X_test) - euclid(X_trian)).^2));
assert(sum(error) < 1e-12, 'BUG on triangulate function')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Reconstruction from two views

%% Read images
Irgb{1} = imread('Data/0001_s.png');
Irgb{2} = imread('Data/0002_s.png');
I{1} = sum(double(Irgb{1}), 3) / 3 / 255;
I{2} = sum(double(Irgb{2}), 3) / 3 / 255;
[h,w] = size(I{1});


%% Compute keypoints and matches.
points = cell(2,1);
descr = cell(2,1);
for i = 1:2
    [points{i}, descr{i}] = sift(I{i}, 'Threshold', 0.01);
    points{i} = points{i}(1:2,:);
end

matches = siftmatch(descr{1}, descr{2});

% Plot matches.
figure();
plotmatches(I{1}, I{2}, points{1}, points{2}, matches, 'Stacking', 'v');


%% Fit Fundamental matrix and remove outliers.
x1 = points{1}(:, matches(1, :));
x2 = points{2}(:, matches(2, :));
[F, inliers] = ransac_fundamental_matrix(homog(x1), homog(x2), 2.0);

% Plot inliers.
inlier_matches = matches(:, inliers);
figure;
plotmatches(I{1}, I{2}, points{1}, points{2}, inlier_matches, 'Stacking', 'v');

x1 = points{1}(:, inlier_matches(1, :));
x2 = points{2}(:, inlier_matches(2, :));

vgg_gui_F(Irgb{1}, Irgb{2}, F');




%% Compute candidate camera matrices.

% Camera calibration matrix
K = [2362.12 0 1520.69; 0 2366.12 1006.81; 0 0 1];
scale = 0.3;
H = [scale 0 0; 0 scale 0; 0 0 1];
K = H * K;

% ToDo: Compute the Essential matrix from the Fundamental matrix
% from F = inv(K)' * E * inv(K) we derive expression for E
E = K' * F * K;

% ToDo: write the camera projection matrix for the first camera
P1 = K * [eye(3) zeros(3,1)];

% ToDo: write the four possible matrices for the second camera
[U, D, V] = svd(E);
W = [ 0 -1 0; 1 0 0; 0 0 1 ];
R1 = W  * V';
R2 = W' * V';

% HINT: You may get improper rotations; in that case you need to change
%       their sign.
% Let R be a rotation matrix, you may check:
% if det(R) < 0
%     R = -R;
% end
if det(R1) < 0
    R1 = -R1;
end
if det(R2) < 0
    R2 = -R2;
end

Pc2 = {};
Pc2{1} = K * [U * R1        U(:, end)];
Pc2{2} = K * [U * R1   -1 * U(:, end)];
Pc2{3} = K * [U * R2        U(:, end)];
Pc2{4} = K * [U * R2   -1 * U(:, end)];

% plot the first camera (red) and the four possible solutions for the 
%second camera (blue)
figure;
subplot(2,2,1);
plot_camera(P1,w,h,1,'r');
plot_camera(Pc2{1},w,h);
title('Solution1')

subplot(2,2,2);
plot_camera(P1,w,h,1,'r');
plot_camera(Pc2{2},w,h);
title('Solution2')

subplot(2,2,3);
plot_camera(P1,w,h,1,'r');
plot_camera(Pc2{3},w,h);
title('Solution3')

subplot(2,2,4);
plot_camera(P1,w,h,1,'r');
plot_camera(Pc2{4},w,h);
title('Solution4')

%% Reconstruct structure
% ToDo: Choose a second camera candidate by triangulating a match.

% we know camera 2 is on the right of camera 1 and both point the same
% direction (parallel configuration)
O1 = optical_center(P1);
D1 = view_direction(P1, x1(:,1));
C2idx = 0;
for i = 1:4
    O2 = optical_center(Pc2{i});
    D2 = view_direction(Pc2{i}, x1(:,1));
    
    if O2(1) > 0 % camera 2 is on the right of camera 1
        if isequal(sign(D1), sign(D2)) % both on the same direction
            C2idx = i;
        end
    end 
end

assert(C2idx > 0, "No candidate found");
P2 = Pc2{C2idx};

% Highlight picked solution
subplot(2,2,C2idx); box on; ax = gca; ax.LineWidth = 2;

% Triangulate all matches.
N = size(x1,2);
X = zeros(4,N);
for i = 1:N
    X(:,i) = triangulate(x1(:,i), x2(:,i), P1, P2, [w h]);
end



%% Plot with colors
r = interp2(double(Irgb{1}(:,:,1)), x1(1,:), x1(2,:));
g = interp2(double(Irgb{1}(:,:,2)), x1(1,:), x1(2,:));
b = interp2(double(Irgb{1}(:,:,3)), x1(1,:), x1(2,:));
Xe = euclid(X);
figure; hold on;
plot_camera(P1,w,h);
plot_camera(P2,w,h);
for i = 1:length(Xe)
    scatter3(Xe(1,i), Xe(3,i), -Xe(2,i), 5^2, [r(i) g(i) b(i)]/255, 'filled');
end
axis equal;
hold off;


%% Compute reprojection error.

% ToDo: compute the reprojection errors
%       plot the histogram of reprojection errors, and
%       plot the mean reprojection error

x1_hat = euclid(P1 * X);
x2_hat = euclid(P2 * X);
RE1 = sqrt(sum((x1 - x1_hat) .^ 2));
RE2 = sqrt(sum((x2 - x2_hat) .^ 2));

figure;
subplot(1,2,1);
h1 = histogram(RE1);
m1 = mean(RE1);
ylim=get(gca,'ylim');
ax = line([m1 m1], ylim, 'Color', 'r');
title({"Histogram of", "reprojection errors for Image 1"});
xlabel("pixels");
legend(ax, 'mean');

subplot(1,2,2);
h2 = histogram(RE2);
m2 = mean(RE2);
ylim=get(gca,'ylim');
ax = line([m2 m2], ylim, 'Color', 'r');
title({"Histogram of", "reprojection errors for Image 2"});
xlabel("pixels");
legend(ax, 'mean');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Depth map computation with local methods (SSD)

% Data images: 'scene1.row3.col3.ppm','scene1.row3.col4.ppm'
% Disparity ground truth: 'truedisp.row3.col3.pgm'

% Write a function called 'stereo_computation' that computes the disparity
% between a pair of rectified images using a local method based on a matching cost 
% between two local windows.
% 
% The input parameters are 5:
% - left image
% - right image
% - minimum disparity
% - maximum disparity
% - window size (e.g. a value of 3 indicates a 3x3 window)
% - matching cost (the user may able to choose between SSD and NCC costs)
%
% In this part we ask to implement only the SSD cost
%
% Evaluate the results changing the window size (e.g. 3x3, 9x9, 20x20,
% 30x30) and the matching cost. Comment the results.
%
% Note 1: Use grayscale images
% Note 2: For this first set of images use 0 as minimum disparity 
% and 16 as the the maximum one.


left  = imread('Data/scene1.row3.col3.ppm');
right = imread('Data/scene1.row3.col4.ppm');
max_d = 16;
min_d = 0;

figure;
subplot(2,2,1);
D = stereo_computation(left, right, min_d, max_d, 3, 'ssd');
imshow(D);
title("3x3 (SSD)");

subplot(2,2,2);
D = stereo_computation(left, right, min_d, max_d, 9, 'ssd');
imshow(D);
title("9x9 (SSD)");

subplot(2,2,3);
D = stereo_computation(left, right, min_d, max_d, 20, 'ssd');
imshow(D);
title("20x20 (SSD)");

subplot(2,2,4);
D = stereo_computation(left, right, min_d, max_d, 30, 'ssd');
imshow(D);
title("30x30 (SSD)");



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Depth map computation with local methods (NCC)

% Complete the previous function by adding the implementation of the NCC
% cost.
%
% Evaluate the results changing the window size (e.g. 3x3, 9x9, 20x20,
% 30x30) and the matching cost. Comment the results.

figure();
subplot(2,2,1);
D = stereo_computation(left, right, min_d, max_d, 3, 'ncc');
imshow(D);
title("3x3 (NCC)");

subplot(2,2,2);
D = stereo_computation(left, right, min_d, max_d, 9, 'ncc');
imshow(D);
title("9x9 (NCC)");

subplot(2,2,3);
D = stereo_computation(left, right, min_d, max_d, 20, 'ncc');
imshow(D);
title("20x20 (NCC)");

subplot(2,2,4);
D = stereo_computation(left, right, min_d, max_d, 30, 'ncc');
imshow(D);
title("30x30 (NCC)");



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONAL:  Depth map computation with local methods (SAD)

% Complete the previous function by adding the implementation of the SAD
% cost.
%
% Evaluate the results changing the window size (e.g. 3x3, 9x9, 20x20,
% 30x30) and the matching cost. Comment the results.

figure();
subplot(2,2,1);
D = stereo_computation(left, right, min_d, max_d, 3, 'sad');
imshow(D);
title("3x3 (SAD)");

subplot(2,2,2);
D = stereo_computation(left, right, min_d, max_d, 9, 'sad');
imshow(D);
title("9x9 (SAD)");

subplot(2,2,3);
D = stereo_computation(left, right, min_d, max_d, 20, 'sad');
imshow(D);
title("20x20 (SAD)");

subplot(2,2,4);
D = stereo_computation(left, right, min_d, max_d, 30, 'sad');
imshow(D);
title("30x30 (SAD)");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Depth map computation with local methods

% Data images: '0001_rectified_s.png','0002_rectified_s.png'

% Test the functions implemented in the previous section with the facade
% images. Try different matching costs and window sizes and comment the
% results.
% Notice that in this new data the minimum and maximum disparities may
% change.


left  = imread('Data/0001_rectified_s.png');
right = imread('Data/0002_rectified_s.png');

max_d = 30;
min_d = 0;
figure;
subplot(2,2,1);
D = stereo_computation(left, right, min_d, max_d, 3, 'ssd');
imshow(D);
title("3x3 (SSD)");

subplot(2,2,2);
D = stereo_computation(left, right, min_d, max_d, 9, 'ssd');
imshow(D);
title("9x9 (SSD)");

subplot(2,2,3);
D = stereo_computation(left, right, min_d, max_d, 20, 'ssd');
imshow(D);
title("20x20 (SSD)");

subplot(2,2,4);
D = stereo_computation(left, right, min_d, max_d, 30, 'ssd');
imshow(D);
title("30x30 (SSD)");

figure();
subplot(2,2,1);
D = stereo_computation(left, right, min_d, max_d, 3, 'ncc');
imshow(D);
title("3x3 (NCC)");

subplot(2,2,2);
D = stereo_computation(left, right, min_d, max_d, 9, 'ncc');
imshow(D);
title("9x9 (NCC)");

subplot(2,2,3);
D = stereo_computation(left, right, min_d, max_d, 20, 'ncc');
imshow(D);
title("20x20 (NCC)");

subplot(2,2,4);
D = stereo_computation(left, right, min_d, max_d, 30, 'ncc');
imshow(D);
title("30x30 (NCC)");

figure();
subplot(2,2,1);
D = stereo_computation(left, right, min_d, max_d, 3, 'sad');
imshow(D);
title("3x3 (SAD)");

subplot(2,2,2);
D = stereo_computation(left, right, min_d, max_d, 9, 'sad');
imshow(D);
title("9x9 (SAD)");

subplot(2,2,3);
D = stereo_computation(left, right, min_d, max_d, 20, 'sad');
imshow(D);
title("20x20 (SAD)");

subplot(2,2,4);
D = stereo_computation(left, right, min_d, max_d, 30, 'sad');
imshow(D);
title("30x30 (SAD)");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Bilateral weights

% Modify the 'stereo_computation' so that you can use bilateral weights (or
% adaptive support weights) in the matching cost of two windows.
% Reference paper: Yoon and Kweon, "Adaptive Support-Weight Approach for Correspondence Search", IEEE PAMI 2006
%
% Comment the results and compare them to the previous results (no weights).
%
% Note: Use grayscale images (the paper uses color images)


left  = imread('Data/scene1.row3.col3.ppm');
right = imread('Data/scene1.row3.col4.ppm');
max_d = 16;
min_d = 0;

figure();
subplot(2,2,1);
D = stereo_computation_BW(left, right, min_d, max_d, 3, 'BW');
imshow(D);
title("3x3 (BW)");

subplot(2,2,2);
D = stereo_computation_BW(left, right, min_d, max_d, 9, 'BW');
imshow(D);
title("9x9 (BW)");

subplot(2,2,3);
D = stereo_computation_BW(left, right, min_d, max_d, 20, 'BW');
imshow(D);
title("20x20 (BW)");

subplot(2,2,4);
D = stereo_computation_BW(left, right, min_d, max_d, 30, 'BW');
imshow(D);
title("30x30 (BW)");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONAL:  Stereo computation with Belief Propagation

% Use the UGM library used in module 2 and implement a  
% stereo computation method that minimizes a simple stereo energy with 
% belief propagation. 
% For example, use an L2 or L1 pixel-based data term (SSD or SAD) and 
% the same regularization term you used in module 2. 
% Or pick a stereo paper (based on belief propagation) from the literature 
% and implement it. Pick a simple method or just simplify the method they propose.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONAL:  Depth computation with Plane Sweeping

% Implement the plane sweeping method explained in class.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONAL:  Depth map fusion 

% In this task you are asked to implement the depth map fusion method
% presented in the following paper:
% B. Curless and M. Levoy. A Volumetric Method for Building Complex
% Models from Range Images. In Proc. SIGGRAPH, 1996.
%
% 1. Use the set of facade images 00xx_s.png to compute depth maps 
% corresponding to different views (and optionally from different pairs of 
% images for the same view).
% 2. Then convert each depth map to a signed distance function defined in 
% a disretized volume (using voxels).
% 3. Average the different signed distance functions, the resulting 
% signed distance is called D.
% 4. Set as occupied voxels (those representing the surface) those 
% where D is very close to zero. The rest of voxels will be considered as 
% empty.
%
% For that you need to compute a depth map from a pair of views in general
% position (non rectified). Thus, you may either use the plane sweep
% algorithm (if you did it) or the local method for estimating depth
% (mandatory task) together with the following rectification method which 
% has an online demo available: 
% http://demo.ipol.im/demo/m_quasi_euclidean_epipolar_rectification/


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONAL:  New view synthesis

% In this task you are asked to implement part of the new view synthesis method
% presented in the following paper:
% S. Seitz, and C. Dyer, View morphing, Proc. ACM SIGGRAPH 1996.

% You will use a pair of rectified stereo images (no need for prewarping
% and postwarping stages) and their corresponding ground truth disparities
% (folder "new_view").
% Remember to take into account occlusions as explained in the lab session.
% Once done you can apply the code to the another pair of rectified images 
% provided in the material and use the estimated disparities with previous methods.
clc; clear variables; close all;

addpath('../sift'); 
addpath('../week3'); 
addpath('../week2'); 

%% Load images and disparity files
Irgb_0 = imread('Data/new_view/im0.png');
Irgb_1 = imread('Data/new_view/im1.png');

% Resize to speed up debbuging
% Irgb_0 =imresize(Irgb_0, 0.1);
% Irgb_1 =imresize(Irgb_1, 0.1);

im_0 = im2double(rgb2gray(Irgb_0));
im_1 = im2double(rgb2gray(Irgb_1));
[h,w] = size(im_0);

% Load disparity maps (ground truth)
di_0 = parsePfm('Data/new_view/disp0.pfm');
di_1 = parsePfm('Data/new_view/disp1.pfm');

% Resize to speed up debbuging
% di_0 = imresize(di_0, 0.1);
% di_1 = imresize(di_1, 0.1);

%% Look for disparity consistency 
% apply some transf to disparity to match same object on both images ?
disparity_th = 20;
di_consistent_on = abs(di_0 - di_1) < disparity_th;
di_inconsistent_on = abs(di_0 - di_1) >= disparity_th;
% figure; imshow(di_consistent_on)
% Keep only consistent values (set others to 0)
im_0_true = di_0 .* di_consistent_on;
figure; imshow(di_consistent_on)

%% Show data
figure(1); 
subplot(2,2,1); imshow(im_0);
subplot(2,2,2); imshow(im_1);
subplot(2,2,3); imshow(di_0, []);
subplot(2,2,4); imshow(di_1, []);

%% Compute SIFT keypoints
[points_0, desc_0] = sift(im_0, 'Threshold', 0.01);
[points_1, desc_1] = sift(im_1, 'Threshold', 0.01);

%% Match SIFT keypoints between Image 0 and Image 1
matches = siftmatch(desc_0, desc_1);
% Show matches
% figure;
% plotmatches(im0, im1, points_0(1:2,:), points_1(1:2,:), ...
%             matches, 'Stacking', 'v');

% keypoints are the matches' coord (homogeneous)
p0 = [points_0(1:2, matches(1,:)); ones(1, length(matches))];
p1 = [points_1(1:2, matches(2,:)); ones(1, length(matches))];

% Estimate fundamental matrix
[F, inliers] = ransac_fundamental_matrix(p0, p1, 2); 

% Show inliers
% figure;
% plotmatches(im0, im1, points_0(1:2,:), points_1(1:2,:), ...
%             matches(:,inliers), 'Stacking', 'v');
% title('Inliers');
% vgg_gui_F(Irgb0, Irgb1, F');

%% Define projection matrices
Cx = -1; % x-axis translation of camera 1 wr camera 0
PR_0 = [eye(3) zeros(3,1)];
PR_1 = [eye(3); Cx 0 0]';

% Compute projection matrix of Image S
s = 0.5;
PR_s = (1-s)*PR_0 + s*PR_1;

% Compute matched points on Image S
ps = (1-s)*p0 + s*p1;

%% Show how keypoints are moved from both real images
% It can be seen the objects are moved to the left on im0 and the right 
% of im1. Indicating the transformed image is inbetween both given cameras.
ps_e = euclid(ps);

figure; 
subplot(1,2,1);
imshow(im_0,[]); hold on;
scatter(p0(1,:), p0(2,:),'g.');
scatter(ps_e(1,:), ps_e(2,:),'r.'); hold off;
title('Translation from camera0')
legend('from camera0', 'new view')
subplot(1,2,2);
imshow(im_1,[]); hold on;
scatter(p1(1,:), p1(2,:),'g.');
scatter(ps_e(1,:), ps_e(2,:),'r.'); hold off;
title('Translation from camera1')
legend('from camera1', 'new view')

%% Find matching points pixel-by-pixel
% Build points on im0
[x,y] = find(ones(size(im_0,1), size(im_0,2)));
p_0 = [x y]';
% Init p1
p_1 = p_0 .* 0;
l = 0;
for j=1:size(im_0,1)
% search foreach scanline
  for k=1:size(im_0,2)
    l = l + 1;
    % Get pixel from image 0
    pixel_0 = im_0(j, k);
    % Find the most similar one on image 1 (same row)
    [I_SSD,I_NCC,Idata] = template_matching(pixel_0, im_1(j,:));
    sim_xy = find(I_SSD == max(I_SSD(:)));
    % Build points
    p_1(:,l) = [k; sim_xy(1)];
    % similiar pixel can be computed as: pixel_1 = im1(i, sim_xy(1));
  end
end

%% Find matching points block-by-block
% Build points on im0
[x,y] = find(ones(size(im_0,1), size(im_0,2)));
p_0 = [x y]';
% Init p1
p_1 = p_0 .* 0;
l = 0;
block_size = 20;

% Pad image to avoid border control
pad_size_y = block_size - mod(size(im_0, 1), block_size);
pad_size_x = block_size - mod(size(im_0, 2), block_size);
im_0_padded = padarray(im_0, [pad_size_x pad_size_y], 0, 'post');
im_1_padded = padarray(im_1, [pad_size_x pad_size_y], 0, 'post');

size(im_0), size(im_0_padded)
% figure; imshow(im_0)
% figure; imshow(im_0_padded)

% foreach block_size-rows (20 rows)
figure;
for j=1:block_size:(size(im_0,1)-block_size)
    % Get some rows from image 0
    image_0 = im_0_padded(j:(j+block_size-1),:);
%     figure; imshow(image_0);
    image_1 = im_1_padded(j:(j+block_size-1),:);
%     figure; imshow(image_1);
  
  for k=1:block_size:(size(im_0,2)-block_size)
    % get one block of 20 columns (20x20 patch)
%     [k, (k+block_size-1)]
    template_0 = im_0_padded(:,k:(k+block_size-1));
%     size(template_0)
%     size(image_1)
    [I_SSD,I_NCC,Idata] = template_matching(template_0, image_1);
%     fun = @(block_struct) template_matching(block_struct.data, block20rows_1);
%     result = blockproc(im_0, [block_size block_size], fun);
%     size(result)
    % Find the most similar one on image 1 (same row)
%     [I_SSD,I_NCC,Idata] = template_matching(block20rows_0, im_1(j,:));
    I_NCC = sum(I_NCC, 1); % acumulate SSD of each column (sum diff of all rows)
    norm_I_NCC =  sum(I_NCC,1)/block_size;
    % Remove padding
    norm_I_NCC = norm_I_NCC(1:stopper);
    plot(norm_I_NCC); hold on;
    stopper = size(im_1,2) - block_size;
    sim_xy = find(norm_I_NCC == max(norm_I_NCC(:)));
%     sim_xy
    % Build points
    p_1(:,l) = [k; sim_xy(1)];
    % similiar pixel can be computed as: block_1 = ... ;
%     break
  end
end
axis tight;
hold off;

%% From the position difference on x-axis, compute disparity
% use p0 and p1
disparity = p_0(1,:) - p_1(1,:);
disparity = disparity - min(disparity(:));
disparity = disparity/max(disparity(:));
disparity_map = im_0 * 0;
for pos = 1:size(disparity,2)
  x = p_0(1,pos);
  y = p_0(2,pos);
  disparity_map(x,y) = disparity(pos);
end

imshow(disparity_map,[])
%% Depth is inv proportional to disparity

%% Use disparity or depth to interpolate new view

%% Replace disparity GT by estimating it using stereo_computation()

%% other stuff
% %% Epipolar rectification (external processing)
% % http://demo.ipol.im/demo/m_quasi_euclidean_epipolar_rectification/result?key=9F251CF6634AABDBED6EF711FAA6CC78
% % Input: (im0s.png, im1s.png)  
% % Output: (rec0.png, rec1.png, correspondences.txt) 
% % Log:
% %   best matching found:  183 points  log(nfa)=-420.636  (500 iterations)
% %   F= [ -2.20981e-12 -4.99157e-11 2.39569e-06; 7.85206e-10 1.8295e-10 -0.000789925; -2.84774e-06 0.00078936 0.000384152 ]
% %   Geometric error threshold: 0.845051
% %   LM iterations: 3 f=1823
% %   K_left: [ 1823 0 541.557; 0 1823 371.5; 0 0 1 ]
% %   K_right: [ 1823 0 540.864; 0 1823 371.5; 0 0 1 ]
% %   Initial rectification error: 0.198191 pix
% %   Final rectification error: 0.192938 pix
% %   Disparity: -208 -13
% rec0 = imread('Data/new_view/rec0.png');
% rec1 = imread('Data/new_view/rec1.png');
% 
% % Correspondences as column vectors as [x1,y1,x2,y2]' coordinates
% corresp = load('Data/new_view/correspondences.txt')';
% 
% figure(2);
% imshow(rec0);
% hold on; scatter(corresp(1,:), corresp(2,:)); hold off;
% figure(3);
% imshow(rec1);
% hold on; scatter(corresp(3,:), corresp(4,:)); hold off;
% 
% %% Following paper
% corresp_0 = homog(corresp(1:2,:));
% corresp_1 = homog(corresp(3:4,:));
% 
% F = [ -2.20981e-12 -4.99157e-11 2.39569e-06; ...
%       7.85206e-10 1.8295e-10 -0.000789925; ...
%       -2.84774e-06 0.00078936 0.000384152 ];
%     
% K_left = [ 1823 0 541.557; 0 1823 371.5; 0 0 1 ];
% K_right = [ 1823 0 540.864; 0 1823 371.5; 0 0 1 ];
% 
% % It should be something like this (only x-movement)
% P1 = eye(3,4);
% R = eye(3,3);
% t = [10 0 1]';
% P2 = [R t];
% 
% PI0 = [K_left K_left*C0];
% s = 0.5
% (1-s) * p0 + s * p0
% 
% %% Computing disparity
% % x_minus_xp = corr(1,:)-corr(3,:);
% 
% %% Select a template in image 1
% % w = 20;
% % row = 250;
% % column = 350;
% % srow1 = rec0(row:(row+w), :, :);
% % srow2 = rec1(row:(row+w), :, :);
% template = im0(560:659, 1530:1669,:);
% [h,w] = size(template);
% image = im1(560:660, :,:);
% % patch1 = rec0(row:(row+w), column:(column+w), :);
% figure(4);
% subplot(2,1,1);
% imshow(template)
% title('Template');
% subplot(2,1,2);
% imshow(image)
% title('Image to look in');
% 
% %% Do template matching over image 2
% [I_SSD,I_NCC, ~] = template_matching(template, image);
% % [err, minErrorPosition] = patch_scan(patch1, patch2);
% % mep = minErrorPosition;
% % figure;
% % plot(err);
% % p2 = rec1(row:(row+w), (mep-w/2):(mep+w/2-1), :);
% 
% %% Plot similarity score of template over image
% figure(5); 
% subplot(2,1,1); plot(I_SSD(:)); title('SSD distance');
% subplot(2,1,2); plot(I_NCC(:)); title('NCC distance');
% [x,y] = find(I_NCC==max(I_NCC(:)))
% die
% %% Get region on image 2 that corresponds to template
% patch2 = rec1((x-w/2):(x+w/2-1), (y-w/2):(y+w/2-1), :);
% % p2 = rec1((x-w/2):(x+w/2-1), (y-w/2):(y+w/2-1), :);
% figure(6);
% subplot(4,1,4);
% imshow(p1);
% 
% % subplot(1,3,1);
% % imshow(patch1);
% % subplot(1,3,2)
% % imshow(p2);
% % subplot(1,3,3)
% 
% 
% % I1_hat = H1 * homog(I1);
% % I2_hat = H2 * I2;
% % X0(:,:,1) = I1;
% % X0(:,:,2) = disp0_gt;
% % p = (1-s) .* [x,y] + s .* [x-d(x,y), y];
% 
% %%
% % Find maximum response
%  I = im2double(imread('lena.png'));
% % Template of Eye Lena
%  T=I(124:140,124:140,:);
%  T=I(250:280 ,250:290,:);
% % Calculate SSD and NCC between Template and Image
%  [I_SSD,I_NCC]=template_matching(T,I);
% % Find maximum correspondence in I_SDD image
%  [x,y]=find(I_SSD==max(I_SSD(:)));
%  [x,y]=find(I_NCC==max(I_NCC(:)));
% % Show result
%  figure, 
%  subplot(2,2,1), imshow(I); hold on; plot(y,x,'r*'); title('Result')
%  subplot(2,2,2), imshow(T); title('The eye template');
%  subplot(2,2,3), imshow(I_SSD); title('SSD Matching');
%  subplot(2,2,4), imshow(I_NCC); title('Normalized-CC');