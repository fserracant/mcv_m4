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
