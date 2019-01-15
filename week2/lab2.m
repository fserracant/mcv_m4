%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 2: Image mosaics

addpath('../sift');
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Compute image correspondences

%% Open images

fprintf(1, 'Opening images\n');
imargb = imread('Data/llanes/llanes_a.jpg');
imbrgb = imread('Data/llanes/llanes_b.jpg');
imcrgb = imread('Data/llanes/llanes_c.jpg');

% imargb = imread('Data/castle_int/0016_s.png');
% imbrgb = imread('Data/castle_int/0015_s.png');
% imcrgb = imread('Data/castle_int/0014_s.png');
% 
% imargb = imread('Data/aerial/site13/frame00000.png');
% imbrgb = imread('Data/aerial/site13/frame00002.png');
% imcrgb = imread('Data/aerial/site13/frame00003.png');
 
ima = sum(double(imargb), 3) / 3 / 255;
imb = sum(double(imbrgb), 3) / 3 / 255;
imc = sum(double(imcrgb), 3) / 3 / 255;

% imargb = double(imread('Data/aerial/site22/frame_00001.tif'))/255;
% imbrgb = double(imread('Data/aerial/site22/frame_00018.tif'))/255;
% imcrgb = double(imread('Data/aerial/site22/frame_00030.tif'))/255;
% ima = imargb;
% imb = imbrgb;
% imc = imcrgb;

%% Compute SIFT keypoints
% From 'sift' doc: Note that the X,Y center coordinates are (0,0) based, 
% contrary to the standard MATLAB convention that uses (1,1) as the 
% top-left image coordinate
fprintf(1, 'Computing SIFT keypoints\n');
[points_a, desc_a] = sift(ima, 'Threshold', 0.01);
[points_b, desc_b] = sift(imb, 'Threshold', 0.01);
[points_c, desc_c] = sift(imc, 'Threshold', 0.01);

figure('Name', 'Image A SIFT descriptors');
imshow(imargb);%image(imargb)
hold on;
plot(points_a(1,:), points_a(2,:),'+y');
figure('Name', 'Image B SIFT descriptors');
imshow(imbrgb);%image(imbrgb);
hold on;
plot(points_b(1,:), points_b(2,:),'+y');
figure('Name', 'Image C SIFT descriptors');
imshow(imcrgb);%image(imcrgb);
hold on;
plot(points_c(1,:), points_c(2,:),'+y');

%% Match SIFT keypoints 

fprintf(1, 'Matching SIFT keypoints\n');
% between a and b
matches_ab = siftmatch(desc_a, desc_b);
figure('Name', 'Matches A-B');
plotmatches(ima, imb, points_a(1:2,:), points_b(1:2,:), matches_ab, 'Stacking', 'v');

% between b and c
matches_bc = siftmatch(desc_b, desc_c);
figure('Name', 'Matches B-C');
plotmatches(imb, imc, points_b(1:2,:), points_c(1:2,:), matches_bc, 'Stacking', 'v');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;  % To control the amount of total figures being displayed!

%% 2. Compute the homography (DLT algorithm) between image pairs

%% Compute homography (normalized DLT) between a and b, play with the homography
fprintf(1, 'Computing homography for images a and b (normalized DLT)\n');
th = 3;
xab_a = [points_a(1:2, matches_ab(1,:)); ones(1, length(matches_ab))];
xab_b = [points_b(1:2, matches_ab(2,:)); ones(1, length(matches_ab))];
[Hab, inliers_ab] = ransac_homography_adaptive_loop(xab_a, xab_b, th, 1000); % ToDo: complete this function

figure('Name', 'Matches A-B (only inliers)');
plotmatches(ima, imb, points_a(1:2,:), points_b(1:2,:), ...
    matches_ab(:,inliers_ab), 'Stacking', 'v');

vgg_gui_H(imargb, imbrgb, Hab);


%% Compute homography (normalized DLT) between b and c, play with the homography
fprintf(1, 'Computing homography for images b and c (normalized DLT)\n');
xbc_b = [points_b(1:2, matches_bc(1,:)); ones(1, length(matches_bc))];
xbc_c = [points_c(1:2, matches_bc(2,:)); ones(1, length(matches_bc))];
[Hbc, inliers_bc] = ransac_homography_adaptive_loop(xbc_b, xbc_c, th, 1000); 

figure('Name', 'Matches B-C (only inliers)');
plotmatches(imb, imc, points_b(1:2,:), points_c(1:2,:), ...
    matches_bc(:,inliers_bc), 'Stacking', 'v');

vgg_gui_H(imbrgb, imcrgb, Hbc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;  % To control the amount of total figures being displayed!

%% 3. Build the mosaic
fprintf(1, 'Building the Mosaic...\n');
corners = [-400 1200 -100 650];
iwb = apply_H_v2(imbrgb, eye(3), corners);   % ToDo: complete the call to the function
iwa = apply_H_v2(imargb, Hab, corners);    % ToDo: complete the call to the function
iwc = apply_H_v2(imcrgb, inv(Hbc), corners);    % ToDo: complete the call to the function

figure('Name', 'Mosaic ''llanes''');
imshow(max(iwc, max(iwb, iwa)));%image(max(iwc, max(iwb, iwa)));axis off;
title('Llanes Mosaic A-B-C');

% In order to reuse code, we created a function compute_mosaic
% of three images with the code from points 1 and 2

% % ToDo: compute the mosaic with castle_int images
im1 = imread('Data/castle_int/0016_s.png');
im2 = imread('Data/castle_int/0015_s.png');
im3 = imread('Data/castle_int/0014_s.png');
mosaic = compute_mosaic(im1, im2, im3, 0.01, true, false);
figure('Name', 'Mosaic ''castle_int''');
imshow(mosaic);
title('Castle Mosaic A-B-C');

% % ToDo: compute the mosaic with aerial images set 13
im1 = imread('Data/aerial/site13/frame00000.png');
im2 = imread('Data/aerial/site13/frame00002.png');
im3 = imread('Data/aerial/site13/frame00003.png');
figure('Name', 'Mosaic ''site13''');
imshow(compute_mosaic(im1, im2, im3, 0.01, true, false));
title('Site13 Mosaic A-B-C');
% 
% % ToDo: compute the mosaic with aerial images set 22
im1 = double(imread('Data/aerial/site22/frame_00001.tif'));
im2 = double(imread('Data/aerial/site22/frame_00018.tif'));
im3 = double(imread('Data/aerial/site22/frame_00030.tif'));
figure('Name', 'Mosaic ''site22''');
imshow(compute_mosaic(im1, im2, im3, 0.05, false, false));
title('Site22 Mosaic A-B-C');

% ToDo: comment the results in every of the four cases: say why it works or
%       does not work

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;  % To control the amount of total figures being displayed!

%% 4. Refine the homography with the Gold Standard algorithm

% Homography ab
fprintf(1, 'Refining homography for images a and b (Gold Standard algorithm)\n');
% x are the xy-coordinates of points in the image 'a' that made match with image 'b'
% xp are the xy-oordinates of points in the image 'b' that made match with image 'a'
x = points_a(1:2, matches_ab(1,inliers_ab));  %ToDo: set the non-homogeneous point coordinates of the 
xp = points_b(1:2, matches_ab(2,inliers_ab)); %      point correspondences we will refine with the geometric method
Xobs = [ x(:) ; xp(:) ];     % The column vector of observed values (x and x')
P0 = [ Hab(:) ; x(:) ];      % The parameters or independent variables

% ToDo: create this function that we need to pass to the lsqnonlin function
% NOTE: gs_errfunction should return E(X) and not the sum-of-squares E=sum(E(X).^2)) that we want to minimize. 
% (E(X) is summed and squared implicitly in the lsqnonlin algorithm.) 
Y_initial = gs_errfunction( P0, Xobs );
err_initial = sum( sum( Y_initial.^2 ));

options = optimset('Algorithm', 'levenberg-marquardt', 'Diagnostics', 'on', 'Display', 'iter');
P = lsqnonlin(@(t) gs_errfunction(t, Xobs), P0, [], [], options);

Hab_r = reshape( P(1:9), 3, 3 );
f = gs_errfunction( P, Xobs ); % lsqnonlin does not return f
err_final = sum( sum( f.^2 ));

% we show the geometric error before and after the refinement
fprintf(1, 'Gold standard reproj error initial %f, final %f\n', err_initial, err_final);
fprintf(1, 'Error has been reduced %f times\n', err_initial/err_final);


%% See differences in the keypoint locations

% ToDo: compute the points xhat and xhatp which are the correspondences
% returned by the refinement with the Gold Standard algorithm

xhat = [reshape(P(9+1:end),2,[]);ones(1,length(x))];
xhatp = Hab_r*xhat;

figure('Name', 'ImA: Keypoints vs refined GS')
imshow(imargb);%image(imargb);
hold on;
plot(x(1,:), x(2,:),'oy');
plot(xhat(1,:), xhat(2,:),'+c');
hold off;

figure('Name', 'ImB: Keypoints vs refined GS')
imshow(imbrgb);%image(imbrgb);
hold on;
plot(xp(1,:), xp(2,:),'oy');
plot(xhatp(1,:), xhatp(2,:),'+c');
hold off;

%%  Homography bc

fprintf(1, 'Refineing homography for images b and c (Gold Standard algorithm)\n');
% ToDo: refine the homography bc with the Gold Standard algorithm
% x are the xy-coordinates of points in the image 'b' that made match with image 'c'
% xp are the xy-oordinates of points in the image 'c' that made match with image 'b'
x = points_b(1:2, matches_bc(1,inliers_bc));  %ToDo: set the non-homogeneous point coordinates of the 
xp = points_c(1:2, matches_bc(2,inliers_bc)); %      point correspondences we will refine with the geometric method
Xobs = [ x(:) ; xp(:) ];     % The column vector of observed values (x and x')
P0 = [ Hbc(:) ; x(:) ];      % The parameters or independent variables

Y_initial = gs_errfunction(P0, Xobs); % ToDo: create this function that we need to pass to the lsqnonlin function
% NOTE: gs_errfunction should return E(X) and not the sum-of-squares E=sum(E(X).^2)) that we want to minimize. 
% (E(X) is summed and squared implicitly in the lsqnonlin algorithm.) 
err_initial = sum( sum( Y_initial.^2 ));

options = optimset('Algorithm', 'levenberg-marquardt', 'Diagnostics', 'on', 'Display', 'iter');
P = lsqnonlin(@(t) gs_errfunction(t, Xobs), P0, [], [], options);

Hbc_r = reshape( P(1:9), 3, 3 );
f = gs_errfunction( P, Xobs ); % lsqnonlin does not return f
err_final = sum( sum( f.^2 ));

% we show the geometric error before and after the refinement
fprintf(1, 'Gold standard reproj error initial %f, final %f\n', err_initial, err_final);
fprintf(1, 'Error has been reduced %f times\n', err_initial/err_final);

%% See differences in the keypoint locations

% ToDo: compute the points xhat and xhatp which are the correspondences
% returned by the refinement with the Gold Standard algorithm

xhat = [reshape(P(9+1:end),2,[]);ones(1,length(x))];
xhatp = Hbc_r*xhat;

figure('Name', 'ImB: Keypoints vs refined GS')
imshow(imbrgb);%image(imbrgb);
hold on;
plot(x(1,:), x(2,:),'+y');
plot(xhat(1,:), xhat(2,:),'+c');

figure('Name', 'ImC: Keypoints vs refined GS')
imshow(imcrgb);%image(imcrgb);
hold on;
plot(xp(1,:), xp(2,:),'+y');
plot(xhatp(1,:), xhatp(2,:),'+c');

%% Build mosaic
fprintf(1, 'Building the Mosaic\n');
corners = [-400 1200 -100 650];
iwb = apply_H_v2(imbrgb, eye(3), corners); % ToDo: complete the call to the function
iwa = apply_H_v2(imargb, Hab_r, corners); % ToDo: complete the call to the function
iwc = apply_H_v2(imcrgb, inv(Hbc_r), corners); % ToDo: complete the call to the function

figure('Name', 'Refined mosaic (GS)')
imshow(max(iwc, max(iwb, iwa)));%image(max(iwc, max(iwb, iwa)));axis off;
title('Mosaic A-B-C');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;  % To control the amount of total figures being displayed!

%% 5. OPTIONAL: Calibration with a planar pattern

addpath('sift');
clear;

%% Read template and images.
T     = imread('Data/calib/template.jpg');
I{1}  = imread('Data/calib/graffiti1.tif');
I{2}  = imread('Data/calib/graffiti2.tif');
I{3}  = imread('Data/calib/graffiti3.tif');
%I{4}  = imread('Data/calib/graffiti4.tif');
%I{5}  = imread('Data/calib/graffiti5.tif');
Tg = sum(double(T), 3) / 3 / 255;
Ig{1} = sum(double(I{1}), 3) / 3 / 255;
Ig{2} = sum(double(I{2}), 3) / 3 / 255;
Ig{3} = sum(double(I{3}), 3) / 3 / 255;

N = length(I);

%% Compute keypoints.
fprintf('Computing sift points in template... ');
[pointsT, descrT] = sift(Tg, 'Threshold', 0.05);
fprintf(' done\n');

points = cell(N,1);
descr = cell(N,1);
for i = 1:N
    fprintf('Computing sift points in image %d... ', i);
    [points{i}, descr{i}] = sift(Ig{i}, 'Threshold', 0.05);
    fprintf(' done\n');
end

%% Match and compute homographies.
H = cell(N,1);
for i = 1:N
    % Match against template descriptors.
    fprintf('Matching image %d... ', i);
    matches = siftmatch(descrT, descr{i});
    fprintf('done\n');

    % Fit homography and remove outliers.
    x1 = pointsT(1:2, matches(1, :));
    x2 = points{i}(1:2, matches(2, :));
    H{i} = 0;
    [H{i}, inliers] =  ransac_homography_adaptive_loop(homog(x1), homog(x2), 3, 1000);

    % Plot inliers.
    figure;
    plotmatches(Tg, Ig{i}, pointsT(1:2,:), points{i}(1:2,:), matches(:, inliers));

    % Play with the homography
    %vgg_gui_H(T, I{i}, H{i});
end

%% Compute the Image of the Absolute Conic

% w = ... % ToDo
V = [];
for i = 1:N
    
    V = [V; computeVij(H{i},1,2)';(computeVij(H{i},1,1)'-computeVij(H{i},2,2)') ];

end
 
[~, ~, Vt] = svd(V);
omega = Vt(:, end);

w = [omega(1), omega(2), omega(3);...
    omega(2), omega(4), omega(5);...
    omega(3), omega(5), omega(6)];


%% Recover the camera calibration.

% K = ... % ToDo
K = inv(chol(w));

% ToDo: in the report make some comments related to the obtained internal
%       camera parameters and also comment their relation to the image size

%% Compute camera position and orientation.
R = cell(N,1);
t = cell(N,1);
P = cell(N,1);
figure;hold;
for i = 1:N
    % ToDo: compute r1, r2, and t{i}
%     r1 = ...
%     r2 = ...
%     t{i} = ...

    r1 = K\H{i}(:,1)./(norm(K\H{i}(:,1)));
    r2 = K\H{i}(:,2)./(norm(K\H{i}(:,2)));
    t{i} = K\H{i}(:,3)./(norm(K\H{i}(:,1)));
    
    
    % Solve the scale ambiguity by forcing r1 and r2 to be unit vectors.
    s = sqrt(norm(r1) * norm(r2)) * sign(t{i}(3));
    r1 = r1 / s;
    r2 = r2 / s;
    t{i} = t{i} / s;
    R{i} = [r1, r2, cross(r1,r2)];
    
    % Ensure R is a rotation matrix
    [U, S, V] = svd(R{i});
    R{i} = U * eye(3) * V';
   
    P{i} = K * [R{i} t{i}];
    plot_camera(P{i}, 800, 600, 200);
end

% ToDo: in the report explain how the optical center is computed in the
%       provided code

[ny,nx] = size(T);
p1 = [0 0 0]';
p2 = [nx 0 0]';
p3 = [nx ny 0]';
p4 = [0 ny 0]';
% Draw planar pattern
vgg_scatter_plot([p1 p2 p3 p4 p1], 'g');
% Paint image texture
surface('XData',[0 nx; 0 nx],'YData',[0 0; 0 0],'ZData',[0 0; -ny -ny],'CData',T,'FaceColor','texturemap');
colormap(gray);
axis equal;

%% Plot a static camera with moving calibration pattern.
figure; hold;
plot_camera(K * eye(3,4), 800, 600, 200);
% ToDo: complete the call to the following function with the proper
%       coordinates of the image corners in the new reference system
for i = 1:N
    i1 =  R{i}*(p1 + t{i});
    i2 =  R{i}*(p2 + t{i});
    i3 =  R{i}*(p3 + t{i});
    i4 =  R{i}*(p4 + t{i});
    vgg_scatter_plot( [i1(1:3)   i2(1:3)   i3(1:3)   i4(1:3)   i1(1:3)], 'r');
end

%% Augmented reality: Plot some 3D points on every camera.
[Th, Tw] = size(Tg);
cube = [0 0 0; 1 0 0; 1 0 0; 1 1 0; 1 1 0; 0 1 0; 0 1 0; 0 0 0; 0 0 1; 1 0 1; 1 0 1; 1 1 1; 1 1 1; 0 1 1; 0 1 1; 0 0 1; 0 0 0; 1 0 0; 1 0 0; 1 0 1; 1 0 1; 0 0 1; 0 0 1; 0 0 0; 0 1 0; 1 1 0; 1 1 0; 1 1 1; 1 1 1; 0 1 1; 0 1 1; 0 1 0; 0 0 0; 0 1 0; 0 1 0; 0 1 1; 0 1 1; 0 0 1; 0 0 1; 0 0 0; 1 0 0; 1 1 0; 1 1 0; 1 1 1; 1 1 1; 1 0 1; 1 0 1; 1 0 0 ]';

X = (cube - .5) * Tw / 4 + repmat([Tw / 2; Th / 2; -Tw / 8], 1, length(cube));

for i = 1:N
    figure; colormap(gray);
    imagesc(Ig{i});
    hold on;
    x = euclid(P{i} * homog(X));
    vgg_scatter_plot(x, 'g');
end



%% 6. OPTIONAL: Detect the UPF logo in the two UPF images using the 
%%              DLT algorithm (folder "logos").
%%              Interpret and comment the results.
fprintf(1, 'Opening images\n');
logoUPF = imread('Data/logos/logoUPF.png');  % Load UPF logo
[src_height, src_width, ~] = size(logoUPF);

% 1st image: UPF stand
UPFstandrgb = imread('Data/logos/UPFstand.jpg');

logoUPF = sum(double(logoUPF), 3) / 3 / 255;
UPFstand = sum(double(UPFstandrgb), 3) / 3 / 255;
threshold = .000001;

% Detect logo in target image
[H_stand, dstcorners_stand] = detectLogos(UPFstand, logoUPF, threshold);

figure('Name', 'Detect logo in UPFstand');
imshow(UPFstandrgb);
hold on;
plot(dstcorners_stand(1,:), dstcorners_stand(2,:),'+y');

% Draw BB
% Corner 1-2
line([dstcorners_stand(1,1) dstcorners_stand(1,2)],...
  [dstcorners_stand(2,1) dstcorners_stand(2,2)], 'Color', 'y', 'LineStyle', '--');

% Corner 1-3
line([dstcorners_stand(1,1) dstcorners_stand(1,3)],...
  [dstcorners_stand(2,1) dstcorners_stand(2,3)], 'Color', 'y', 'LineStyle', '--');

% Corner 2-4
line([dstcorners_stand(1,2) dstcorners_stand(1,4)],...
  [dstcorners_stand(2,2) dstcorners_stand(2,4)], 'Color', 'y', 'LineStyle', '--');

% Corner 3-4
line([dstcorners_stand(1,3) dstcorners_stand(1,4)],...
  [dstcorners_stand(2,3) dstcorners_stand(2,4)], 'Color', 'y', 'LineStyle', '--');

% 2nd image: UPF building
UPFbuildingrgb = imread('Data/logos/UPFbuilding.jpg');
UPFbuilding = sum(double(UPFbuildingrgb), 3) / 3 / 255;

threshold = 0.05;

% Detect logo in target image
[H_building, dstcorners_building] = detectLogos(UPFbuilding, logoUPF, threshold);

figure('Name', 'Detect logo in UPFbuilding');
imshow(UPFbuildingrgb);
hold on;
plot(dstcorners_building(1,:), dstcorners_building(2,:),'+y');

% Draw BB
% Corner 1-2
line([dstcorners_building(1,1) dstcorners_building(1,2)],...
  [dstcorners_building(2,1) dstcorners_building(2,2)], 'Color', 'y', 'LineStyle', '--');

% Corner 1-3
line([dstcorners_building(1,1) dstcorners_building(1,3)],...
  [dstcorners_building(2,1) dstcorners_building(2,3)], 'Color', 'y', 'LineStyle', '--');

% Corner 2-4
line([dstcorners_building(1,2) dstcorners_building(1,4)],...
  [dstcorners_building(2,2) dstcorners_building(2,4)], 'Color', 'y', 'LineStyle', '--');

% Corner 3-4
line([dstcorners_building(1,3) dstcorners_building(1,4)],...
  [dstcorners_building(2,3) dstcorners_building(2,4)], 'Color', 'y', 'LineStyle', '--');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all;  % To control the amount of total figures being displayed!

%% 7. OPTIONAL: Replace the logo of the UPF by the master logo
%%              in one of the previous images using the DLT algorithm.

% We have to use the homography computed above with a scaling
% transformation to account for the difference between the source and destination
% template sizes.

% Load target logo (MsC in CVC)
master_logo = imread('Data/logos/logo_master.png');

% 1st image: UPF stand
corners = [1, size(UPFstand,2), 1, size(UPFstand,1)];
[master_logo_stand, H2_stand] = insertLogo(master_logo, logoUPF, H_stand,...
  corners);

figure('Name', 'MCV in UPFStand');
imshow(max(UPFstandrgb .* uint8(master_logo_stand < 1), master_logo_stand))

% 2nd image: UPF building
corners = [1, size(UPFbuilding,2), 1, size(UPFbuilding,1)];
[master_logo_building, H2_building] = insertLogo(master_logo, logoUPF,...
  H_building, corners);

figure('Name', 'MCV in UPFbuilding');
imshow(max(UPFbuildingrgb .* uint8(master_logo_building < 1), master_logo_building))
