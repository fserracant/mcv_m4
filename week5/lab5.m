%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 5: Reconstruction from uncalibrated viewas


addpath('sift'); % ToDo: change 'sift' to the correct path where you have the sift functions


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. Create synthetic data

% 3D points
X = [];
Xi = [0 0 10]'; % middle facade
X = [X Xi];
Xi = [0 20 10]';
X = [X Xi];
Xi = [0 20 0]';
X = [X Xi];
Xi = [0 0 0]';
X = [X Xi];
Xi = [10 20 10]'; % right facade
X = [X Xi];
Xi = [10 20 0]';
X = [X Xi];
Xi = [30 0 10]';  % left facade
X = [X Xi];
Xi = [30 0 0]'; 
X = [X Xi];
Xi = [0 5 8]';   % middle squared window
X = [X Xi];
Xi = [0 8 8]';
X = [X Xi];
Xi = [0 5 5]';
X = [X Xi];
Xi = [0 8 5]';
X = [X Xi];
Xi = [0 12 5]';  % middle rectangular window
X = [X Xi];
Xi = [0 18 5]';
X = [X Xi];
Xi = [0 12 2]';
X = [X Xi];
Xi = [0 18 2]';
X = [X Xi];
Xi = [3 20 7]';   % left squared window
X = [X Xi];
Xi = [7 20 7]';
X = [X Xi];
Xi = [3 20 3]';
X = [X Xi];
Xi = [7 20 3]';
X = [X Xi];
Xi = [5 0 7]';   % right rectangular window
X = [X Xi];
Xi = [25 0 7]';
X = [X Xi];
Xi = [5 0 3]';
X = [X Xi];
Xi = [25 0 3]';
X = [X Xi];


% cameras

K = [709 0 450; 0 709 300; 0 0 1];
Rz = [cos(0.88*pi/2) -sin(0.88*pi/2) 0; sin(0.88*pi/2) cos(0.88*pi/2) 0; 0 0 1];
Ry = [cos(0.88*pi/2) 0 sin(0.88*pi/2); 0 1 0; -sin(0.88*pi/2) 0 cos(0.88*pi/2)];
R1 = Rz*Ry;
t1 = -R1*[40; 10; 5];

Rz = [cos(0.8*pi/2) -sin(0.8*pi/2) 0; sin(0.8*pi/2) cos(0.8*pi/2) 0; 0 0 1];
Ry = [cos(0.88*pi/2) 0 sin(0.88*pi/2); 0 1 0; -sin(0.88*pi/2) 0 cos(0.88*pi/2)];
Rx = [1 0 0; 0 cos(-0.15) -sin(-0.15); 0 sin(-0.15) cos(-0.15)];
R2 = Rx*Rz*Ry;
t2 = -R2*[45; 15; 5];

P1 = zeros(3,4);
P1(1:3,1:3) = R1;
P1(:,4) = t1;
P1 = K*P1;

P2 = zeros(3,4);
P2(1:3,1:3) = R2;
P2(:,4) = t2;
P2 = K*P2;

w = 900;
h = 600;

% visualize as point cloud
figure; hold on;
plot_camera2(P1,w,h);
plot_camera2(P2,w,h);
for i = 1:length(X)
    scatter3(X(1,i), X(2,i), X(3,i), 5^2, [0.5 0.5 0.5], 'filled');
end;
axis equal;
axis vis3d;

%% visualize as lines
figure;
hold on;
X1 = X(:,1); X2 = X(:,2); X3 = X(:,3); X4 = X(:,4);
plot3([X1(1) X2(1)], [X1(2) X2(2)], [X1(3) X2(3)]);
plot3([X3(1) X4(1)], [X3(2) X4(2)], [X3(3) X4(3)]);
X5 = X(:,5); X6 = X(:,6); X7 = X2; X8 = X3;
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = X(:,7); X6 = X(:,8); X7 = X1; X8 = X4;
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = X(:,9); X6 = X(:,10); X7 = X(:,11); X8 = X(:,12);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = X(:,13); X6 = X(:,14); X7 = X(:,15); X8 = X(:,16);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = X(:,17); X6 = X(:,18); X7 = X(:,19); X8 = X(:,20);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = X(:,21); X6 = X(:,22); X7 = X(:,23); X8 = X(:,24);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
axis vis3d
axis equal
plot_camera2(P1,w,h);
plot_camera2(P2,w,h);

%% Create homogeneous coordinates

% homogeneous 3D coordinates
Xh=[X; ones(1,length(X))];

% homogeneous 2D coordinates
x1 = P1*Xh;
x1(1,:) = x1(1,:)./x1(3,:);
x1(2,:) = x1(2,:)./x1(3,:);
x1(3,:) = x1(3,:)./x1(3,:);
x2 = P2*Xh;
x2(1,:) = x2(1,:)./x2(3,:);
x2(2,:) = x2(2,:)./x2(3,:);
x2(3,:) = x2(3,:)./x2(3,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Projective reconstruction (synthetic data)

% ToDo: create the function 'factorization_method' that computes a
% projective reconstruction with the factorization method of Sturm and
% Triggs '1996
% This function returns an estimate of:
%       Pproj: 3*Ncam x 4 matrix containing the camera matrices
%       Xproj: 4 x Npoints matrix of homogeneous coordinates of 3D points
% 
% As a convergence criterion you may compute the Euclidean
% distance (d) between data points and projected points in both images 
% and stop when (abs(d - d_old)/d) < 0.1 where d_old is the distance
% in the previous iteration.

%% Check projected points (estimated and data points)

for i=1:2
    x_proj{i} = euclid(Pproj(3*i-2:3*i, :)*Xproj);
end
x_d{1} = euclid(P1*Xh);
x_d{2} = euclid(P2*Xh);

% image 1
figure;
hold on
plot(x_d{1}(1,:),x_d{1}(2,:),'r*');
plot(x_proj{1}(1,:),x_proj{1}(2,:),'bo');
axis equal

% image 2
figure;
hold on
plot(x_d{2}(1,:),x_d{2}(2,:),'r*');
plot(x_proj{2}(1,:),x_proj{2}(2,:),'bo');


%% Visualize projective reconstruction
Xaux(1,:) = Xproj(1,:)./Xproj(4,:);
Xaux(2,:) = Xproj(2,:)./Xproj(4,:);
Xaux(3,:) = Xproj(3,:)./Xproj(4,:);
X=Xaux;

figure;
hold on;
X1 = X(:,1); X2 = X(:,2); X3 = X(:,3); X4 = X(:,4);
plot3([X1(1) X2(1)], [X1(2) X2(2)], [X1(3) X2(3)]);
plot3([X3(1) X4(1)], [X3(2) X4(2)], [X3(3) X4(3)]);
X5 = X(:,5); X6 = X(:,6); X7 = X2; X8 = X3;
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = X(:,7); X6 = X(:,8); X7 = X1; X8 = X4;
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = X(:,9); X6 = X(:,10); X7 = X(:,11); X8 = X(:,12);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = X(:,13); X6 = X(:,14); X7 = X(:,15); X8 = X(:,16);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = X(:,17); X6 = X(:,18); X7 = X(:,19); X8 = X(:,20);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = X(:,21); X6 = X(:,22); X7 = X(:,23); X8 = X(:,24);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
axis vis3d
axis equal



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Affine reconstruction (synthetic data)

% ToDo: create the function 'vanishing_point' that computes the vanishing
% point formed by the line that joins points xo1 and xf1 and the line 
% that joins points x02 and xf2
%
% [v1] = vanishing_point(xo1, xf1, xo2, xf2)

% Compute the vanishing points in each image
v1 = vanishing_point(x1(:,21),x1(:,22),x1(:,23),x1(:,24));
v2 = vanishing_point(x1(:,21),x1(:,23),x1(:,22),x1(:,24));
v3 = vanishing_point(x1(:,1),x1(:,2),x1(:,4),x1(:,3));

v1p = vanishing_point(x2(:,21),x2(:,22),x2(:,23),x2(:,24));
v2p = vanishing_point(x2(:,21),x2(:,23),x2(:,22),x2(:,24));
v3p = vanishing_point(x2(:,1),x2(:,2),x2(:,4),x2(:,3));

% ToDo: use the vanishing points to compute the matrix Hp that 
%       upgrades the projective reconstruction to an affine reconstruction


%% check results

Xa = euclid(Hp*Xproj);
figure;
hold on;
X1 = Xa(:,1); X2 = Xa(:,2); X3 = Xa(:,3); X4 = Xa(:,4);
plot3([X1(1) X2(1)], [X1(2) X2(2)], [X1(3) X2(3)]);
plot3([X3(1) X4(1)], [X3(2) X4(2)], [X3(3) X4(3)]);
X5 = Xa(:,5); X6 = Xa(:,6); X7 = X2; X8 = X3;
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = Xa(:,7); X6 = Xa(:,8); X7 = X1; X8 = X4;
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = Xa(:,9); X6 = Xa(:,10); X7 = Xa(:,11); X8 = Xa(:,12);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = Xa(:,13); X6 = Xa(:,14); X7 = Xa(:,15); X8 = Xa(:,16);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = Xa(:,17); X6 = Xa(:,18); X7 = Xa(:,19); X8 = Xa(:,20);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = Xa(:,21); X6 = Xa(:,22); X7 = Xa(:,23); X8 = Xa(:,24);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
axis vis3d
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric reconstruction (synthetic data)

% ToDo: compute the matrix Ha that 
%       upgrades the projective reconstruction to an affine reconstruction
% Use the following vanishing points given by three pair of orthogonal lines
% and assume that the skew factor is zero and that pixels are square

v1 = vanishing_point(x1(:,2),x1(:,5),x1(:,3),x1(:,6));
v2 = vanishing_point(x1(:,1),x1(:,2),x1(:,3),x1(:,4));
v3 = vanishing_point(x1(:,1),x1(:,4),x1(:,2),x1(:,3));

%% check results

Xa = euclid(Ha*Hp*Xproj);
figure;
hold on;
X1 = Xa(:,1); X2 = Xa(:,2); X3 = Xa(:,3); X4 = Xa(:,4);
plot3([X1(1) X2(1)], [X1(2) X2(2)], [X1(3) X2(3)]);
plot3([X3(1) X4(1)], [X3(2) X4(2)], [X3(3) X4(3)]);
X5 = Xa(:,5); X6 = Xa(:,6); X7 = X2; X8 = X3;
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = Xa(:,7); X6 = Xa(:,8); X7 = X1; X8 = X4;
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = Xa(:,9); X6 = Xa(:,10); X7 = Xa(:,11); X8 = Xa(:,12);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = Xa(:,13); X6 = Xa(:,14); X7 = Xa(:,15); X8 = Xa(:,16);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = Xa(:,17); X6 = Xa(:,18); X7 = Xa(:,19); X8 = Xa(:,20);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
X5 = Xa(:,21); X6 = Xa(:,22); X7 = Xa(:,23); X8 = Xa(:,24);
plot3([X5(1) X6(1)], [X5(2) X6(2)], [X5(3) X6(3)]);
plot3([X7(1) X8(1)], [X7(2) X8(2)], [X7(3) X8(3)]);
plot3([X5(1) X7(1)], [X5(2) X7(2)], [X5(3) X7(3)]);
plot3([X6(1) X8(1)], [X6(2) X8(2)], [X6(3) X8(3)]);
axis vis3d
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Projective reconstruction (real data)

%% read images
Irgb{1} = double(imread('Data/0000_s.png'))/255;
Irgb{2} = double(imread('Data/0001_s.png'))/255;

I{1} = sum(Irgb{1}, 3) / 3; 
I{2} = sum(Irgb{2}, 3) / 3;

Ncam = length(I);

% ToDo: compute a projective reconstruction using the factorization method

% ToDo: show the data points (image correspondences) and the projected
% points (of the reconstructed 3D points) in images 1 and 2. Reuse the code
% in section 'Check projected points' (synthetic experiment).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Affine reconstruction (real data)

% ToDo: compute the matrix Hp that updates the projective reconstruction
% to an affine one
%
% You may use the vanishing points given by function 'detect_vps' that 
% implements the method presented in Lezama et al. CVPR 2014
% (http://dev.ipol.im/~jlezama/vanishing_points/)

% This is an example on how to obtain the vanishing points (VPs) from three
% orthogonal lines in image 1

img_in =  'Data/0000_s.png'; % input image
folder_out = '.'; % output folder
manhattan = 1;
acceleration = 0;
focal_ratio = 1;
params.PRINT = 1;
params.PLOT = 1;
[horizon, VPs] = detect_vps(img_in, folder_out, manhattan, acceleration, focal_ratio, params);


%% Visualize the result

% x1m are the data points in image 1
% Xm are the reconstructed 3D points (projective reconstruction)

r = interp2(double(Irgb{1}(:,:,1)), x1m(1,:), x1m(2,:));
g = interp2(double(Irgb{1}(:,:,2)), x1m(1,:), x1m(2,:));
b = interp2(double(Irgb{1}(:,:,3)), x1m(1,:), x1m(2,:));
Xe = euclid(Hp*Xm);
figure; hold on;
[w,h] = size(I{1});
for i = 1:length(Xe)
    scatter3(Xe(1,i), Xe(2,i), Xe(3,i), 2^2, [r(i) g(i) b(i)], 'filled');
end;
axis equal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Metric reconstruction (real data)

% ToDo: compute the matrix Ha that updates the affine reconstruction
% to a metric one and visualize the result in 3D as in the previous section

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. OPTIONAL: Projective reconstruction from two views

% ToDo: compute a projective reconstruction from the same two views 
% by computing two possible projection matrices from the fundamental matrix
% and one of the epipoles.
% Then update the reconstruction to affine and metric as before (reuse the code).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8. OPTIONAL: Projective reconstruction from more than two views

% ToDo: extend the function that computes the projective reconstruction 
% with the factorization method to the case of three views. You may use 
% the additional image '0002_s.png'
% Then update the reconstruction to affine and metric.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 9. OPTIONAL: Any other improvement you may incorporate 

% Add a 4th view, incorporate new 3D points by triangulation, 
% incorporate new views by resectioning, 
% apply any kind of processing on the point cloud, ...)

