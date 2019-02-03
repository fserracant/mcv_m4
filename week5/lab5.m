%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 5: Reconstruction from uncalibrated viewas


addpath('../sift'); % ToDo: change 'sift' to the correct path where you have the sift functions
addpath('../vanishing_points_v0.9');
addpath('../vanishing_points_v0.9/lib');
addpath('../vanishing_points_v0.9/mex_files');
mkdir('output');
clearvars
close all

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
% Use the t1 that generates a positive definite matrix on Part 3
% t1 = -R1*[40; 10; 5];   
t1 = -R1*[42; 5; 10];

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
title('Scene vertices as 3D point cloud')

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
title('3D scene with lines')

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

matches = [x1; x2];

% Run the Factorization method to retrieve Pproj and Xproj
[Pproj, Xproj] = factorization_method(matches, 'SturmAndTriggs');

%% Check projected points (estimated and data points)

for i=1:2
  x_proj{i} = euclid(Pproj(3*i-2:3*i, :)*Xproj);
end
x_d{1} = euclid(P1*Xh);
x_d{2} = euclid(P2*Xh);

% image 1
figure;
subplot(1,2,1);
hold on
plot(x_d{1}(1,:),x_d{1}(2,:),'r*');
plot(x_proj{1}(1,:),x_proj{1}(2,:),'bo');
axis equal
hold off
title('Projected points on View1')

% image 2
subplot(1,2,2);
hold on
plot(x_d{2}(1,:),x_d{2}(2,:),'r*');
plot(x_proj{2}(1,:),x_proj{2}(2,:),'bo');
hold off
title('Projected points on View2')

proj_error = [x_d{1} - x_proj{1} x_d{2} - x_proj{2}];
proj_error = proj_error.^2;
proj_error = sum(proj_error(:));

msj = sprintf('Error on projecting points: %g', proj_error);
disp(msj);
assert(proj_error < 1e-18, 'Error projected data points is not very low')

%% Visualize projective reconstruction
Xaux(1,:) = Xproj(1,:)./Xproj(4,:);
Xaux(2,:) = Xproj(2,:)./Xproj(4,:);
Xaux(3,:) = Xproj(3,:)./Xproj(4,:);
X=Xaux;

% Note: result has a projective ambiguity (is up to a projective
% transformation)
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
title('3D scene projective reconstructed')

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
% Compute plane at infinity as intersection of 3 lines (v1, v2, v3) for
% image 1 and, v1p, v2p and v3p in image 2.

Xv1 = triangulate(euclid(v1), euclid(v1p), Pproj(1:3,:), Pproj(4:6,:), [h,w]);
Xv2 = triangulate(euclid(v2), euclid(v2p), Pproj(1:3,:), Pproj(4:6,:), [h,w]);
Xv3 = triangulate(euclid(v3), euclid(v3p), Pproj(1:3,:), Pproj(4:6,:), [h,w]);

% Solve Ap = 0 ==> p is the null vector of A ==> SVD and pick last column
% (eigenvectors associated to lowest eigenvalue)
A = [Xv1'; Xv2'; Xv3'];
[~,~,V] = svd(A);
p = V(:,end);
p = p / p(end);

Hp = eye(4,4);
Hp(4,:) = p';

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
title('3D scene affine rectified')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric reconstruction (synthetic data)

% ToDo: compute the matrix Ha that
%       upgrades the affine reconstruction to a metric reconstruction
% Use the following vanishing points given by three pair of orthogonal lines
% and assume that the skew factor is zero and that pixels are square

v1 = vanishing_point(x1(:,2),x1(:,5),x1(:,3),x1(:,6));
v2 = vanishing_point(x1(:,1),x1(:,2),x1(:,3),x1(:,4));
v3 = vanishing_point(x1(:,1),x1(:,4),x1(:,2),x1(:,3));

w = inv(K)'*inv(K);
assert(w(1,2) == 0, 'Camera has non-zero skew')
assert(w(1,1) == w(2,2), 'Pixels are not square')

% Change nomenclature to match the one in lectures
u = v1;
v = v2;
z = v3;

A = [
  u(1)*v(1), u(1)*v(2) + u(2)*v(1), u(1)*v(3) + u(3)*v(1), u(2)*v(2), u(2)*v(3)+u(3)*v(2), u(3)*v(3);
  u(1)*z(1), u(1)*z(2) + u(2)*z(1), u(1)*z(3) + u(3)*z(1), u(2)*z(2), u(2)*z(3)+u(3)*z(2), u(3)*z(3);
  v(1)*z(1), v(1)*z(2) + v(2)*z(1), v(1)*z(3) + v(3)*z(1), v(2)*z(2), v(2)*z(3)+v(3)*z(2), v(3)*z(3);
  0 1 0 0 0 0;
  1 0 0 -1 0 0;
  ];

[~,~,V] = svd(A);
w_v = V(:,end);
w_v = w_v / w_v(end);
w = [ w_v(1) w_v(2) w_v(3); 
      w_v(2) w_v(4) w_v(5); 
      w_v(3) w_v(5) w_v(6) ];

P = Pproj(4:6,:);
M = P(1:3,1:3);
A = chol(inv(M' * w * M));

% Build Ha
Ha = zeros(4,4);
Ha(4, 4) = 1;
Ha(1:3,1:3) = inv(A);

% in case you don't get a positive definite matrix in the metric 
% reconstruction in the synthetic case please change the translation 
% vector of the first camera to this new one:  t1 = -R1*[42; 5; 10];

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
title('3D scene metric rectified');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Projective reconstruction (real data)

% read images
Irgb{1} = double(imread('Data/0000_s.png'))/255;
Irgb{2} = double(imread('Data/0001_s.png'))/255;

I{1} = sum(Irgb{1}, 3) / 3;
I{2} = sum(Irgb{2}, 3) / 3;

Ncam = length(I);

% ToDo: compute a projective reconstruction using the factorization method
% Compute keypoints and matches for all pairs of views
points = cell(Ncam,1);
descr = cell(Ncam,1);
% Descriptors
for i = 1:Ncam
  [points{i}, descr{i}] = sift(I{i}, 'Threshold', 0.01);
  points{i} = points{i}(1:2,:);
end

% Matches
matches_mtx = [];
for i = 1:Ncam-1
  matches = siftmatch(descr{i}, descr{i+1});
  x1 = points{i}(:, matches(1, :));
  x2 = points{i+1}(:, matches(2, :));
  matches_mtx = [matches_mtx; homog(x1); homog(x2)];
end

% Run the Factorization method to retrieve Pproj and Xproj
[Pproj, Xproj] = factorization_method(matches_mtx, 'SturmAndTriggs');

% ToDo: show the data points (image correspondences) and the projected
% points (of the reconstructed 3D points) in images 1 and 2. Reuse the code
% in section 'Check projected points' (synthetic experiment).

% Check projected points (estimated and data points)
for i=1:Ncam
  x_proj{i} = euclid(Pproj(3*i-2:3*i, :)*Xproj);
end

% Note: ideally this should be generalised to Ncam views
figure('Name', 'View 1: correspondences vs reconstructed 3D points');
imshow(Irgb{1},[]);
hold on;
% Plot matches
x1 = euclid(matches_mtx(1:3,:));  % Added 'euclid' to be formal, matches(1:2,:)

plot(x1(1,:), x1(2,:), 'r*');
plot(x_proj{1}(1,:),x_proj{1}(2,:),'bo');
axis equal
hold off;

figure('Name', 'View 2: correspondences vs reconstructed 3D points');
imshow(Irgb{2},[]);
hold on;

% Plot matches
x2 = euclid(matches_mtx(4:6,:));  % will do the same since the third row=1

plot(x2(1,:), x2(2,:), 'r*');
plot(x_proj{2}(1,:),x_proj{2}(2,:),'bo');
axis equal
hold off;


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
folder_out = 'output'; % output folder
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

