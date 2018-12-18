function [I_metrect] = compute_metricRectification(I_affrect)
% Two step metric rectification
% Assumption: the input image has already been affinely rectified.

% Define two set of non-parallel orthogonal lines (l1, m1) (l2, m2) 
% l = (l1, l2, l3), m = (m1, m2, m3), l and m are orthogonal
lines = getLinesCoordinates(I_affrect, 1);
l1 = [lines(1,1), lines(1,2)];
m1 = [lines(2,1), lines(2,2)];
l2 = [lines(3,1), lines(3,2)];
m2 = [lines(4,1), lines(4,2)];

% Step 1: define system of equations determined by I^T Mm
M = [l1(1) * m1(1), (l1(1) * m1(2) + l1(2) * m1(1));
  l2(1) * m2(1), (l2(1) * m2(2) + l2(2) * m2(1))];
b = [-l1(2) * m1(2); -l2(2) * m2(2)];

% Step 2: solve the above system of equations
x = linsolve(M, b);

% Define matrix S (2x2 symmetric)
S = eye(2);
S(1,1) = x(1);
S(1,2) = x(2);
S(2,1) = x(2);
% S(2,2) = 1

% Use the Cholesky decomposition of S to compute an upper triangular matrix
% K such that S = K K'
% PROBLEM: function 'chol' needs positive definite matrices and 'cholcov'
% outputs an empty K

% Found implementation using a decomposition via SVD.
% Discussion:
% from: https://stackoverflow.com/questions/12449695/metric-rectification-using-the-dual-degenerate-conic-in-matlab
% "Theoretically, since the distorted conic is determined by C*'=HC*H' (C*
% is the dual degenerate conic, ' is transpose, H is my homography)..."
% Actual repo with implementation:
% https://github.com/syedhope/Image_Rectification/blob/master/metric_rect1.m

[U, D, V] = svd(S);
sqrtD = sqrt(D);
U_t = transpose(U);
A = U_t * sqrtD;
A = A * V;

% Define homography
H2 = eye(3);
H2(1,1) = A(1,1);
H2(1,2) = A(1,2);
H2(2,1) = A(2,1);
H2(2,2) = A(2,2);
if H2(1,1) < 0
  H2(1,1) = -H2(1,1);
  
elseif H2(2,2) < 0
  H2(2,2) = -H2(2,2);
end

% ---- COMMENTS ----
% This implementation uses Matlab built-in deprecated 'maketform'
% which transposes the homography w.r.t. to the affine model
% Still wondering why, double check!
% TODO: test 'apply_H' with both H and H2... compare w. reference

H = H2';
I_metrect = apply_H(I_affrect, H2);