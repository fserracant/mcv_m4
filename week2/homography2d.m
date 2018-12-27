function H = homography2d(x1, x2)
% Compute 2D homography between a set of correspondences x1 <==> x2

% x1 and x2 are sets of 4 points randomly sampled out of all the
% correspondences with homogeneous coordinates: (xi, yi, 1)
% Algorithm: Normalized DLT

% Check that we have enough correspondences
n_corresp = size(x1, 2);
assert(n_corresp >= 4, 'Not enough correspondences given (n>=4)');

%% Step 1: Normalization of x1
[x1_tilde, T_x1] = normalise_matches(x1, n_corresp);

%% Step 2: normalization of x2
[x2_tilde, T_x2] = normalise_matches(x2, n_corresp);

%% Step 3: apply the DLT algorithm
% Define the system of equations Ah = 0
% Compute and assemble A; size 2n x 9 (n = 4 in this case)
%A = zeros(2 * n_corresp, 9);
A = [];
for i = 1:n_corresp
  A = [A;zeros(1,3) (-1) * x1_tilde(:,i)' x2_tilde(2,i) * x1_tilde(:,i)';
    x1_tilde(:,i)' zeros(1,3) -x2_tilde(1,i) * x1_tilde(:,i)']; 
end

% Pick h as the last column of V^T in the SVD decomposition of A
[~, ~, Vt] = svd(A);

h = Vt(:, end);
% Reshape h into a 3x3 matrix
H_tilde = reshape(h, 3, 3);

%% Step 4: undo the translation and the scaling
% Denormalise => H = inv(T_x2) * H_tilde * T_x1
H = T_x2 \ H_tilde * T_x1;

% Normalise H (up to a factor)
% E.g.: divide everything by H(3,3)
H = H / H(end,end);
