function [Pproj, Xproj] = factorization_method(xh, points_views, lambda_init)
% Implements the factorization method from [Sturm96] to do projective
% reconstruction (recover projective structure + motion).
%
% Inputs
%   - xh:             a matrix of 3*Ncam x Npoints corresponding to
%                     each view.
%                     (1-2, 2-3, 3-4,..., m-1-m).
%   - points_views:   location of the keypoints found for each view.
%                     It has size 2 x Npoints (i.e.: non homogeneous).
%   -
%
% Outputs
%   - Pproj:          3*Ncam x 4 matrix containing the camera matrices.
%   - Xproj:          4 x Npoints matrix of homogeneous coordinates of 3D points.
%
% References
%
% [Sturm96]: "A Factorization Based Algorithm for Multi-Image Projective Structure and
% Motion", Peter Sturm and Bill Triggs.
%
% [Har95]: "Multiple view geometry in Computer Vision", Richard Hartley and
% Andrew Zisserman.

%% Step 0: some initializations
Ncam = size(xh, 1)/3;  % a.k.a. 'Ncam'
ransac_thr = 2.0;  % test different values if needed
%% Step 1: normalise image coordinates (centroid at origin, dist. sqrt(2))
% Is this really needed if the Normalized 8-point algorithm already does
% this?
xh_norm = zeros(size(xh));
T = struct(Ncam,1);
for i = 1:3:Ncam
  [xh_norm(i:i+2,:), T{i}] = normalise2dpts(xh(i:i+2,:));
end

%% Step 2: compute fundamental matrix/ces with the Robust 8-point algorithm*
% * maybe to speed up this we'll need to use "only" the Normalised version.
F = struct(Ncam-1,1);
idx_inliers = struct(Ncam-1,1);

for i = 1:Ncam-1
%   x1 = points_views{i}(:, xh{i}(1, :));
%   x2 = points_views{i+1}(:, xh{i}(2, :));
  x1 = points_views{i}(:, xh_norm(3*i-2, :));
  x2 = points_views{i+1}(:, xh_norm(3*i-1, :));
  [F{i}, idx_inliers{i}] = ransac_fundamental_matrix(homog(x1), homog(x2), ransac_thr);
  % F = fundamental_matrix(homog(x1), homog(x2));  % Normalized 8-point algorithm
end

%% Step 3: determine scale factors (lambda's) via eq(3) of [Sturm96]
if strcmpi(lambda_init, 'ones')  % All lambda(i,j) set to 1
  lambda = ones(size(xh));
elseif strcmpi(lambda_init, 'SturmAndTriggs')  % only l1j is initialized to 1.
  % TODO: compute lambdas as in [Sturm96]
else
  error('Non-valid weight-initialization method. ''ones'' or ''sturm'' are allowed\n');
end

% TODO: embed following steps in a while() loop until convergence
coverged = false;

%% Step 4: build the rescaled measurement matrix W
% W has size 3m x n and rank<=4
W = lambda .* xh_norm;

%% Step 5: balance W by column and "triplet-of-rows"-wise normalisation
% [Sturm96] states a 3-step normalisation:
%   1. Rescale each column l so that sum(wrl.^2) = 1
%   2. Rescale each triplet of rows (3k-2,3k-1,3k) so sum(sum(wil.^2))=1
%   3. If the entries of W changed "significantly", repeat 1 and 2.

% Note: we follow the suggested approach in [Har95] (set rows-wise norms
% to 1 and then set columns similarly, in two passes). Notice that by
% doing this, we are not explicitly implementing point 3 above.

% Note 2: in [Sturm96], the authors say that balancing lambda alone would
% be enough as we are working with normalised image coordinates (instead of
% W as a whole). We may want to try this in the future and compare both results.

p = 2;  % Norm used, here L2.
% 1. Normalise rows
rows_norm = vecnorm(W, p, 2);
W_b = W ./ repmat(rows_norm, [1, size(W,2)]);
% 2. Normalise columns
cols_norm = vecnorm(W_b, p, 1);
W_b = W_b ./ repmat(cols_norm, [size(W_b,1), 1]);

%% Step 6: compute the SVD of balanced W
[U,D,V] = svd(W_b);
% Force W_b to have rank 4 by setting eigenvalues


%% Step 7: recover projective motion + shape from SVD decomposition
% Identify motion and shape from SVD

%% Step 8: adapt projective motion (Pi) accounting for Ti => denormalise
% Pi_hat' = inv(Ti) * Pi_hat
Pproj = [];
for i = 1:Ncam
  aux = T{i} \ Pproj_hat;
  Pproj(3*i-2:3*i, :) = aux;
end
% Compute projected 3D points.
Xproj = Pproj' * xh;  % Pproj is 3*Ncam x 4, xh 3*Ncam x Npoints=> Xproj is 4 x Npoints

end