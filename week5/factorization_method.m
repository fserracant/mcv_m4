function [Pproj, Xproj] = factorization_method(matches_views, points_views, weight_init)
% Implements the factorization method from [Sturm96] to do projective
% reconstruction (recover projective structure + motion).
%
% Inputs
%   - matches_views:  a struct of matching matrices for all pairs of views
%                     (1-2, 2-3, 3-4,..., m-1-m).
%   - points_views:   location of the keypoints found for each view.
%
% Outputs
%   - Pproj:          3*Ncam x 4 matrix containing the camera matrices.
%   - Xproj:          4 x Npoints matrix of homogeneous coordinates of 3D points.
%
% [Sturm96]: "A Factorization Based Algorithm for Multi-Image Projective Structure and
% Motion", Peter Strum and Bill Triggs.
%
% Implementation details: 
% the goal would be to implement it for an arbitrary number of views m. 
% However, for the sake of simplicity, we MAY need to start by assuming m = 2.

%% Step 0: some initializations
nViews = length(matches_views);  % a.k.a. 'Ncam'
ransac_thr = 2.0;  % test different values if needed
%% Step 1: normalise image coordinates (centroid at origin, dist. sqrt(2))
% Is this really needed if the Normalized 8-point algorithm already does
% this?
% for i = 1:nViews
%   
% end

%% Step 2: compute fundamental matrix/ces with the Robust 8-point algorithm*
% * maybe to speed up this we'll need to use "only" the Normalised version.
F = struct(nViews-1,1);
idx_inliers = struct(nViews-1,1);

for i = 1:nViews-1
  x1 = points_views{i}(:, matches_views{i}(1, :));
  x2 = points_views{i+1}(:, matches_views{i}(2, :));
  [F{i}, idx_inliers{i}] = ransac_fundamental_matrix(x1, x2, ransac_thr);
end

%% Step 3: determine scale factors (lambda's) via eq(3) of [Sturm96]
if strcmpi(weight_init, 'ones')  % All lambda(i,j) set to 1
%weights = ones(nViews, 
elseif strcmpi(weight_init, 'sturm')  % only l1j is initialized to 1.
  
else
  error('Non-valid weight-initialization method. ''ones'' or ''sturm'' are allowed\n');
end

%% Step 4: build the rescaled measurement matrix W
% W has size 3m x n and rank<=4

%% Step 5: balance W by column and "triplet-of-rows"-wise normalisation

%% Step 6: compute the SVD of balanced W

%% Step 7: recover projective motion + shape from SVD decomposition
% Identify motion and shape from SVD

%% Step 8: adapt projective motion (Pi) accounting for Ti => denormalise
% Pi_hat' = inv(Ti) * Pi_hat