function [x_tilde, T] = normalise_matches(x, n_corresp)
% Normalises the matches in x and returns the transformation matrix used.
%
% Input:
%   - x         vector of matches in one of the two views of the scene.

% Outputs:
%   - x_tilde:   normalised vector of matches.
%   - T:        transformation that maps x => x_norm (scale + translation).

% Compute centroid
x_centroid = mean(x(1:2,:), 2);  % for xi, yi, zi=1, remains fixed

% Subtract the centroid (done as a transformation)
trans_mtx = [1 0 -x_centroid(1); 0 1 -x_centroid(2); 0 0 1];
x_new = trans_mtx * x;

% 1.2. Scale x_new so its average distance is sqrt(2)
% sqrt(2) = 1/n sum(1..n) norm(xi_new)^2
% Compute sqrt distance
sqrt_dist = sum(sqrt(sum(x_new(1:2,:).^2)));
sqrt_dist = sqrt_dist / n_corresp;
% Compute scale
scale = sqrt(2) / sqrt_dist;

% Define (symmetric) scaling transformation
scale_mtx = [scale, 0, 0; 0 scale 0; 0 0 1];
% Apply scale to x_new
x_tilde = scale_mtx * x_new;

% Define 'global' transformation as S * T and its inverse (for denormalization)
T = scale_mtx * trans_mtx;
