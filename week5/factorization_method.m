function [Pproj, Xproj] = factorization_method(xh, lambda_init)
% Implements the factorization method from [Sturm96] to do projective
% reconstruction (recover projective structure + motion).
%
% Inputs
%   - xh:             a matrix of 3*Ncam x Npoints matches corresponding to
%                     each view.
%                     (1-2, 2-3, 3-4,..., m-1-m).
%   - points_views:   location of the keypoints found for each view.
%                     It has size 2 x Npoints (i.e.: non homogeneous).
%   - lambda_init:    how to compute the depth weights 'lambda'. 'ones'
%                     sets all lambdas to 1, 'SturmAndTriggs' uses the method in [Sturm96].
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

% Initializations
Ncam = size(xh, 1) / 3;  % a.k.a. 'Ncam'
Npoints = size(xh, 2);
ransac_thr = 2.0;  % test different values if needed

% Normalization
xh_norm = zeros(size(xh));
T = cell(Ncam,1);
for i = 1:Ncam
  [xh_norm(3*i-2:3*i,:), T{i}] = normalise2dpts(xh(3*i-2:3*i,:));
end

if strcmpi(lambda_init, 'SturmAndTriggs')
  %we apply the same constraint as [Sturm96]: the images are taken
  %   pairwise in sequence, F12, F23, F34,..., Fm-1m
  
  %TODO: first naive implementation, naïve loops for all views and all points
  lambda = ones(Ncam, Npoints);  % Lambda for the first view is 1
  iref = 1;
  
  for i = 2:Ncam
    xref = xh_norm(3*iref-2:3*iref, :);
    xi = xh_norm(3*i-2:3*i, :);
    [Fri, ~, eri, eir] = ransac_fundamental_matrix(xi, xref, ransac_thr);
    
    % Determine scale factors (lambda's) via eq(3) of [Sturm96]
    for p = 1:Npoints
      xrp = xref(:,p);
      xip = xi(:,p);
      cross_eri_xip = cross(eri, xip);
      lambda(i,p) = lambda(iref,p) * ...
        abs((xrp' * Fri * cross_eri_xip) /...
        (cross_eri_xip' * cross_eri_xip));
    end
    iref = i;  % Update reference image.
  end
  
elseif strcmpi(lambda_init, 'ones')  % All lambda(i,j) set to 1
  lambda = ones(Ncam, Npoints);
else
  error('Non-valid weight-initialization method. ''ones'' or ''sturm'' are allowed\n');
end

dist = inf;
iteration_number = 0;
convergence_threshold = 1e-6;
converged = false;

while ~converged
  iteration_number = iteration_number+1;
  sprintf('iteration: %d', iteration_number)
  
  p_norm = 1;  % Norm used
  % Normalise columns and rows twice
  lambda = rescaleUntilConvergence(lambda, p_norm, 1e-9);
  
  
  W_b = zeros(3*Ncam, Npoints);  % Define measurement matrix W

  for i = 1:Ncam
    W_b(3*i-2:3*i,:) = lambda(i,:) .* xh_norm(3*i-2:3*i,:);
  end
  
  % Compute the SVD of balanced W
  [U,D,V] = svd(W_b);
  D4 = diag(D);
  % Sigma = D4(1:4);  % not used as Sigma = Sigma_p * Sigma_pp
  Sigma_p = diag(sqrt(D4(1:4)));  % Force rank 4
  Sigma_pp = Sigma_p;
  % W = U' * Sigma_p * Sigma_pp * V' = U_hat * V_hat (eq(4))
  % where U' contains the first 4 columns of U, V' the first 4 rows of V
  %   W_b = U(:,1:4) * Sigma_p * Sigma_pp * V(:,1:4)';
  
  % Recover projective motion + shape from SVD decomposition
  P_hat = U(:,1:4) * Sigma_p;
  X_hat = Sigma_pp * V(:,1:4)';
  
  % Compute reprojection error
  dist_old = dist;
  % Euclidean distance, make sure that we are using Euclidean coordinates!
  dist = sum(sqrt(sum((xh_norm - (P_hat * X_hat)) .^2)));
  converged = (abs(dist - dist_old) / dist) < convergence_threshold;
  if ~converged
    % Update lambda
    xp = P_hat * X_hat;
    for i = 1:Ncam
      lambda(i,:) = xp(3*i,:);
    end
  end
  
end

Pproj = P_hat;
Xproj = X_hat;

% Adapt projective motion (Pi) accounting for Ti => denormalise
for i = 1:Ncam
  aux = T{i} \ P_hat(3*i-2:3*i, :);
  Pproj(3*i-2:3*i, :) = aux;
end

end

function [lambda_r] = rescaleLambda(lambda, p_norm)
  % Re-scale each column
  lambda_r = lambda./vecnorm(lambda,p_norm,1);
  
  % Re-scale triplets of rows
  for k = 1:3:size(lambda,1)-2
    triplet = lambda_r(k:k+2, :);
    lambda_r(k:k+2,:) = triplet ./ vecnorm(triplet(:), p_norm, 1);
  end
end

function [lambda_r] = rescaleUntilConvergence(lambda, p_norm, th)
  tmp = lambda;
  lambda_r = lambda;
  convergence = false;
  while ~convergence
    tmp = rescaleLambda(lambda_r, p_norm);
    sum(abs(lambda_r(:) - tmp(:)))
    if sum(abs(lambda_r(:) - tmp(:))) < th
      convergence = true;
    end
    lambda_r = tmp;
  end
end
