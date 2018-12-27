function [d2] = compute_symmetricGeometricError(x1, x2, H)
% Compute the symmetric geometric error (backprojection error)
%
% Inputs:
%   - x1:     matches from first image.
%   - x2:     matches from second image.
%   - H:      homography that relates x1 ==> x2 (or x1 <== x2 via inv(H))
%
% Ouput:
%   - d2:     symmetric geometric error (d2 = d(Hx1,x2) + d(x1, inv(H)x2))

% Forward-projected d(H*x1,x2)
forward = H * x1;
forward_n = forward(1:2,:) ./ forward(3,:);

d2_fwd = sqrt(sum((forward_n - x2(1:2,:)) .^2));

% Backward-projected d(x1,inv(H)*x2)
backward = H \ x2;
backward_n = backward(1:2,:) ./ backward(3,:);
d2_bwd = sqrt(sum((x1(1:2,:) - backward_n) .^2));

% Combine both terms
d2 = d2_fwd + d2_bwd;