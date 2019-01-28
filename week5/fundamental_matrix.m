% FUNDMATRIX - computes fundamental matrix from 8 or more points
%
% Function computes the fundamental matrix from 8 or more matching points in
% a stereo pair of images.  The normalised 8 point algorithm given by
% Hartley and Zisserman p265 is used.
%
% Usage:   F = fundmatrix(x1, x2)
%         
% Arguments:
%          x1, x2 - Two sets of corresponding 3xN set of homogeneous
%          points.
%         
% Returns:
%          F      - The 3x3 fundamental matrix such that x2'*F*x1 = 0.
%

% Based on Peter Kovesi's code

function F = fundamental_matrix(x1, x2)

    npts = size(x1,2); 
    
    % Normalise each set of points so that the origin 
    % is at centroid and mean distance from origin is sqrt(2). 
    % normalise2dpts also ensures the scale parameter is 1.
    [x1, T1] = normalise2dpts(x1);
    [x2, T2] = normalise2dpts(x2);
    
    % Build the constraint matrix
    A = [x2(1,:)'.*x1(1,:)'   x2(1,:)'.*x1(2,:)'  x2(1,:)' ...
         x2(2,:)'.*x1(1,:)'   x2(2,:)'.*x1(2,:)'  x2(2,:)' ...
         x1(1,:)'             x1(2,:)'            ones(npts,1) ];       

	[U,D,V] = svd(A,0); % use the economy decomposition

    % Extract fundamental matrix from the column of V corresponding to
    % smallest singular value.
    F = reshape(V(:,9), 3, 3)';
    
    % Enforce constraint that fundamental matrix has rank 2 by performing
    % a svd and then reconstructing with the two largest singular values.
    [U,D,V] = svd(F);
    F = U * diag([D(1,1) D(2,2) 0]) * V';
    
    % Denormalise
    F = T2'*F*T1;
    