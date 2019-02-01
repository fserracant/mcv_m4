function [Hp] = affineReprojection(v1,v1p,v2,v2p,v3,v3p,Pproj,w,h)
% AFFINEREPROJECTION computes the homography Hp that transforms a
% projective constructed image into an affine one.
% It works by mapping by triangulation the image of 3 pairs of vanishing
% points found in two views of the same scene.
%
% Inputs
%   - v1:           first vanishing point in the first view.
%   - v1p:          first vanishing point in the second view.
%   - v2:           second vanishing point in the first view.
%   - v2p:          second vanishing point in the second view.
%   - v3:           third vanishing point in the first view.
%   - v3p:          third vanishing point in the second view.
%   - Pproj:        contains the camera matrices after projective
%                   reconstruction.
%   - w:            width of the images.
%   - h:            height of the images.
%
% Outputs
%   - Hp:           homography that updates the projective reconstruction
%                   into an affine one.

if size(v1,1) == 3  || size(v1p,1) == 3 % assuming same size for one image vp's.
  v1 = euclid(v1);
  v1p = euclid(v1p);
  v2 = euclid(v2);
  v2p = euclid(v2p);
  v3 = euclid(v3);
  v3p = euclid(v3p);
end

Xv1 = triangulate(v1, v1p, Pproj(1:3,:), Pproj(4:6,:), [h,w]);
Xv2 = triangulate(v2, v2p, Pproj(1:3,:), Pproj(4:6,:), [h,w]);
Xv3 = triangulate(v3, v3p, Pproj(1:3,:), Pproj(4:6,:), [h,w]);

% Solve Ap = 0 ==> p is the null vector of A ==> SVD and pick last column
% (eigenvectors associated to lowest eigenvalue)
A = [Xv1'; Xv2'; Xv3'];
[~,~,V] = svd(A);
p = V(:,end);
p = p / p(end);

Hp = eye(4,4);
Hp(4,:) = p';

end