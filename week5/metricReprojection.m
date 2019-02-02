function [Ha] = metricReprojection(v1,v2,v3,Pproj,Hp)
% METRICEPROJECTION computes the homography Ha that transforms a
% affine image into a metric one.
% It works by mapping by finding the absolute conic in one of the images
%
% Inputs
%   - v1:           first vanishing point in the first view.
%   - v2:           second vanishing point in the first view.
%   - v3:           third vanishing point in the first view.
%   - Pproj:        contains the camera matrices after projective
%                   reconstruction.
%   - Hp:            homography that updates the projective reconstruction
%                   into an affine one.
%
% Outputs
%   - Ha:           homography that updates the projective reconstruction
%                   into an affine one.

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

end