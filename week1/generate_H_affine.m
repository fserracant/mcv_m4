function H = generate_H_affine(theta, phi, lambda1, lambda2, t_x, t_y)
% ToDo: generate a matrix H which produces an affine transformation
% H will be of the shape: H = [a11 a12 tx
%                              a21 a22 ty
%                              0   0   1 ]
%
% if det(A) > 0 then, A has a SVD decomposition as follows:
% A = UDV^T = UV^T(VDV^T) = R(theta) R(-phi) D R(phi) where R(angle) stands
% for a rotation matrix 2x2 by the angle 'angle'. D is a diagonal matrix
% 2x2 like D = [ lambda1 0
%                0       lambda2]

% Compute rotation matrices (R)
R_theta = [cos(theta), -sin(theta); sin(theta), cos(theta)];
R_phi = [cos(phi), -sin(phi); sin(phi), cos(phi)];
R_phi_minus = [cos(phi), sin(phi); -sin(phi), cos(phi)];

% Compute diagonal matrix (D)
D = diag([lambda1, lambda2]);

% Compute A following the formula above:
A = R_theta * R_phi_minus * D * R_phi; % Double-check order

% Generate H by properly concatenating A, the translations (t_x, t_y) and
% the homogeneous coordinates [0, 0, 1]
H = [A(1, :), t_x; A(2, :), t_y; 0, 0, 1];
