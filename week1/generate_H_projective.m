function Hp = generate_H_projective(theta, phi, lambda1, lambda2, t_x, t_y, v)
% Generates a 2D projective transformation. It does so creating an affine
% transformation and then adding the vector v = (v1,v2).

Ha = generate_H_affine(theta, phi, lambda1, lambda2, t_x, t_y);

Hp = [Ha(1,:); Ha(2,:); v, 1];