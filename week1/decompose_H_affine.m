function [rotation1, rotation2, scale, translation] = decompose_H_affine(H)
% DECOMPOSE_H_AFFINE Decomposes an affine transformation in four transformations two rotations, a scale and a translation
%   [r1, r2, sc, tr] = DECOMPOSE_H_AFFINE(H) decomposes an affine
%   transformation in two rotations, an scaling and a translation

  tolerance = 1e-6;

  % Jon: Any normalization ?
  A = H(1:2, 1:2);
  t = H(1:2, 3);
  v = H(3, 1:2);
  assert(v(1)==0 & v(2)==0 & H(3,3)==1, ...
    'Error: Matrix does not correspond to an affine transformation')

  [U, D, Vt] = svd(A);
  R_phi = Vt;
  R_mphi = Vt';
  R_theta = U*Vt';

  A_prime = R_theta * R_mphi * D * R_phi;
  H_prime = [A_prime(1,:), t(1); A_prime(2,:), t(2); 0 0 1];
  assert(sum(H(:) - H_prime(:) < tolerance)==9, ...
    'Bug1: cannot recreate input from SVD result')

  rotation1 = [R_phi(1,:) 0; R_phi(2,:) 0; 0 0 1];
  rotation2 = [R_theta(1,:) 0; R_theta(2,:) 0; 0 0 1];
  scale = [D(1,1) 0 0; 0 D(2,2) 0; 0 0 1];
  translation = [1 0 t(1); 0 1 t(2); 0 0 1];
  H_composed = translation * rotation2 * rotation1' * scale * rotation1;
  assert(sum(H(:) - H_composed(:) < tolerance)==9, ...
    'Bug2: cannot recreate input from 4 transformations')
end
