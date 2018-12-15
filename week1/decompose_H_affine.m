function H_dec = decompose_H_affine(A, t_x, t_y)
% ToDo: decompose the affinity in four transformations: two
% rotations, a scale, and a translation
% (probably) we should decompose A into its SVD decomposition of 3 matrices
% (2 rotations + 1 scale) and then add the translation at the end.
% H_prime = Transl * SVD_A = Transl * Rotation1 * Rotation2 * scaling

% Define translation matrix
Transl = [1, 0, t_x; 0, 1, t_y; 0, 0, 1];

% Compute SVD decomposition of A (1 diagonal mtx => scale, 2 square
% matrices => rotations)
[Rot1, Scaling, Rot2] = svd(A);

% Obtain the decomposition of H into 4 matrices as follows:
H_dec = Transl * Rot1 * Rot2 * Scaling;