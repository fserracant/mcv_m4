% Translation
H = [1 0 0; 0 1 0; 3 5 1]';
[r1, r2, sc, tr] = decompose_H_affine(H);
disp('Translation OK')

% Scale
H = [2 0 0; 0 2 0; 0 0 1]';
[r1, r2, sc, tr] = decompose_H_affine(H);
disp('Scale OK')

% Rotation
theta = pi/4;
H = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]';
[r1, r2, sc, tr] = decompose_H_affine(H);
disp('Rotation OK')

% Roration + Translation + Scale
theta = pi/4; 
H = [2*cos(theta) -sin(theta) 0; sin(theta) 2*cos(theta) 0; 3 5 1]';
[r1, r2, sc, tr] = decompose_H_affine(H);
disp('Roration + Translation + Scale OK')
