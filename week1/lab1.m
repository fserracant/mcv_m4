%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 1: Image rectification


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Applying image transformations

% ToDo: create the function  "apply_H" that gets as input a homography and
% an image and returns the image transformed by the homography.
% The size of the transformed image has to be automatically set so as to 
% contain the whole transformed image.
% At some point you will need to interpolate the image values at some points,
% you may use the Matlab function "interp2" for that.


%% 1.1. Similarities
close all;
clear;
I=imread('Data/0005_s.png'); % we have to be in the proper folder

% ToDo: generate a matrix H which produces a similarity transformation
s = 0.5;
theta = pi/4;
t_x = 0; t_y = 0;
H = [s*cos(theta) s*-sin(theta) t_x; s*sin(theta) s*cos(theta) t_y; 0 0 1];

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));


%% 1.2. Affinities

% ToDo: generate a matrix H which produces an affine transformation

lambda1 = 1;
lambda2 = 0.5;
theta =pi/4;
phi = pi/4;
t_x = 0; t_y = 0;

H = generate_H_affine(theta, phi, lambda1, lambda2, t_x, t_y);



I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));

% ToDo: decompose the affinity in four transformations: two
% rotations, a scale, and a translation
[rotation1, rotation2, scale, translation] = decompose_H_affine(H);



% ToDo: verify that the product of the four previous transformations
% produces the same matrix H as above

H_comp = translation * rotation2 * rotation1' * scale * rotation1;
assert(verify_product_H_affine(H, H_comp), 'H matrices not equal');

% ToDo: verify that the proper sequence of the four previous
% transformations over the image I produces the same image I2 as before

I_comp = uint8(apply_H(I,rotation1));
I_comp = uint8(apply_H(I_comp,scale));
I_comp = uint8(apply_H(I_comp,rotation2 * rotation1'));
I_comp = apply_H(I_comp,translation);


figure; imshow(uint8(I2));figure; imshow(uint8(I_comp));

%assert(verify_images(I2, I_comp), 'Result images not equal');



%% 1.3 Projective transformations (homographies)

% ToDo: generate a matrix H which produces a projective transformation
lambda1 = 1;
lambda2 = 0.5;
theta =pi/4;
phi = pi/4;
t_x = 0; t_y = 0;
v = [0.000075, -0.0005];
Hp = generate_H_projective(theta, phi, lambda1, lambda2, t_x, t_y, v)


I2 = apply_H(I, Hp);
figure; imshow(I); figure; imshow(uint8(I2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Affine Rectification

% choose the image points
I=imread('Data/0000_s.png');
A = load('Data/0000_s_info_lines.txt');

% crop the image to get only right facade
col_crop = 270;
[cols,rows] = size(I);
I = I(:,col_crop:cols,:);

% indices of lines
i = 424;
p1 = [A(i,1)-col_crop A(i,2) 1]';
p2 = [A(i,3)-col_crop A(i,4) 1]';
i = 240;
p3 = [A(i,1)-col_crop A(i,2) 1]';
p4 = [A(i,3)-col_crop A(i,4) 1]';
i = 712;
p5 = [A(i,1)-col_crop A(i,2) 1]';
p6 = [A(i,3)-col_crop A(i,4) 1]';
i = 565;
p7 = [A(i,1)-col_crop A(i,2) 1]';
p8 = [A(i,3)-col_crop A(i,4) 1]';

% ToDo: compute the lines l1, l2, l3, l4, that pass through the different pairs of points
l1 = create_line(p1, p2); l1 = l1 / l1(3);
l2 = create_line(p3, p4); l2 = l2 / l2(3);
l3 = create_line(p5, p6); l3 = l3 / l3(3);
l4 = create_line(p7, p8); l4 = l4 / l4(3);

% show the chosen lines in the image
figure;imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');

% ToDo: compute the homography that affinely rectifies the image
% compute vanishing points where lines cross at ininity
v1 = cross(l1, l2)
v2 = cross(l3, l4)

% compute line that passes through vanishing points: line at infinite 
l_inf = cross(v1, v2); 
l_inf = l_inf / l_inf(3)

% define H for affine rectification and apply to image
H = [1 0 0; 0 1 0; l_inf(1) l_inf(2) 1];
I2 = apply_H(I, H);
figure; imshow(uint8(I2));

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4
% by using l' = H^(-T) * l
lr1 = H' \ l1; lr1 = lr1 / lr1(3);
lr2 = H' \ l2; lr2 = lr2 / lr2(3);
lr3 = H' \ l3; lr3 = lr3 / lr3(3);
lr4 = H' \ l4; lr4 = lr4 / lr4(3);

% show the transformed lines in the transformed image
figure;imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'y');

% ToDo: to evaluate the results, compute the angle between the different pair 
% of lines before and after the image transformation

% calculate slopes of original lines
s1 = l1(1) / l1(2); 
s2 = l2(1) / l2(2);
s3 = l3(1) / l3(2);
s4 = l4(1) / l4(2);
% calculate slope of rectified lines
sr1 = lr1(1) / lr1(2); 
sr2 = lr2(1) / lr2(2);
sr3 = lr3(1) / lr3(2);
sr4 = lr4(1) / lr4(2);

% angles between "parallel" original lines
a12 = rad2deg( atan(s1) - atan(s2) ) % a12 = 0.0992
a34 = rad2deg( atan(s3) - atan(s4) ) % a34 = -1.3435
% angles of really parallel rectified lines
ar12 = rad2deg( atan(sr1) - atan(sr2) ) % ar12 = -7.9514e-16 (almost 0)
ar34 = rad2deg( atan(sr3) - atan(sr4) ) % ar34 = 0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric Rectification

%% 3.1 Metric rectification after the affine rectification (stratified solution)

% ToDo: Metric rectification (after the affine rectification) using two non-parallel orthogonal line pairs
%       As evaluation method you can display the images (before and after
%       the metric rectification) with the chosen lines printed on it.
%       Compute also the angles between the pair of lines before and after
%       rectification.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Affine and Metric Rectification of the left facade of image 0001

% ToDo: Write the code that rectifies the left facade of image 0001 with
%       the stratified method (affine + metric). 
%       Crop the initial image so that only the left facade is visible.
%       Show the (properly) transformed lines that use in every step.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Metric Rectification in a single step
% Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)



