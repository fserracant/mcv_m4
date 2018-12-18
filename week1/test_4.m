close all;
clear;

% choose the image points
I=imread('Data/0001_s.png');
A = load('Data/0001_s_info_lines.txt');

% crop the image to get only left facade
col_crop = 550;
I = I(:,1:col_crop,:);

% indices of lines
i = 159;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 614;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 541;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 645;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';
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
a12 = rad2deg( atan(s1) - atan(s2) ) % a12 = -4.4208
a34 = rad2deg( atan(s3) - atan(s4) ) % a34 = 0.1249
% angles of really parallel rectified lines
ar12 = rad2deg( atan(sr1) - atan(sr2) ) % ar12 = 0
ar34 = rad2deg( atan(sr3) - atan(sr4) ) % ar34 = 0


% Orthogonal pair of lines 
l1 = lr1;
m1 = lr3;
%  Orthogonal pair of lines 
l2 = lr2;
m2 = lr4;

% A = [l1(1)*m1(1),   l1(1)*m1(2)+l1(2)*m1(1),    l1(2)*m1(2);
%      l2(1)*m2(1),   l2(1)*m2(2)+l2(2)*m2(1),    l2(2)*m2(2)];
 
 A = [l2(1)*m2(1),   l2(1)*m2(2)+l2(2)*m2(1),    l2(2)*m2(2)];

s = null(A);


S = [s(1,2),  s(2,2); 
    s(2,2),   s(3,2)];

K = chol(S); 

H = eye(3);
K = inv(K); 
H(1:2,1:2) = K;
H = H';

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4
% by using l' = H^(-T) * l

lr1 = H' \ l1; lr1 = lr1 / lr1(3);
lr2 = H' \ l2; lr2 = lr2 / lr2(3);
m1r = H' \ m1; m1r = m1r / m1r(3);
m2r = H' \ m2; m2r = m2r / m2r(3);



% show the transformed lines in the transformed image
I3 = apply_H(uint8(I2), H);
figure; imshow(uint8(I3));
hold on;
t=1:0.1:10000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(m1r(1)*t + m1r(3)) / m1r(2), 'y');
plot(t, -(m2r(1)*t + m2r(3)) / m2r(2), 'y');
