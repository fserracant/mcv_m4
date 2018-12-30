function [ d ] = geometric_dist(a,b)
% Computes the geometric distance between points in homogeneous coord
% a is a matrix 3xN of N points. The result is the same as computing the 
% euclidean distances in the inhomogenous coordinates.

assert(size(a,1) == 3)
assert(size(b,1) == 3)

d = (a(1,:) ./ a(3,:) - b(1,:) ./ b(3,:)).^2 + ...
    (a(2,:) ./ a(3,:) - b(2,:) ./ b(3,:)).^2;
