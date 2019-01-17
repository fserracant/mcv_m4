function [ o ] = optical_center( P )
%OPTICAL_CENTER Computes the optical center of a camera.
%   Detailed explanation goes here

[U,S,V] = svd(P);
o = V(:,end);
o = o(1:3) / o(4);

end

