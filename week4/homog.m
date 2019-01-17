function [ xh ] = homog( xe )
%HOMOG Homogeneous coordinates from euclidean

xh = [xe; ones(1, size(xe, 2))];

end

