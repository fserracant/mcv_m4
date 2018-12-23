function [ xe ] = euclid( xh )
%EUCLID Euclidean coordinates from homogeneous
D = size(xh, 1) - 1;
xe = xh(1:D,:) ./ repmat(xh(end,:), D, 1);

end

