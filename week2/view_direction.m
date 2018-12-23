function [ v ] = view_direction( P, x )
%VIEW_DIRECTION Vector pointing to the viewing direction of a pixel.
%   We solve x = P v with v(3) = 0;

v = P(:,1:3)\[x(1); x(2); 1];

end

