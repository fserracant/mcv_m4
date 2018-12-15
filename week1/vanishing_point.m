function [point] = vanishing_point(line1, line2)
  % Computes the vanishing point of two paralel lines in homogeneous coordinates

  assert(size(line1,1) == 3 && size(line1,2) == 1, ...
    'line1 is not a 3-dim column vector')
  assert(size(line2,1) == 3 && size(line2,2) == 1, ...
    'line2 is not a 3-dim column vector')

  point = cross(line1,line2);

  assert(dot(point,line1)==0 & dot(point,line2)==0, ...
    'Bug: computed point is not perpendicular to given lines')
end
