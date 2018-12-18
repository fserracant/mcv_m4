function [line] = create_line(point1, point2)
  % Computes a line that passes through two points in homogeneous coordinates

  assert(size(point1,1) == 3 && size(point1,2) == 1, ...
      'point1 is not a 3-dim column vector')
  assert(size(point2,1) == 3 && size(point2,2) == 1, ...
      'point2 is not a 3-dim column vector')

  line = cross(point1,point2);

  % allow some tolerance because some numbers are not exactly zero but very low
  % instead
  tol = 1e-6;
  assert(dot(line,point1) < tol & dot(line,point2) < tol, ...
   'Bug: computed line is not perpendicular to given points')
end
