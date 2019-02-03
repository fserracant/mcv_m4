function [vp] = vanishing_point(xo1, xf1, xo2, xf2)
  %VANISHING_POINT computes the vanishing point 'vp' formed by the line
  % that joins points xo1 and xf1 and the line that joins points x02 and xf2.
  %
  % Inputs
  %   - xo1:      first point of the first line.
  %   - xf1:      second point of the first line.
  %   - xo2:      first point of the second line.
  %   - xf2:      second point of the second line.
  %
  % Outputs
  %   - vp:       vanishing point formed by the two lines defined by x01-xf2.
  %

  tol = 1e-6;

  % Check input vectors dimensions
  assert(size(xo1,1) == 3 && size(xf1,2) == 1, ...
    'first points are not 3-dim column vectors')

  assert(size(xo2,1) == 3 && size(xf2,2) == 1, ...
    'second points are not 3-dim column vectors')

  % Compute line that joins first set of points xo1,xf1
  l1 = cross(xo1, xf1);
  assert(dot(l1,xo1) < tol & dot(l1,xf1) < tol, ...
    'Bug: line ''l1'' is not perpendicular to given points')
  l1 = l1 / l1(3);

  % Compute line that joins the second set of points xo2,xf2
  l2 = cross(xo2, xf2);
  assert(dot(l2,xo2) < tol & dot(l2,xf2) < tol, ...
    'Bug: line ''l1'' is not perpendicular to given points')
  l2 = l2 / l2(3);

  % Compute vanishing point as the intersection of both lines
  vp = cross(l1, l2);
end