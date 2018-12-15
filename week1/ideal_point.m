function [point] = ideal_point(line)
  % Computes the ideal point or point at inf of a line

  assert(size(line,1) == 3 & size(line,2) == 1, ...
    'line is not a 3-dim column vector')
  assert(line(1) ~= 0 | line(2) ~=0, ...
    'ideal point from this line would not be legal')

  point = [-line(2), line(1), 0]';

  assert(dot(point,line)==0 | dot(point, [0,0,1]')==0, ...
    'Bug: computed point is not perpendicular to given line or to ideal line')
end
