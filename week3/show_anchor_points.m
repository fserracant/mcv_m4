function [] = show_anchor_points(points)
  %SHOW_ANCHOR_POINTS Show descriptor's anchor points into ref image
  % It needs a figure open
  idx = 1:size(points, 2);
  label = cellstr(num2str(idx'));
  scatter(points(1,:), points(2,:), 200, 'filled', 'r')
  %text(points(1,:)+1, points(2,:), label);
end
