function [points] = getLinesCoordinates(I, method)
% Stratified metric rectification method:
% Draw two pairs of non-parallel orthogonal lines on the image to rectify I
% It is very important than lp1 is perpendicular to lp2 and lp3
% perpendicular to lp4 but the pairs should NOT be parallel.

% Metric rectification in one step (including affine + metric)
% We need 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)

% VERY IMPORTANT: draw the pairs in order, that is:
% First and second lines perpendicular to one another. Third and fourth
% also perpendicular to each other.

[height, width, ~] = size(I);
center_x = round(width/2);
center_y = round(height/2);
line_def_length = 50;
figure('Name', 'Select lines (''Enter'' to finish)');
imshow(I,[]);

if method == 1
  % Default, two-step method.
  n_points = 4;
  points = zeros(n_points, 4);
elseif method == 2
  % One-step method.  
  n_points = 10;
  points = zeros(n_points, 4);
  
else
  error('Not valid method, we do not know how many orthogonal lins to seek\n');
end

for i = 1:n_points
  % Let the user draw parallel lines (more straightforward than manually
  % defining them!)
  lineROI = images.roi.Line(gca, 'Position', [center_x, center_y; ...
    center_x, center_y + line_def_length], 'StripeColor','r');
  draw(lineROI);  % Draw interactively
% lineROI = drawline('Position',[center_x, center_y; center_x, center_y + line_def_length],...
%   'StripeColor','r');
  
  % Save drawn line parameters (two points):
  p1 = lineROI.Position(1,:);
  p2 = lineROI.Position(2,:);
  points(i,:) = [p1, p2];
  
  % Compute line from points (y = mx + b) NOT NEEDED
%   m = (y2 - y1) / (x2 - x1);
%   b = y1 - m * x1;
%   line_vec(i,:) = [m, b];
%   
end
