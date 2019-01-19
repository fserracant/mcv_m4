function plot_camera( P, w, h, scale, colour )
%PLOT_CAMERA Plot a camera as a pyramid.
% colour: colour of the drawing: see plot or plot3 function help

if nargin == 3
  scale = 1; 
  colour = 'b';
end

if nargin == 4
  colour = 'b';
end

o = optical_center(P);
p1 = o + view_direction(P, [0 0]) * scale;
p2 = o + view_direction(P, [w 0]) * scale;
p3 = o + view_direction(P, [w h]) * scale;
p4 = o + view_direction(P, [0 h]) * scale;

vgg_scatter_plot([o,p1],colour);
hold on;
vgg_scatter_plot([o,p2],colour);
vgg_scatter_plot([o,p3],colour);
vgg_scatter_plot([o,p4],colour);
vgg_scatter_plot([o,(p1+p2)/2],colour);
vgg_scatter_plot([p1,p2],colour);
vgg_scatter_plot([p2,p3],colour);
vgg_scatter_plot([p3,p4],colour);
vgg_scatter_plot([p4,p1],colour);
axis equal;

end

