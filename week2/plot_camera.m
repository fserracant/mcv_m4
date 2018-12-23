function plot_camera( P, w, h, scale )
%PLOT_CAMERA Plot a camera as a pyramid.

if nargin == 3, scale = 1; end

o = optical_center(P);
p1 = o + view_direction(P, [0 0]) * scale;
p2 = o + view_direction(P, [w 0]) * scale;
p3 = o + view_direction(P, [w h]) * scale;
p4 = o + view_direction(P, [0 h]) * scale;

vgg_scatter_plot([o,p1],'b');
hold on;
vgg_scatter_plot([o,p2],'b');
vgg_scatter_plot([o,p3],'b');
vgg_scatter_plot([o,p4],'b');
vgg_scatter_plot([o,(p1+p2)/2],'b');
vgg_scatter_plot([p1,p2],'b');
vgg_scatter_plot([p2,p3],'b');
vgg_scatter_plot([p3,p4],'b');
vgg_scatter_plot([p4,p1],'b');
axis equal;

end

