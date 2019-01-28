function plot_camera2( P, w, h, scale )
%PLOT_CAMERA Plot a camera as a pyramid.

if nargin == 3, scale = 1; end

o = optical_center(P);
p1 = o + view_direction(P, [0 0]) * scale;
p2 = o + view_direction(P, [w 0]) * scale;
p3 = o + view_direction(P, [w h]) * scale;
p4 = o + view_direction(P, [0 h]) * scale;

vgg_scatter_plot2([o,p1],'k');
hold on;
vgg_scatter_plot2([o,p2],'k');
vgg_scatter_plot2([o,p3],'k');
vgg_scatter_plot2([o,p4],'k');
vgg_scatter_plot2([o,(p1+p2)/2],'k');
vgg_scatter_plot2([p1,p2],'k');
vgg_scatter_plot2([p2,p3],'k');
vgg_scatter_plot2([p3,p4],'k');
vgg_scatter_plot2([p4,p1],'k');
axis equal;

end

