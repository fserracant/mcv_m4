% Plot the homogeneous line l (3 dimensional array) in yellow

function plot_homog_line(l)

if(abs(l(1)) > abs(l(2))) % line is more vertical
    ylim = get(get(gcf,'CurrentAxes'), 'Ylim');
    l1 = [0 1 ylim(1)]'; % line y = ylim(1)    
    l2 = [0 1 ylim(2)]'; % line y = ylim(2)    
else
    xlim = get(get(gcf,'CurrentAxes'), 'Xlim');
    l1 = [1 0 -xlim(1)]'; % line x = xlim(1)    
    l2 = [1 0 -xlim(2)]'; % line x = xlim(2)
end

P1 = cross(l, l1); % intersection point of l with l1
P2 = cross(l, l2); % intersection point of l with l2

P1 = P1/P1(3);
P2 = P2/P2(3);
    
% plot the line that joins points P1 and P2
hl = line([P1(1); P2(1)],[P1(2); P2(2)]);
set(hl,'color','y');
