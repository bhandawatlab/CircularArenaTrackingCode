function [xPts,yPts] = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xPts = r * cos(th) + x;
yPts = r * sin(th) + y;

end