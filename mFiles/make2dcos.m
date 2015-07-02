% make2dcos.m
% TCS - 4/3/12

function z = make2dcos(x,y,x_center,y_center,r,pow)
% x, y is a meshgrid of x, y values at which to compute the 2d cos
% x_center, y_center is the center of the function
% r is the distance from center to 0 (T/2) - function will go from z = 0 to
% 0 across 2*r at widest point
% pow is power to which cosine function is raised (7, usually)

myr = ((x-x_center).^2+(y-y_center).^2).^0.5;

z = ((0.5*(1 + cos(myr*pi/r) )).^pow) .* (myr<=r)  ;


return