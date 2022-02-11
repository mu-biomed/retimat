function [Rc, x_cen, y_cen, tilt] = retinal_curvature(x_bm, y_bm)

% Avoid nan
mask = ~isnan(x_bm) & ~isnan(y_bm);

% Fit a circle
[x_cen, y_cen, Rc] = circle_fit(x_bm(mask)',y_bm(mask)');

% Tilting angle
tilt = atand(x_cen/y_cen);    
