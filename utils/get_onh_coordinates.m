function [X_fun,Y_fun,X_oct,Y_oct,x_onh,y_onh] = get_onh_coordinates(header,scale_x_fun, scale_y_fun)
% get onh coordinates - Compute A-Scan coordinates of onh from .vol data
%
% [X_fun,Y_fun,X_oct,Y_oct,x_onh,y_onh] = get_ascan_coords(header)
%
% Input arguments:
%   header: struct with .vol header
%
% Output arguments:
%   X_fun:
%   Y_fun:
%   X_oct:
%   Y_oct:
%   x_onh:
%   y_onh:
%
% David Romero-Bascones 
% dromero@mondragon.edu
% 2021, Mondragon Unibertsitatea, Biomedical Engineering Department
% ------------------------------------------------------------------------

% Global properties of fundus image 
nX = double(header.SizeXSlo);
nY = double(header.SizeYSlo);

scaleXSlo = scale_x_fun;
scaleYSlo = scale_y_fun;

% Adjusting factor to correct Spectralis scaling
adjust_x = scale_x_fun/header.ScaleXSlo;
adjust_y = scale_y_fun/header.ScaleYSlo;

nAScan = header.NumAScans;

eye = header.Eye(1:2);


% -----------------------------------------------------------------
% 1. Get original coordinates (as in .vol file)
% -----------------------------------------------------------------        

% Starting A-Scan
startX = header.StartX(1)*adjust_x;
startY = header.StartY(1)*adjust_y;

% Fundus coordinates
[X_fun,Y_fun] = meshgrid(linspace(0,scaleXSlo*(nX-1),nX),linspace(0,scaleYSlo*(nY-1),nY));

% ONH center
x_onh = header.EndX(1)*adjust_x;
y_onh = header.EndY(1)*adjust_y;

% Revert Y axis
Y_fun = -Y_fun;
y_onh = -y_onh;
startY = -startY;

% Compute distance to the fundus image centre
x_offset = scaleXSlo*(nX-1)/2;
y_offset = -scaleYSlo*(nY-1)/2;

% Translate to image centre
startX = startX - x_offset;
startY = startY - y_offset;
X_fun = X_fun - x_offset;
Y_fun = Y_fun - y_offset;
x_onh = x_onh - x_offset;
y_onh = y_onh - y_offset;

% Compute B-Scan coordinates

% Compute onh radius
radius = sqrt((y_onh-startY)^2 + (x_onh-startX)^2);
% Generate a circle (clockwise or anticlockwise depending on eye type)
if strcmp(eye,'OS')
    [X_oct,Y_oct] = pol2cart(linspace(0,2*pi,nAScan),repmat(radius,1,nAScan));
elseif strcmp(eye,'OD')
    [X_oct,Y_oct] = pol2cart(linspace(pi,-pi,nAScan),repmat(radius,1,nAScan));
end

% Translate the circle to the onh center
X_oct = X_oct + x_onh;
Y_oct = Y_oct + y_onh;
% 
% surf(X_fun,Y_fun,fundus,'EdgeColor','none')
% view(0,90)
% colormap(gray)
% hold on;
% scatter3(x_onh,y_onh,max(double(fundus(:))),'r');
% scatter3(X_oct,Y_oct,max(double(fundus(:)))*ones(1,nAScan),'b');