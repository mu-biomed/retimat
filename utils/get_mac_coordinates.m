function [X_fun,Y_fun,X_oct,Y_oct] = get_mac_coordinates(header, scale_x_fun, scale_y_fun)
% getAscanCoordinates - Compute A-Scan coordinates from .vol data
%
% [X_fun,Y_fun,X_ascan,Y_ascan] = get_ascan_coords(header, scale_x_fun, scale_y_fun);
%
% Input arguments:
%   header: struct with .vol header
%   scale_x_fun, scale_y_fun: fundus (slo) scaling calculated previously
% Output arguments:
%   X,Y: matrixes with cartesian A-Scan coordinates
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

% -----------------------------------------------------------------
% 1. Get original coordinates (as in .vol file)
% -----------------------------------------------------------------    
% Define fundus coordinate grid. Initial origin is at upper-left
% corner of the image and B-Scan delimiting values are in mm values
[X_fun,Y_fun] = meshgrid(linspace(0,scaleXSlo*(nX-1),nX),linspace(0,scaleYSlo*(nY-1),nY));

% Get B-Scan star-end coordinates
startX = header.StartX*adjust_x;
startY = header.StartY*adjust_y;

endX = header.EndX*adjust_x;
endY = header.EndY*adjust_y;        

% -----------------------------------------------------------------
% 2. Revert Y Axis (we want y axis in inferior-superior direction)
% -----------------------------------------------------------------        
Y_fun = -Y_fun;
startY = -startY;
endY = -endY;

% -----------------------------------------------------------------
% 3. Set coordinate origin at the center of the fundus image
% -----------------------------------------------------------------

% Distance to fundus image center
xOffset = scaleXSlo*(nX-1)/2;
yOffset = -scaleYSlo*(nY-1)/2;

% Translate fundus coordinates
X_fun = X_fun - xOffset;
Y_fun = Y_fun - yOffset;

% Translate start-end B-Scan coordinates
startX = startX - xOffset;
startY = startY - yOffset;
endX = endX -xOffset;
endY = endY -yOffset;

% Compute A-Scan coordinates
nAScans = header.NumAScans;
nBScans = header.NumBScans;

X_oct = nan(nBScans,nAScans);
Y_oct = nan(nBScans,nAScans);

for iBScan = 1:nBScans
    X_oct(iBScan,:) = linspace(startX(iBScan),endX(iBScan),nAScans);
    Y_oct(iBScan,:) = linspace(startY(iBScan),endY(iBScan),nAScans);
end

% -----------------------------------------------------------------
% 4. Set coordinate origin at the acquisition center (B-Scan 13)
% -----------------------------------------------------------------
% Get central B-Scan
midBScan = (nBScans+1)/2; % y nBScans=25 -> midBScan=13

% Check nBScans is odd
if mod(nBScans,2) == 0
    error('Even number of B-Scans in raster');
end

% Get Y coordinate of central B-Scan (average for small
% differences)
y_oct_center = mean(Y_oct(midBScan,:));

% Get X coordinate of central A-Scan             
midAScan = (nAScans+1)/2;
if mod(nAScans,2) ~=0 % If odd get the central directly
    x_oct_center = X_oct(midBScan,midAScan);
else % if odd there is no center point (find it between previous and posterior)
    xc_prev = X_oct(midBScan,floor(midAScan));
    xc_post = X_oct(midBScan,ceil(midAScan));            
    x_oct_center = (xc_prev + xc_post)/2;
end

% Translate all coordinates
X_fun = X_fun - x_oct_center;
Y_fun = Y_fun - y_oct_center;

X_oct = X_oct - x_oct_center;
Y_oct = Y_oct - y_oct_center;
  



