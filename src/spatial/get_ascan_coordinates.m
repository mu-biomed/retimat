function [X, Y] = get_ascan_coordinates(header)
% Compute A-scan coordinates (X, Y) based on the scanning protocol
%
%
% Input arguments
% --------------- 
% * **header**: File header obtained when reading an OCT file.
% 
%
% Output arguments
% ---------------- 
% * **X**:      2D matrix with X coordinates pointing temporal to nasal.
%
% * **Y**:      2D matrix with Y coordinates pointing inferior to superior.

switch header.bscan_pattern
    case 'raster'
        x_max = (header.scale_x * (header.n_ascan - 1)) /2;
        y_max = (header.scale_y * (header.n_bscan - 1)) /2;

        x_range = linspace(-x_max, x_max, header.n_ascan);
        % We asume segmentation and bscan data is already coded so that
        % superior is the first row.
        y_range = linspace(y_max, -y_max, header.n_bscan);
        
        % Deprecated: OS flipping to have x > 0 pointing nasal
        % Better leave this step for later processing if needed.
        %  if strcmp(header.eye, 'OS')
        %      x_range = -x_range;
        %  end
        
        [X, Y] = meshgrid(x_range, y_range);         
    case 'star'
        % Assumptions:
        % - 1st bscan vertical
        % - Next bscans clockwise
        theta = linspace(pi/2, -pi/2, header.n_bscan + 1);
        theta = theta(1:end-1).' * ones(1, header.n_ascan);
        
        size_x = (header.scale_x * (header.n_ascan - 1)) / 2;
        rho = repmat(-size_x:header.scale_x:size_x, header.n_bscan, 1);
        
        [X, Y] = pol2cart(theta, rho);
    otherwise
        X = [];
        Y = [];
        warning('Unable to compute A-scan coordinates for this pattern');
end