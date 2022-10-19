function [X, Y] = get_ascan_coordinates(header)
%GET_ASCAN_COORDINATES Compute A-scan coordinates (X, Y)
%
%   [X, Y] = get_ascan_coordinates(header)
%   Uses the information regarding the scanning protocol to compute A-scan
%   coordinate grid.
%
%   Input arguments:
%  
%   'header'         File header obtained when reading an OCT file.
%  
%
%   Output arguments:
%  
%   'X'              X coordinates pointing temporal to nasal.
%
%   'Y'              Y coordinates pointing inferior to superior.
%
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2022

switch header.bscan_pattern
    case 'raster'
        x_max = (header.scale_x * (header.n_ascan - 1)) /2;
        y_max = (header.scale_y * (header.n_bscan - 1)) /2;

        % We asume segmentation and bscan data is already coded so that
        % superior is the first row.
        y_range = linspace(y_max, -y_max, header.n_bscan);

        % OD: already pointing nasal
        % OS: need to flip
        x_range = linspace(-x_max, x_max, header.n_ascan);
        if strcmp(header.eye, 'OD')
            x_range = -x_range;
        end
        [X, Y] = meshgrid(x_range, y_range);            
    otherwise
        X = [];
        Y = [];
        warning('Unable to compute A-scan coordinates for this pattern');
end