function [X1, Y1, Z1] = resample_map(X, Y, Z, grid_type, varargin)
%RESAMPLE_MAP Resample a 2D map into a new grid by interpolation/extrapolation
%
%   [X1, Y1, Z1] = resample_map(X, Y, Z, grid_type, varargin)
%   Rarange a thickness map into a new grid by means of interpolation and
%   extrapolation.
%
%   Input arguments (mandatory):
%  
%   'X'              Original X coordinates.
%
%   'Y'              Original Y coordinates.
%
%   'Z'              Original Z values (e.g., thickness values).
%            
%   'grid_type'      String specifying the resampling grid type. Depending on 
%                    grid type it must be followed by specific extra arguments.
%                    Options: ['regular', 'star']
%
%
%   Input arguments (name/value pair)
%
%   'n_point'        Meaning for each grid_type:
%                    - 'regular': number of points in X and Y directions
%                    - 'star': number points from 0 to max_d
%
%   'max_d'          Distance to cover. Meaning for each grid_type:
%                    - 'regular': maximum X and Y values
%                    - 'star': maximum radius
%
%   'n_angle'        Number of directions to use (when the grid_type is 'star')
%                    
%
%   'interp_method'  interpolation method ['linear','cubic']
%                    Default: 'linear'
%   
%   'extrapolate'    If true then points with NaN outside the original range 
%                    are extrapolated.
%                    Default: false  
%  
%
%   Output arguments:
%  
%   'X1'             New X coordinates.
%
%   'Y1'             New Y coordinates.
%
%   'Z1'             New Z values.
%  
%
%   
%   Notes
%   -----
%   If the region to cover is bigger than the original data region it might not
%   be possible to extrapolate values accurately.
%
%
%
%   Example
%   ---------      
%   % Resample original a-scans to a regular grid
%
%   [header, seg, ~, ~] = read_vol(myfile.vol, 'coordinates');
%   Thickness = compute_thickness(seg, 'TRT', header.scale_z);
%
%   [X, Y, TRT] = resample_map(header.X_oct, header.Y_oct, Thickness.TRT, ...
%   'regular', 'n_point', 100, 'max_d', 2.5);
%
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2022.

if mod(nargin,2) ~=0
    error('Number of arguments must be even');
end

n_point = nan;
max_d = nan;
n_angle = nan;
interp_method = 'linear';
extrapolate = false;

% Process varargin
for i=1:2:length(varargin)-1
    switch varargin{i}
        case 'n_point'
            n_point = varargin{i+1};
        case 'max_d'
            max_d = varargin{i+1};
        case 'n_angle'
            n_angle = varargin{i+1};
        case 'interp_method'
            interp_method = varargin{i+1};
        case 'extrapolate'
            extrapolate = varargin{i+1};
        otherwise
            error('Unknown name/value argument provided');
    end            
end

if isnan(n_point)
    error("You must provide a valid value for 'n_point'");
end
if isnan(max_d)
    error("You must provide a valid value for 'max_d'");
end

%  Generate the resampling grid in both cartesian and polar coordinates
switch grid_type    
    case 'regular'
        x_grid = linspace(-max_d, max_d, n_point);
        y_grid = x_grid;
        [X1, Y1] = meshgrid(x_grid, y_grid);  
        
    case 'star'        
        if isnan(n_angle)
            error("You must provide a valid value for 'n_angle'");
        end
        
        Rho = linspace(0, max_d, n_point);
        Theta = linspace(0, 2*pi, n_angle+1);
        Theta(end) = [];
        Rho = repmat(Rho, n_angle,1);
        Theta = repmat(Theta', 1, size(Rho, 2));
        
        [X1, Y1] = pol2cart(Theta, Rho);   
    otherwise
        error("Uknown grid type. Use 'regular' or 'star'.");
end

% Interpolate over new grid
is_num = ~isnan(Z); % mask of not nan     
    
% Interpolation (griddata does not extrapolate)  
Z1 = griddata(X(is_num), Y(is_num), Z(is_num), X1(:), Y1(:), interp_method);
Z1 = reshape(Z1, size(X1));

% Extrapolation of remaining NaNs (for outside sampled region)
% Caution not to extrapolate too far from the sampled region as
% the accuracy might be very poor
% scatteredInterpolant does not support cubic so it is linear           
if extrapolate
    is_num = ~isnan(Z1); % mask of not nan                
    if sum(~is_num(:))>0
        warning('Extrapolation necessary');
        interpol = scatteredInterpolant(X1(is_num), Y1(is_num), Z1(is_num), 'linear');
        Z1 = reshape(interpol(X1(:), Y1(:)), size(X1));                 
    end
end