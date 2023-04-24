function Zs = sectorize_peripapillar(X, Y, Z, sect_type, varargin)
% Sectorize thickness values of a peripapillar circular scan
%
%
% Input arguments
% --------------- 
% * **X**:              X coordinates of map points
%
% * **Y**:              Y coordinates of map points
%
% * **Z**:              Z coordinates of map points
%
% * **sector_type**:    String defining the sectorization type:
%
%   - 'average'
%   - '4_quadrant'
%   - 'qustom': must be followed by the number of angles and the initial angle (see examples).
%   
% * **varargin**:       Extra arguments to define the sectorization. 
%                    
%  
% Output arguments
% ---------------- 
% * **Zs**:             Sectorized values for each sector.          
%
%   
% Notes
% -----
% X and Y coordinates must follow the temporal-nasal and inferior-superior
% conventions.
% 
% '4_quadrant' sectorization returns values in the order: 'nasal',
% 'superior', 'temporal' and 'inferior'.
%
%
% Examples
% --------      
% 4 quadrants sectorization
% ^^^^^^^^^^^^^^^^^^^^^^^^^
% .. code-block:: matlab
%
%   [header, seg, ~, ~] = read_vol(file,'verbose', 'get_coordinates');
%   Thickness = compute_thickness(seg, 'TRT', header.scale_z);
%   Z = sectorize_peripapillar(X, Y, TRT, '4_quadrant');
%
% Qustom number of quadrants sectorization
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% .. code-block:: matlab
%
%   [header, seg, ~, ~] = read_vol(file,'verbose', 'get_coordinates');
%   Thickness = compute_thickness(seg, 'TRT', header.scale_z);
%   Z = sectorize_peripapillar(X, Y, TRT, 'qustom', 8, -pi/8);

switch sect_type
    case 'average'
        Zs = mean(Z, 'omitnan');
    case '4_quadrant'
        n_angle = 4;
        theta_0 = -pi/4;
        Zs = average_sectors(X, Y, Z, n_angle, theta_0);
    case 'qustom'
        n_angle = varargin{1};
        theta_0 = varargin{2};
        Zs = average_sectors(X, Y, Z, n_angle, theta_0);
    otherwise
        error("Unknown sectorization");
end

function Zs = average_sectors(X, Y, Z, n_angle, theta_0)
    [theta, ~] = cart2pol(X, Y);
    
    % Rotate all angles to set to position initial theta as theta = 0
    theta = theta - theta_0;

    % Reconfigure to have only positive values
    theta(theta < 0) = theta(theta < 0) + 2*pi;
    
    theta_grid = linspace(0, 2*pi, n_angle + 1);
    Zs = nan(1, n_angle);
    for i_angle=1:n_angle
        mask = (theta >= theta_grid(i_angle)) & (theta < theta_grid(i_angle+1));
        Zs(i_angle) = mean(Z(mask), 'omitnan');
    end