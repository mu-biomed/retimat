function X = get_morph_params(rho, Z, parameters, average)
% Compute morphological parameters of the foveal pit
%
%
% Input arguments
% --------------- 
% * **rho**:         Matrix with rho coordinates (polar coordinates). Each row is a different angular direction.
%  
% * **Z**:           Matrix with thickness profile values. Each row is a different angular direction.
%            
% * **parameters**:  String or cell array of strings with the name of the parameters to computed. By default all available parameters are returned.
%  
% * **average**:     If true then radial parameters such as max slope are averaged across angular directions. Default is true.
%
%
% Output arguments
% ---------------- 
% * **X**:           Struct with computed features.          
%  
%  
% Notes
% -----
% This parameters are usually computed for the TRT profile but can be 
% computed as well for layers with a convex thickness profile such as GCL,
% IPL or INL
%
%
% References
% ----------
% [1] Scheibe P. et al., Parametric Model for the 3D Reconstruction of
% Individual Fovea Shape from OCT Data, Experimental Eye Research, 2014.
%
% [2] Yadav S. et al., CuBe: Parametric Modeling of 3D Foveal Shape Using 
% Cubic Bézier, Biomedical Optics Express, 2017.

% Feature names with flags for necessary previous steps
% Name    [cft, rim_height, rim_radius, slope, max_slope, max_slope_rad]

param_data = {'cft',                   1, 0, 0, 0, 0, 0;...
              'max_slope',             0, 1, 1, 1, 1, 0;...
              'max_slope_radius',      0, 1, 0, 1, 1, 1;...   
              'max_slope_height',      0, 1, 0, 1, 1, 0;...   
              'max_slope_disk_area',   0, 1, 0, 1, 1, 1;...   
              'max_slope_disk_perim',  0, 1, 0, 1, 1, 1;...   
              'mean_slope',            0, 1, 0, 1, 0, 0;...   
              'min_height',            0, 0, 0, 0, 0, 0;...
              'pit_area',              1, 1, 0, 0, 0, 0;...
              'pit_depth',             1, 1, 0, 0, 0, 0;...    
              'pit_volume',            1, 1, 1, 0, 0, 0;...
              'rim_height',            0, 1, 0, 0, 0, 0;...
              'rim_radius',            0, 1, 1, 0, 0, 0;...
              'rim_disk_area',         0, 1, 1, 0, 0, 0;...
              'rim_disk_perim',        0, 1, 1, 0, 0, 0};             

if nargin == 1
    error("The function expects at least 2 input arguments");
end
if nargin < 3
    parameters = param_data(:,1);
end
if nargin < 4
    average = true;
end

if ischar(parameters)
     if isequal(parameters, 'all')
         parameters = param_data(:,1);
     else
         parameters = {parameters}; 
     end
elseif ~iscell(parameters)
     error("features must be either a string or a cell of strings");
end

flags = zeros(1, 6);
for i=1:length(parameters)
    idx = find(strcmp(parameters{i}, param_data(:,1)));
    
    if isempty(idx)
        error(['Unknown parameter: ' parameters{i}']);
    end    
    flags = flags + cell2mat(param_data(idx,2:end));
end

% Check both axis in mm
Z = convert_mm_um(Z, 'mm');
rho = convert_mm_um(rho, 'mm');

n_angle = size(Z, 1);
step = rho(1,2) - rho(1,1);
theta = linspace(0, 2*pi, n_angle + 1); theta(end) = [];
        
% Common computations to avoid repetitions
if flags(1) ~= 0
    cft = Z(:, 1)';    
end
if flags(2) ~= 0  
    [rim_height, idx_rim] = max(Z.');
end
if flags(3) ~= 0  
    rim_radius = get_2d_points(rho, 1:n_angle, idx_rim);
end
if flags(4) ~= 0
    Zd = diff(Z, [], 2)./step;  % 1st derivative
end
if flags(5) ~= 0
    idx_maxslope = nan(1, n_angle);
     for n=1:n_angle
        [~, idx_maxslope(n)] = max(Zd(n, 1:idx_rim(n)-1));  
     end     
end
if flags(6) ~= 0
    % when computing the derivative we lose one unit (idx_maxslope + 1)
     max_slope_radius = get_2d_points(rho, 1:n_angle, idx_maxslope + 1);             
end

% Feature computation
for i=1:length(parameters)
    switch parameters{i}
        case 'cft'  % See [1]
            % Central foveal thickness. Should be equal for all directions
            X.cft = 1e3*cft(1);

        case 'max_slope_height'  % See [2]
            % When computing the derivative we lose one unit (idx_maxslope + 1)
            X.max_slope_height = 1e3*get_2d_points(Z, 1:n_angle, idx_maxslope + 1);                    

        case 'max_slope_radius'  % See [2]
            X.max_slope_radius = max_slope_radius;

        case 'max_slope'  % See [1]
            max_Zd = get_2d_points(Zd, 1:n_angle, idx_maxslope);
            X.max_slope = rad2deg(atan(max_Zd));
        
        case 'max_slope_disk_area'  % See [2]    
            [x, y] = pol2cart(theta, max_slope_radius);
            X.max_slope_disk_area = polyarea(x, y);

        case 'max_slope_disk_perim'  % See [2]
            [x, y] = pol2cart(theta, max_slope_radius);
            X.max_slope_disk_perim = perimeter(x, y);

        case 'mean_slope'        
            for n=1:n_angle
                slope = rad2deg(atan(Zd(n,1:idx_rim(n)-1))); % In degrees
                X.mean_slope(n) = mean(slope);
            end          
    %                 X.slope_mean(n) = 180*atan(mean(slope))/pi;            
        case 'min_height'
            % Might not correspond to the pit if there is a smaller minima at
            % large eccentricities
            X.min_height = min(1e3*Z, [], 2)';
                
        case 'pit_depth'  % See [1]
            X.pit_depth = 1e3*(rim_height - cft);
                    
        case 'pit_area'   % See [1]
            area_square = step * rim_height .* (idx_rim - 2);
            AUC = cft*step/2 + sum(Z(2:idx_rim-1)*step) + Z(idx_rim)*step/2;
            X.pit_area = area_square - AUC;

        case 'pit_volume' % See [2] (not exactly the same implementation)
            % Computed as the difference between all the volume under the
            % rim (Vt) and the volume under the foveal surface (Vs).
            Vt = nan(1,n_angle);
            Vs = nan(1,n_angle);
            
            % Total rim volume 
            % Computed for each angle as a 3d polygon defined by a base and
            % top triangles with coordinates (x,y,z_t)
            % We need to add 1 to compute last sector
            z = [rim_height rim_height(1)];
            rad = [rim_radius rim_radius(1)];
            th = [theta theta(1)];
            
            for n=1:n_angle
                % Lower triangle + upper triangle
                z_t = [0 0 0 mean(z) z(n) z(n+1)];
                rho2 = [0 rad(n) rad(n+1) 0 rad(n) rad(n+1)];
                th2 = [th(n) th(n) th(n+1) th(n) th(n) th(n+1)];
                [x, y] = pol2cart(th2, rho2);
                [k, Vt(n)] = convhull(x, y, z_t);
                
%                 trisurf(k,x,y,z_t);hold on;
            end
            
            % Volume under the foveal surface (Vs)
            % The formula comes from the diskreted volume integral defined
            % in equation 30-31 in [2]. We substitute c by Z and r by rho.
            for n=1:n_angle
                Vs(n) = sum(Z(n, 1:idx_rim(n)) .* (2*pi/n_angle).*rho(1:idx_rim(n)));
            end
                                    
            % Pit volume
            X.pit_volume = sum(Vt) - sum(Vs);
            
        case 'rim_height'  % See [1]
            X.rim_height = 1e3 * rim_height;
            
        case 'rim_radius'  % See [1]   
            X.rim_radius = rim_radius;

        case 'rim_disk_area'  % See [2]       
            [x_rim, y_rim] = pol2cart(theta, rim_radius);
            X.rim_disk_area = polyarea(x_rim, y_rim);
            
        case 'rim_disk_perim'  % See [2]                   
            [x_rim, y_rim] = pol2cart(theta, rim_radius);
            X.rim_disk_perim = perimeter(x_rim, y_rim);

        otherwise
            error(['Unknown parameter: ' parameters{i}]);
    end
end

if average
    for i=1:length(parameters)
        X.(parameters{i}) = mean(X.(parameters{i}));
    end
end

function p = perimeter(x, y)
% Compute perimeter based on pairwise distances between vertices
d = sqrt((x - [x(2:end) x(1)]).^2 + (y - [y(2:end) y(1)]).^2);
p = sum(d);        

function vals = get_2d_points(A, rows, cols)
% Access a set of 2d points efficiently
idx = sub2ind(size(A), rows, cols);  
vals = A(idx);   
