function X = get_morph_params(rho, Z, features, average)
% GET_MORPH_PARAMS - compute morphological parameters of the foveal pit
% based on thickness profile values
%
% Output arguments:
%   Paramorf: struct with morphological parameters
%
%   Paramorf = get_morph_params(rho, Z, features, average)
%   Detail explanation goes here
%
%   Input arguments:
%  
%   'rho'            Matrix with rho coordinates (polar coordinates). Each row
%                    is a different angular direction.
%  
%   'Z'              Matrix with thickness profile values. Each row is a
%                    different angular direction.
%            
%   'features'       String or cell array of strings with the name of the
%                    parameters to computed. By default all available 
%                    parameters are returned.
%  
%   'average'        If true then radial parameters such as max slope are
%                    averaged across angular directions. Default is true.
%
%
%   Output arguments:
%  
%   'X'              Struct with computed features.          
%  
%
%   
%   Notes
%   -----
%   This parameters are usually computed for the TRT profile but can be 
%   computed as well for layers with a convex thickness profile such as GCL,
%   IPL or INL
%
%
%   References
%   ----------
%   [1] 
%
%   Example 1
%   ---------      
%   % Example description
%
%     I = [1 1 5 6 8 8;2 3 5 7 0 2; 0 2 3 5 6 7];
%     [GLCMS,SI] = graycomatrix(I,'NumLevels',9,'G',[])
%     
%
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2021

all_features =  {'cft',...
                'height_min',...
                'height_rim',...
                'height_maxslope',...
                'pit_area',...
                'pit_depth',...
                'radius_rim',...
                'radius_maxslope',...
                'slope_max',...
                'slope_mean'};
                         
if nargin == 1
    error("The function expects at least 2 input arguments");
end
if nargin < 3
    features = all_features;
end
if nargin < 4
    average = true;
end

if ischar(features)
     if isequal(features, 'all')
         features = all_features;
     end
elseif ~iscell(features)
     error("features must be either a string or a cell of strings");
end
    
% Get step
n_angle = size(Z, 1);

step = rho(1,2) - rho(1,1);

% Check both axis in mm
if max(Z(1))>1
   Z = 1e-3*Z;
   warning("Input thickness seems to be in um. Converting to mm during computation");
end

% Common computations
% if intersect(features, {'height_rim', 'height_max_slope', 'pit_depth',...}
        
cft = 1e3*Z(:, 1);
[Z_max, ind_max] = max(Z, [], 2);  % find rim height to limit search roi
Zd = diff(Z, [], 2)./step;  % 1st derivative

for i=1:length(features)
    switch features{i}
        case 'cft'
            % Central foveal thickness. Should be equal for all directions
            X.cft = cft;
            
        case 'height_rim'
            % Rim height
            X.height_rim = 1e3 * Z_max;
            
        case 'height_min'
            % Minimum height
            X.height_min = min(1e3*Z, [], 2);
            
        case 'height_maxslope'
            for n=1:n_angle
                [~, ind_maxd] = max(Zd(n, 1:ind_max(n)));  % Steepest point
                ind_maxd = ind_maxd(1);  % In case there are two maxima
                % when computing the derivative we lose one unit (ind_maxd + 1)
                X.max_slope_height(n) = 1e3*Z(n, ind_maxd + 1);                
            end
            
        case 'pit_depth'
            X.pit_depth = Z_max - cft;
            
        case 'pit_area'           
            area_square = step * Z_max .* (ind_max - 2);
            AUC = cft*step/2 + sum(Z(2:ind_max-1)*step) + Z(ind_max)*step/2;
            X.pit_area = area_square - AUC;
    
        case 'radius_rim'
            for n=1:n_angle
                X.radius_rim(n) = rho(n, ind_max(n));
            end
            
        case 'radius_maxslope'
            for n=1:n_angle
                [~, ind_maxd] = max(Zd(n, 1:ind_max(n)));  % Steepest point
                ind_max = ind_maxd(1);  % In case there are two maxima
                % when computing the derivative we lose one unit (ind_maxd + 1)
                X.max_slope_radius(n) = rho(n, ind_maxd + 1);
            end   
            
        case 'slope_max'
            for n=1:n_angle
                max_Zd = max(Zd(n, 1:ind_max(n)));  % Steepest point
                max_Zd = max_Zd(1);  % In case there are two maxima
                X.slope_max(n) = rad2deg(atan(max_Zd));
            end
            
        case 'slope_mean'
            for n=1:n_angle
                slope = rad2deg(atan(Zd(n,1:ind_max(n)))); % In degrees
                X.slope_mean(n) = mean(slope);
            end          
%                 X.slope_mean(n) = 180*atan(mean(slope))/pi;            
              
        otherwise
            error("Unknown feature provided");
    end
end

if average
    for i=1:length(features)
        if length(X.(features{i})) > 1
            X.(features{i}) = mean(X.(features{i}));
        end
    end
end
