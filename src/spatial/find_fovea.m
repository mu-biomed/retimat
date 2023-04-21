function [x_fovea, y_fovea] = find_fovea(X, Y, TRT, method, max_d)
% Find foveal center based on total retinal thickness (TRT) map.
%
%
% Input arguments
% --------------- 
% * **X**:              Matrix with X coordinates of each A-Scan.
%
% * **Y**:              Matrix with Y coordinates of each A-Scan.
%
% * **TRT**:            Matrix with total retinal thickness values.            
%
% * **method**:         Method used to find the foveal center. Options ['none', 'min', 'resample_min', 'smooth_min']
%  
% * **max_d**:          Maximum alignment error. Default: 0.85
%
%
% Output arguments
% ---------------- 
% * **x_fovea**:        X coordinate of the foveal center.          
%  
% * **y_fovea**:        Y coordinate of the foveal center.
%   
%
% Notes
% -----
% The ``smooth+min`` method is based on the ``foveafinder.m`` function of
% AURA tools. If you use it, please provide appropriate credit to the
% original work (https://www.nitrc.org/projects/aura_tools/).
%
%
% References
% ----------
% [1] Romero-Bascones et al., Foveal Pit Morphology Characterization: A 
% Quantitative Analysis of the Key Methodological Steps, Entropy, 2021
% https://doi.org/10.3390/e23060699
%   
% Example
% -------      
%
% .. code-block:: matlab
% 
%   file = '../data/raster.vol';
%   [header, seg, ~, ~] = read_vol(file, 'coordinates');
%   Thickness = compute_thickness(seg, 'TRT', header.scale_z);
%   [x_fovea, y_fovea] = find_fovea(header.X, header.Y, Thickness.TRT)

if nargin == 3
    method = 'smooth_min';
end
if nargin <= 4
    max_d = 0.85;
end

switch method
    case 'flood'
        n_seed   = 3000;
        max_step = 30;
        n_point  = 128;
        
        Z = TRT;
        % not needed
        % Z = round(255 * (Z - min(Z(:))) ./ (max(Z(:)) - min(Z(:))));
        
        hist2d = zeros(n_point, n_point);  % count of end pixels
        
        Z_flat = Z(:); % for indexing efficiency
        for i_seed=1:n_seed
            i = randi(n_point, 1);
            j = randi(n_point, 1);

            for i_step=1:max_step
                fprintf('%d/%d step %d\n',i_seed, n_seed, i_step);

                % Define suitable neighbors
                neigh_i = [i-1   i i+1 i-1 i+1 i-1   i i+1];
                neigh_j = [j-1 j-1 j-1   j   j j+1 j+1 j+1];

                include = neigh_i > 0 & neigh_j > 0 & neigh_i <= n_point & neigh_j <=n_point;
                if sum(include) == 0
                    break;
                end
                
                neigh_i = neigh_i(include);
                neigh_j = neigh_j(include);
                
                % Compute gradient to neighbors
                idx = sub2ind([n_point n_point], neigh_i, neigh_j);
                d = Z_flat(idx) - Z(i, j);

                % If no neighbors with smaller thickness
                if all(d(~isnan(d)) >= 0)
                    break;
                end

                % Select neighbor as new step
                [~, idx_min] = min(d);
                i = neigh_i(idx_min);
                j = neigh_j(idx_min);
            end

            % If end point is an edge --> remove
            if any([i j] == 1) || any([i j] == n_point)
                continue;
            end

            % Increase count for the reached end point
            idx = sub2ind([n_point n_point], i, j);

            hist2d(i, j) = hist2d(i, j) + 1;    
            Z(i, j) = Z(i, j) + 1;
            Z_flat(idx) = Z_flat(idx) + 1;
        end

        [~, idx_fov] = max(hist2d(:));
        [i_fov, j_fov] = ind2sub([n_point n_point], idx_fov);
        x_fovea = X(i_fov, j_fov);
        y_fovea = Y(i_fov, j_fov);
        
    case 'none'
        x_fovea = 0;
        y_fovea = 0;
    
    case 'min'
        % Search only in the central region
        roi_radius = max_d; 
        [~, rho] = cart2pol(X, Y);
        TRT(rho > roi_radius) = Inf;        
        [x_fovea, y_fovea] = find_min(X, Y, TRT);
        
     case 'resample_min'
        % Interpolate to a small regular grid 
        [X, Y, TRT] = resample_map(X, Y, TRT, 'regular', 'max_d', max_d, ...
                                   'n_point', 50);        
        [x_fovea, y_fovea] = find_min(X, Y, TRT);
         
    case 'smooth_min'        
        % Interpolate to small regular grid 
        [X, Y, TRT] = resample_map(X, Y, TRT, 'regular', 'max_d', max_d, ...
                                   'n_point', 50);
        
        % Smooth the TRT map with a circular kernel of 0.25 mm radius 
        kernel_rad = 0.25;
        [~, rho] = cart2pol(X, Y);
        kernel = rho < kernel_rad;
        TRT_filt = imfilter(TRT, double(kernel)/sum(kernel(:)),'replicate');

        [x_fovea, y_fovea] = find_min(X, Y, TRT_filt);
        
    otherwise
        error(strcat("Unsupported fovea location method. Valid options: ",...
            "'flood','none','min','resample_min','smooth_min'"));
end

function [x_min, y_min] = find_min(X, Y, Z)
[~, ind_min] = min(Z(:));
[ind_x, ind_y] = ind2sub(size(X), ind_min);        
x_min = X(ind_x, ind_y);
y_min = Y(ind_x, ind_y);        