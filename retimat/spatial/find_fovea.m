function [x_fovea, y_fovea] = find_fovea(X, Y, TRT, method, max_d)
%FIND_FOVEA Find foveal center based on total retinal thickness (TRT) map
%
%   [x_fovea, y_fovea] = template(X, Y, TRT, method)
%   Returns x and y coordinates of the foveal center
%
%   Input arguments:
%  
%   'X'              Matrix with X coordinates of each A-Scan.
%
%   'Y'              Matrix with Y coordinates of each A-Scan.
%
%   'TRT'            Matrix with total retinal thickness values.            
%  
%   'method'         Method used to find the foveal center.
%                    Default: 'smooth_min'
%                    Options: ['none', 'min', 'resample_min', 'smooth_min']
%   
%   'max_d'          Maximum alignment error.
%                    Default: 0.85
%
%
%   Output arguments:
%  
%   'x_fovea'        X coordinate of the foveal center.          
%  
%   'y_fovea'        Y coordinate of the foveal center.
%   
%
%   Notes
%   -----
%   The smooth+min method is based on the foveafinder.m function
%   of AURA Tools. If you use it, please provide appropriate credit to the
%   original work (https://www.nitrc.org/projects/aura_tools/).
%
%
%   References
%   ----------
%   [1] Romero-Bascones et al., Foveal Pit Morphology Characterization: A 
%   Quantitative Analysis of the Key Methodological Steps, Entropy, 2021
%   https://doi.org/10.3390/e23060699
%   
%   Example 
%   ---------      
%   % Example description
%
%     file = '../data/raster.vol';
%     [header, seg, ~, ~] = read_vol(file, 'coordinates');
%     Thickness = compute_thickness(seg, 'TRT', header.scale_z);
%     [x_fovea, y_fovea] = find_fovea(header.X, header.Y, Thickness.TRT)
%  
%
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2021

if nargin == 3
    method = 'smooth_min';
end
if nargin <= 4
    max_d = 0.85;
end

switch method
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
        
    case 'robust'
        % Interpolate to small regular grid 
        [X, Y, TRT] = resample_map(X, Y, TRT, 'regular', 'max_d', max_d, ...
                                   'n_point', 100);
        
        % Find the ring with the highest mean thickness and define a search
        % region as any point inside it.
        interpol = scatteredInterpolant(X(:), Y(:), TRT(:));

        n_r = 5;
        n_x = 10;
        n_y = 10;
        
        r_range = linspace(1,  1.5, n_r);
        x_range = linspace(-1, 1,   n_x);
        y_range = linspace(-1, 1,   n_y);
        
        metric = nan(n_r, n_x, n_y);
        
        n_point = 20;
        theta = linspace(0, 2*pi, n_point);
        
        for i=1:n_x
            for j=1:n_y
                for k=1:n_r
                    [x, y] = pol2cart(theta, r_range(k)*ones(1,n_point));
                    x = x - x_range(i);
                    y = y - y_range(j);            
                    metric(i, j, k) = mean(interpol(x,y),'omitnan');            
                    
        %             x = X - x_range(i);
        %             y = Y - y_range(j);
        %             mask = sqrt(x.^2 + y.^2) < r_range(k);
        %             metric(i,j,k) = mean(Z(mask),'omitnan');
                end
            end
        end

        [~, idx] = max(metric(:));
        [i,j,k] = ind2sub(size(metric), idx);

        x_opt = x_range(i);
        y_opt = y_range(j);
        r_opt = r_range(k);
        
%         close all;
%         surf(X,Y,TRT,'EdgeColor','none');view(0,90);hold on;
%         [x,y] = pol2cart(linspace(0,2*pi,20), r_opt*ones(1,20));
%         x = x - x_opt;
%         y = y - y_opt;
%         plot3(x,y,max(TRT(:))*ones(1,20),'--r');

        search_mask = sqrt((X + x_opt).^2 + (Y + y_opt).^2) < r_opt;
        
        % Smooth the TRT map with a circular kernel of 0.25 mm radius 
        kernel_rad = 0.25;
        [~, rho] = cart2pol(X, Y);
        kernel = rho < kernel_rad;
        kernel = double(kernel)/sum(kernel(:));

        TRT_filt = imfilter(TRT, kernel,'replicate');

        TRT_filt(~search_mask) = nan;
        [x_fovea, y_fovea] = find_min(X, Y, TRT_filt);

    otherwise
        error(strcat("Unsupported fovea location method. Valid options: ",...
            "'none','min','resample_min','smooth_min'"));
end

function [x_min, y_min] = find_min(X, Y, Z)
[~, ind_min] = min(Z(:));
[ind_x, ind_y] = ind2sub(size(X), ind_min);        
x_min = X(ind_x, ind_y);
y_min = Y(ind_x, ind_y);        