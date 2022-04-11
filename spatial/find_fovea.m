function [x_fovea, y_fovea] = find_fovea(X, Y, TRT, method)
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
%                    Default is: 'smooth_min'
%                    Options: ['none', 'min', 'resample_min', 'smooth_min']
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

switch method
    case 'none'
        x_fovea = 0;
        y_fovea = 0;
    
    case 'min'
        % Search only in the central region
        roi_radius = 0.85; 
        [~, rho] = cart2pol(X, Y);
        TRT(rho > roi_radius) = Inf;        
        [x_fovea, y_fovea] = find_min(X, Y, TRT);
        
     case 'resample_min'
        % Interpolate to a small regular grid 
        [X, Y, TRT] = resample_map(X, Y, TRT, 'regular', 'max_d', 0.85, ...
                                   'n_point', 50);        
        [x_fovea, y_fovea] = find_min(X, Y, TRT);
         
    case 'smooth_min'        
        % Interpolate to small regular grid 
        [X, Y, TRT] = resample_map(X, Y, TRT, 'regular', 'max_d', 0.85, ...
                                   'n_point', 50);
        
        % Smooth the TRT map with a circular kernel of 0.25 mm radius 
        kernel_rad = 0.25;
        [~, rho] = cart2pol(X, Y);
        kernel = rho < kernel_rad;
        TRT_filt = imfilter(TRT, double(kernel)/sum(kernel(:)),'replicate');

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