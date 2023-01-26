function [Zs, Sectors] = sectorize_map(X, Y, Z, metric, sector_type, varargin)
%SECTORIZE_MAP Sectorize a 2D map into several sectors
%
%   Zs = sectorize_map(X, Y, Z, metric, sector_type)
%   Summarize thickness point values across sectors. For example by taking
%   the average thickness for ETDRS sectorization.
%
%   Input arguments:
%  
%   'X'              X coordinates of map points
%
%   'Y'              Y coordinates of map points
%
%   'Z'              Z coordinates of map points
%
%   'metric'         Metric to be used to sectorize data. 
%                    Options: ['mean', 'std', skewness', 'kurtosis']
%
%   'sector_type'    Sectorization definition. Two options:
%                    - String defining the sectorization type. It must be 
%                    followed by appropriate arguments (see examples).
%                    - Struct with sectorization info created beforehand.
%   
%   'varargin'       Extra arguments to define the sectorization. 
%                    
%  
%   Output arguments:
%  
%   'Zs'             Sectorized values for each sector.          
%  
%   'Sectors'        Struct with spatial data defining the sectorization. 
%                    Useful when the Sectors input variable is a string.
%
%   
%   Notes
%   -----
%   All angles are expected in radians. To convert from degrees to radians 
%   use deg2rad() function.   
%
%
%
%   Example 1
%   ---------      
%   % ETDRS sectorization
%
%   [header, seg, ~, ~] = read_vol(file,'verbose', 'get_coordinates');
%   Thickness = compute_thickness(seg, 'TRT', header.scale_z);
%   [X, Y, TRT] = resample_map(header.X_oct, header.Y_oct, Thickness.TRT, ...
%        'regular', 'n_point', 100, 'max_d', 2.5);
%   [Z, G] = sectorize_map(X, Y, TRT, 'mean', 'etdrs');
%     
%
%   Example 2
%   ---------      
%   % 2 ring sectorization
%
%   [header, seg, ~, ~] = read_vol(file,'verbose', 'get_coordinates');
%   Thickness = compute_thickness(seg, 'TRT', header.scale_z);
%   [X, Y, TRT] = resample_map(header.X_oct, header.Y_oct, Thickness.TRT, ...
%        'regular', 'n_point', 100, 'max_d', 2.5);
%   [Z, G] = sectorize_map(X, Y, TRT, 'mean', 'ring', [0.5 1.5 3]);
%  
%
%
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2021

if nargin < 4
    error("A minimum of 4 input arguments must be provided");
end

if ischar(sector_type)    
    switch sector_type
        case 'regular'
            Sectors.type = sector_type;
            Sectors.n_x = varargin{1};
            Sectors.n_y = varargin{2};
            Sectors.n_sect = Sectors.n_x*Sectors.n_y;
            Sectors.X_edge = linspace(min(X(:)), max(X(:)), Sectors.n_x + 1);
            Sectors.Y_edge = linspace(min(Y(:)), max(Y(:)), Sectors.n_y + 1);        
        case 'disk'
            Sectors.type = sector_type;
            Sectors.radius = varargin{1};        
            Sectors.n_sect = 1;
        case 'ring'
            Sectors.type = sector_type;
            Sectors.radius = varargin{1};
            Sectors.n_sect = length(Sectors.radius)-1;
        case 'pie'
            Sectors.type = sector_type;
            Sectors.n_angle = varargin{1};
            Sectors.theta_0 = varargin{2};
            Sectors.n_sect = Sectors.n_angle;
        case 'wedge'
            Sectors.type = sector_type;
            Sectors.radius = varargin{1};
            Sectors.n_angle = varargin{2};            
            Sectors.theta_0 = varargin{3};
            Sectors.n_sect = 1 + (length(Sectors.radius)-1)*Sectors.n_angle;
        case 'etdrs'
            Sectors.type = 'wedge';  % etdrs is effectively a wedge sectorization
            Sectors.radius = [0.5 1.5 3];
            Sectors.theta_0 = -pi/4;
            Sectors.n_angle = 4;
            Sectors.n_sect = 9;
        otherwise
            error("Incorrect sector type. Options: 'regular','disk','ring','pie','wege','etdrs'");
    end    
elseif isstruct(sector_type)
    Sectors = sector_type;
else
    error("sector_type argument must be either a struct or a string");
end

Zs = nan(1, Sectors.n_sect);

switch metric
    case 'mean'
        fun = @(x) mean(x, 'omitnan');
    case 'median'
        fun = @(x) median(x, 'omitnan');
    case 'std'
        fun = @(x) std(x, 'omitnan');
    case 'skewness'
        fun = @(x) skewness(x);
    case 'kurtosis'
        fun = @(x) kurtosis(x);
    case 'snr'
        fun = @(x) mean(x, 'omitnan')/std(x, 'omitnan');
    otherwise
        error(strcat("Unknown metric. Accepted values: 'mean','median',",...
            "'std','skewness','kurtosis','snr'"));
end

switch Sectors.type
    case 'regular'        
        % Get grid 
        X_edge = Sectors.X_edge;
        Y_edge = Sectors.Y_edge;

        % Add small value to capture all Z points despite the < sign
        X_edge(end) = X_edge(end) + 0.1;
        Y_edge(end) = Y_edge(end) + 0.1;
        
        % Compute average per sector
        i_sect = 1;
        for i=1:Sectors.n_x     
            for j=1:Sectors.n_y
         
                % Compute sector mask
                mask_x = X>=X_edge(i) & X<X_edge(i+1);
                mask_y = Y>=Y_edge(j) & Y<Y_edge(j+1);
                mask = mask_x & mask_y;   

                % Average values
                Zs(i_sect) = fun(Z(mask));            
                i_sect = i_sect + 1;
            end
        end                  
    case 'disk'
        [~, rho] = cart2pol(X, Y);        
        mask = rho <= Sectors.radius;        
        Zs = fun(Z(mask));        
    case 'ring'
        [~, rho] = cart2pol(X, Y);
        radius = Sectors.radius;
        radius(end) = radius(end) + 0.01; % to capture all points
        for i=1:Sectors.n_sect
            mask = rho >= radius(i) & rho < radius(i+1);        
            Zs(i) = fun(Z(mask));
        end
    case 'pie'
        theta_0 = Sectors.theta_0;
        n_angle = Sectors.n_angle;
        
        angle_edge = linspace(theta_0, theta_0 + 2*pi, n_angle + 1);
     
        [theta, ~] = cart2pol(X, Y);
        theta(theta<0) = theta(theta<0) + 2*pi;
        
        for i=1:n_angle
            mask = theta <= angle_edge(i) & theta < angle_edge(i+1);
            Zs(i) = fun(Z(mask));        
        end        
    case 'wedge'        
        % Get grid parameters (for readability)
        radius = Sectors.radius;
        theta_0 = Sectors.theta_0;
        n_angle = Sectors.n_angle;
        
        angle_edge = linspace(0, 2*pi, n_angle + 1);
               
        % Change to polar
        [theta, rho] = cart2pol(X, Y);
        
        % Rotate all angles to set to position initial theta as theta = 0
        theta = theta - theta_0;
        
        % Reconfigure to have only positive values
        theta(theta<0) = theta(theta<0) + 2*pi;
        
        % Central region
        mask = abs(rho) <= radius(1);
        Zs(1) = fun(Z(mask));

        % Loop through rings and angles
        i_sect = 2;
        for i=1:length(radius)-1     
            for t=1:Sectors.n_angle
                % Compute roi mask
                mask_x = abs(rho)>radius(i) & abs(rho)<=radius(i+1);
                mask_y = theta>=angle_edge(t) & theta<angle_edge(t+1);
                mask = mask_x & mask_y;   
                
                Zs(i_sect) = fun(Z(mask));

                i_sect = i_sect + 1;
            end
        end   
end
