function [Zs, Sectors] = sectorize_map(X, Y, Z, metric, Sectors, varargin)
%SECTORIZE_MAP Sectorize a 2D map into several sectors
%
%   Zs = sectorize_map(X, Y, Z, Sectors)
%   Detail explanation goes here
%
%   Input arguments:
%  
%   'X'              X coordinates of map points
%
%   'Y'              Y coordinates of map points
%
%   'Z'              Z coordinates of map points
%
%   'metric'         metric to be used to sectorize data. 
%                    Options: ['mean', 'std', skewness', 'kurtosis']
%
%   'Sectors'        Two options:
%                    - String defining the sectorization type.  It must be 
%                    followed by appropriate arguments.
%                    - Struct with sectorization info created beforehand.
%
%   'varargin'       Input arguments necessary when Sectors is a string. See 
%                    examples to see which arguments are required for each type
%                    of sectorization.
%            
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
%   
%
%
%
%   Example 1
%   ---------      
%   % Example description
%
%     I = [1 1 5 6 8 8;2 3 5 7 0 2; 0 2 3 5 6 7];
%     [Zs, Sectors] = sectorize_map(X, Y, Z, 'mean', 'ring', [1.5 2 2.5 3])
%     
%
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2021


if nargin < 5
    error("A minimum of 4 input arguments must be provided");
end

if ischar(Sectors)
    aux.type = Sectors;
    Sectors = aux;
    
    switch Sectors.type
        case 'regular'
            Sectors.n_x = varargin{1};
            Sectors.n_y = varargin{2};
            Sectors.n_sect = nx*ny;
            Sectors.X_edge = linspace(min(X(:)), max(X(:)), Sectors.n_x + 1);
            Sectors.Y_edge = linspace(min(Y(:)), max(Y(:)), Sectors.n_y + 1);        
        case 'disk'
            Sectors.radius = varargin{1};        
            Sectors.n_sect = 1;
        case 'ring'
            Sectors.radius = varargin{1};
            Sectors.n_sect = length(Sectors.radius)-1;
        case 'pie'
            Sectors.n_angle = varargin{1};
            Sectors.theta_0 = varargin{2};
            Sectors.n_sect = Sectors.n_angle;
            Sectors.theta = linspace(0, 2*pi, Sectors.n_angle + 1);
        case 'wedge'
            Sectors.radius = varargin{1};
            Sectors.n_angle = varargin{2};            
            Sectors.theta_0 = varargin{3};
            Sectors.n_sect = 1 + (length(Sectors.radius)-1)*Sectors.n_angle;
        case 'etdrs'
            Sectors.radius = [0.5 1.5 3];
            Sectors.angle = [0 pi/2 pi 3*pi/2];
            Sectors.n_sect = 9;
    end
    
elseif ~isstruct(Sectors)
    error("Sectors must be either a struct or a string");
end

Zs = nan(1, Sectors.n_sect);

switch metric
    case 'mean'
        fun = @(x) nanmean(x);
    case 'median'
        fun = @(x) median(x);
    case 'std'
        fun = @(x) nanstd(x);
    case 'skewness'
        fun = @(x) skewness(x);
    case 'kurtosis'
        fun = @(x) kurtosis(x);
    case 'snr'
        fun = @(x) mean(x)/std(x);
    otherwise
        error("Unknown metric. Accepted options=['mean','std','skewness','kurtosis','snr']");
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
    case 'disc'
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
        [theta, ~] = cart2pol(X, Y);
        theta = theta - Sectors.theta_0;
        theta(theta<0) = theta(theta<0) + 2*pi;
        
        for i=1:n_angle
            mask = theta <= Sectors.angle(i);
            Zs(i) = fun(Z(mask));        
        end        
    case 'wedge'        
        % Get grid parameters (for readability)
        radius = Sectors.radius;
        theta_0 = Sectors.theta_0;
        theta = Sectors.theta;
               
        % Change to polar
        [theta_s,rho_s] = cart2pol(X, Y);
        
        % Rotate all angles to set to position initial theta as theta = 0
        theta_s = theta_s - theta_0;
        
        % Reconfigure to have only positive values
        theta_s(theta_s<0) = theta_s(theta_s<0) + 2*pi;
        
        % Central region
        mask = abs(rho_s) <= radius(1);
        Z(1) = fun(Z(mask));

        % Loop through rings and angles
        i_sect = 2;
        for i = 1:length(radius)-1     
            for t=1:Sectors.n_angle
                % Compute roi mask
                mask_x = abs(rho_s)>radius(i) & abs(rho_s)<=radius(i+1);
                mask_y = theta_s >= theta(t) & theta_s<theta(t+1);
                mask = mask_x & mask_y;   
                
                Z(i_sect) = fun(Z(mask));

                i_sect = i_sect+1;
            end
        end   
end
