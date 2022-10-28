function [X_radial, Y_radial, Z_radial] = rearange_star_coords(X, Y, Z)
% REARANGE_STAR_CORDS - Organize star coords radially
%
%   [X_rad, Y_radial, Z_radial] = rearange_star_coords(X, Y, Z)
%
%   This functions changes the structure of data organized as matrix of n_bscan
%   x n_ascan into a matrix of n_angles x n_ascan/2, where n_angles is twice
%   the number of bscans (n_bscan).
%
%   Input arguments:
%  
%   'X'              Matrix of n_bscan x n_ascan with A-Scan X coordinates.
%  
%   'Y'              Matrix of n_bscan x n_ascan with A-Scan Y coordinates.
%  
%   'Z   '           Matrix of n_bscan x n_ascan with numeric values
%                    (segmentation, thickness). Alternatively, a struct with
%                    each field containing a different matrix (useful when 
%                    several layer data is stored)
%
%  
%   Output arguments:
%  
%   'X_radial'       Matrix of 2*n_bscan x n_ascan/2 with A-Scan X coordinates 
%                    organized radially.        
%  
%   'Y_radial'       Matrix of 2*n_bscan x n_ascan/2 with A-Scan Y coordinates 
%                    organized radially.  
%
%   'Z_radial'       Matrix of 2*n_bscan x n_ascan/2 with numeric values
%                    (segmentation, thickness) organized radially. 
%                    Alternatively, it can be a struct of matrixes if that is 
%                    what was provided as input.
%
%
%   Notes
%   -----
%   This has been only tested with Spectralis star scans containing 12 bscans
%   and might not work perfectly with other configurations.
%
%
%   References
%   ----------
%
%   Example
%   ---------      
%   % Read a star oct file and convert the coordinates
%     
%     [header, segment]] = read_vol('myoct.vol')
%     [X_radial, Y_radial, Data] = rearange_star_coords(header.X_oct, header.Y_oct, segment)
%     
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2021


[n_bscan, n_ascan] = size(X);

n_angle = 2*n_bscan;

if mod(n_ascan,2) ~= 0
     warning("Number of A-Scans is odd. Values will be repeated");
end
n_ascan_rad = ceil(n_ascan/2); % number of A-scans in every direction

% Remap B-Scans to a radial pattern
X_radial = nan(n_angle, n_ascan_rad);
Y_radial = nan(n_angle, n_ascan_rad);     
        
% Right part of each B-Scan
X_radial(1:n_bscan, :) = X(:, (end - n_ascan_rad + 1):end);
Y_radial(1:n_bscan, :) = Y(:, (end - n_ascan_rad + 1):end);
        
% Left part of each B-Scan (needs to be inverted)
X_radial((n_bscan + 1):end, :) = X(:, n_ascan_rad:-1:1);
Y_radial((n_bscan + 1):end, :) = Y(:, n_ascan_rad:-1:1);  

% Remap also layer data to a radial pattern
if isnumeric(Z)
    Z_radial = nan(n_angle, n_ascan_rad);

    Z_radial(1:n_bscan, :) = Z(:, (end - n_ascan_rad + 1):end);
    Z_radial((n_bscan+1):end, :) = Z(:, n_ascan_rad:-1:1);    
    
elseif istruct(Z)
    data_fields = fields(Z);
    for i=1:length(data_fields)
        field = data_fields{i};
        
        Z_radial.(field) = nan(n_angle, n_ascan_rad);
        
        Z_radial.(field)(1:n_bscan, :) = Z.(field)(:, (end - n_ascan_rad + 1):end);
        Z_radial.(field)((n_bscan+1):end, :) = Z.(field)(:, n_ascan_rad:-1:1);
    end
    
else
    error("Unsupported data type for Z");
end
end
