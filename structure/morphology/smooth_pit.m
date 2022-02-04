function Z_smooth = smooth_pit(rho, Z, span)
%SMOOTH_PIT Smooth the foveal pit based on LOESS
%
%   Z_smooth = smooth_pit(rho, Z, span)
%   Returnes a smoothed version of Z based on LOESS smoothing.
%
%   Input arguments:
%  
%   'rho'            Matrix with radial coordinates of thickness values.          
%  
%   'Z'              Matrix with foveal pit values.           
%    
%   'span'           Percentage of samples included in each window.
%
%
%   Output arguments:
%  
%   'Z_smooth'       Smoothed version of the thickness map.
%
%   
%   Notes
%   -----
%   The function expects the data in radial format: [n_angle x n_points].
%   To use it for entire B-Scans convert data to radial format first.
%
%
%   References
%   ----------
%   [1] Romero-Bascones et al., Foveal Pit Morphology Characterization: A 
%   Quantitative Analysis of the Key Methodological Steps, Entropy, 2021
%   https://doi.org/10.3390/e23060699
%
%
%   Example
%   ---------      
%   % Example description
%
%   [header, seg]  = read_vol('my_file.vol','coordinates');
%   Thickness = compute_thickness(seg, 'TRT');
%   [X, Y, TRT] = resample_map(header.X_oct, header.Y_oct, Thickness.TRT, ...
%   'star', 'n_point', 100, 'max_d', 2.4, 'n_angle', n_angle);     
%   [theta, rho] = cart2pol(X, Y);
%   TRT_smooth = smooth_pit(rho, TRT, 20);     
%
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2021

if nargin ~= 3
    error(['3 input arguments are expected but ' nargin ' where provided']);
end

% Check the presence of nan values
if sum(isnan(Z(:))) > 0
    warning('NaN values in layer');
end

[n_angle, n_point] = size(Z);
n_bscan = n_angle/2;

Z_smooth = nan(n_angle, n_point);

for n=1:n_bscan

    % Reconstruct the B-Scan
    x = [-fliplr(rho(n+n_bscan,:)) rho(n,2:end)];
    z = [fliplr(Z(n+n_bscan,:)) Z(n,2:end)];  %  Dont get the center two times

    % Smooth the signal
    z_smooth = smooth(x, z, span, 'loess')';

    % From fitted B-Scan to radial again
    Z_smooth(n,:) = z_smooth((end-n_point+1):end);  %  Right part of the B-Scan
    Z_smooth(n+n_bscan,:) = z_smooth(n_point:-1:1);
end
