function Paramorf = get_morph_params(rho, Y)
% GET_MORPH_PARAMS - compute morphological parameters of the foveal pit
% based on TRT values
%
% Paramorf = getMorphParams(rho,TRT)
%
% Input arguments:
%   rho: matrix with rho coordinates
%   TRT: matrix with TRT points (each row is an angular direction)
%
% Output arguments:
%   Paramorf: struct with morphological parameters
%
%  2021, Mondragon Unibertsitatea, Biomedical Engineering Department

% Get step
n_angles = size(Y, 1);

step = rho(1,2) - rho(1,1);

% Check both axis in mm
if max(Y(1))>1
   Y = 1e-3*Y;
end

% Radial parameter computation
for n=1:n_angles

    x = rho(n,:);
    y = Y(n,:);

    % Central Foveal Thickness
    Paramorf.cft(n) = 1e3*y(1);

    % Rim height 
    [Paramorf.rim_height(n),ind_max] = max(1e3*y);
    ind_max = ind_max(1); % por si > 1 mínimo

    % Rim Radius
    Paramorf.rim_radius(n) = x(ind_max);

    % Max Slope
    slope = diff(y(1:ind_max))./step;
    [max_slope, ind_max_slope] = max(slope);
    max_slope = max_slope(1);    
    Paramorf.max_slope(n) = 180*atan(max_slope)/pi;
    
    % Max slope radius
    Paramorf.max_slope_radius(n) = x(ind_max_slope);
    
    % Max slope height
    Paramorf.max_slope_height(n) = y(ind_max_slope);
    
    % Mean slope
    Paramorf.mean_slope(n) = 180*atan(mean(slope))/pi;
end