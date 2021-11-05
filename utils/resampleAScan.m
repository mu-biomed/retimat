function [X1,Y1,Z1] = resampleAScan(X0,Y0,Z0,interpMethod,ResampleGrid, extrapolate)
% resampleAScan - resample A-Scan pattern to a different grid
%
% [X1,Y1,Z1] = resampleAScan(X0,Y0,Z0,interpMethod,ResampleGrid)
%
% Input arguments:
%   X0,Y0: matrixes with cartesian coordinates of A-Scan points
%   Z0: matrix of values (thickness)
%   interpMethod: interpolation method to be used
%   ResampleGrid: struct specifying the output grid (regular or star)
%   extrapolate: boolean 
%
% Output arguments:
%   X1,Y1: matrixes with cartesian coordinates of the resampling
%   grid
%   Z1: thickness values interpolated on the new grid
%
%  2021, Mondragon Unibertsitatea, Biomedical Engineering Department


%  Generate the resampling grid in both cartesian and polar coordinates
switch ResampleGrid.type
    
   case 'regular'
        step = ResampleGrid.step;
        maxD = ResampleGrid.maxD;
        [X1,Y1] = meshgrid(-maxD:step:maxD,-maxD:step:maxD);  
        
    case 'star'
        % Get params of grid
        n_angles = ResampleGrid.n_angles;
        step = ResampleGrid.step;
        maxD = ResampleGrid.maxD;
        
        % Define Grid
        Rho_grid = 0:step:maxD;
        Theta_grid = linspace(0,2*pi,n_angles+1);
        Theta_grid(end) = [];
        Rho_grid = repmat(Rho_grid,n_angles,1);
        Theta_grid = repmat(Theta_grid',1,size(Rho_grid,2));
        [X1,Y1] = pol2cart(Theta_grid,Rho_grid);   
end

% Interpolate over new grid
maskNum = ~isnan(Z0); % mask of not nan     
    
% Interpolation (griddata does not extrapolate)    
Z1 = reshape(griddata(X0(maskNum),Y0(maskNum),Z0(maskNum),...
             X1(:),Y1(:),interpMethod),size(X1));

if extrapolate
% EXTRAPOLATION for remaining NaNs (for outside sampled region)
% Caution not to extrapolate too far from the sampled region as
% the quality might be very poor
% scatteredInterpolant does not support cubic so it is linear           
    maskNum = ~isnan(Z1); % mask of not nan                
    if sum(~maskNum(:))>0
        warning('Extrapolation necessary');
        interpol = scatteredInterpolant(X1(maskNum),Y1(maskNum),Z1(maskNum),'linear');
        Z1 = reshape(interpol(X1(:),Y1(:)),size(X1));                 
    end
end