function [tform,fig] = register_fundus(raster, onh, transform, return_fig, visible)
%REGISTER_FUNDUS Register 2 fundus images
%
%   Usage example OUT = template(IN1)
%   Detail explanation goes here
%
%   Input arguments:
%  
%   'ARG1'           Description of the argument. Type and purpose.          
%  
%                    Accepted values
%
%                    Default: 
%            
%  
%  
%   Output arguments:
%  
%   'ARG1'           Description of the argument. Type and purpose.          
%  
%
%   
%   Notes
%   -----
%   Important usage informationAnother name for a gray-level co-occurrence matrix is a gray-level
%   spatial dependence matrix.
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
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2022

% Get data
X_fix = raster.X_fun;
Y_fix = raster.Y_fun;
Z_fix = double(raster.fundus)./255;

X_mov = onh.X_fun;
Y_mov = onh.Y_fun;
Z_mov = double(onh.fundus)./255;

% Equalize histogram
Z_fix = adapthisteq(Z_fix);
Z_mov = adapthisteq(Z_mov);

% Detect SURF Features
fixed_points = detectSURFFeatures(Z_fix,'MetricThreshold',750,...
                                        'NumOctaves',3,...
                                        'NumScaleLevels',5);

moving_points = detectSURFFeatures(Z_mov,'MetricThreshold',750,...
                                         'NumOctaves',3,...
                                         'NumScaleLevels',5);

% Extract features
[fixed_features, fixed_valid_points] = extractFeatures(Z_fix, fixed_points,...
                                                       'Upright',true);
                                                   
[moving_features, moving_valid_points] = extractFeatures(Z_mov,moving_points,...
                                                         'Upright',true);

% Match features
index_pairs = matchFeatures(fixed_features, moving_features, ...
                            'MatchThreshold',1,...
                            'MaxRatio',0.50);

fixed_matched_points = fixed_valid_points(index_pairs(:,1));
moving_matched_points = moving_valid_points(index_pairs(:,2));


% Points pixel to coordinates (raster)
fixed_match_real = nan(length(fixed_matched_points),2);

y_pixel1 = 1;
y_real1 = Y_fix(1,1);
y_pixel2 = 768;
y_real2 = Y_fix(end,1);
ay = (y_real2 - y_real1)/(y_pixel2 - y_pixel1);
by = y_real1 - ay*y_pixel1;

x_pixel1 = 1;
x_real1 = X_fix(1,1);
x_pixel2 = 768;
x_real2 = X_fix(1,end);
ax = (x_real2 - x_real1)/(x_pixel2 - x_pixel1);
bx = x_real1 - ax*x_pixel1;

fixed_match_real(:,1) = ax*fixed_matched_points.Location(:,1) + bx;
fixed_match_real(:,2) = ay*fixed_matched_points.Location(:,2) + by;

% Points pixel to coordinates (raster)
moving_match_real = nan(length(moving_matched_points),2);

y_pixel1 = 1;
y_real1 = Y_mov(1,1);
y_pixel2 = 768;
y_real2 = Y_mov(end,1);
ay = (y_real2 - y_real1)/(y_pixel2 - y_pixel1);
by = y_real1 - ay*y_pixel1;

x_pixel1 = 1;
x_real1 = X_mov(1,1);
x_pixel2 = 768;
x_real2 = X_mov(1,end);
ax = (x_real2 - x_real1)/(x_pixel2 - x_pixel1);
bx = x_real1 - ax*x_pixel1;

moving_match_real(:,1) = ax*moving_matched_points.Location(:,1) + bx;
moving_match_real(:,2) = ay*moving_matched_points.Location(:,2) + by;

% Estimate transformation
tform = estimateGeometricTransform2D(moving_match_real,fixed_match_real,transform);

if return_fig 

    [X_aligned,Y_aligned] = transformPointsForward(tform,X_mov,Y_mov);

    close all;    
    fig = figure('position',[0 0 1300 600],'visible',visible);    
    
    subplot(131); hold on;
    surf(X_fix,Y_fix,Z_fix,'EdgeColor','none');title('raster (fixed)');view(0,90);colormap(gray);daspect([1 1 1]);
    xlim([-5 5]);ylim([-5 5]);
    scatter3(fixed_match_real(:,1),fixed_match_real(:,2),repmat(max(Z_fix(:)),size(fixed_match_real,1),1),'r','filled');

    subplot(132); hold on;
    surf(X_mov,Y_mov,Z_mov,'EdgeColor','none');title('onh (moving)');view(0,90);colormap(gray);daspect([1 1 1]);
    scatter3(moving_match_real(:,1),moving_match_real(:,2),repmat(max(Z_mov(:)),size(fixed_match_real,1),1),'g','filled');

    % Plot both figures together
    x_right_1 = max(X_fix(:));
    x_right_2 = max(X_aligned(:));
    
    [Xr,Yr] = meshgrid(linspace(x_right_1,x_right_2,300),linspace(min(Y_aligned(:)),max(Y_aligned(:)),300));
    step = 4; % to save time
    x = X_aligned(:);
    y = Y_aligned(:);
    z = Z_mov(:);
    interpol = scatteredInterpolant(x(1:step:end),y(1:step:end),z(1:step:end));
    Z_aligned_right = reshape(interpol(Xr(:),Yr(:)),size(Xr));
    
    subplot(133);hold on;
    surf(X_fix,Y_fix,Z_fix,'EdgeColor','none');
    surf(Xr,Yr,Z_aligned_right,'EdgeColor','none');
    view(0,90);colormap(gray);daspect([1 1 1]);caxis([min(Z_fix(:)) max(Z_fix(:))]);
    xlim([-5 x_right_2]);ylim([-5 5]);
    
else
    fig = nan;
end
