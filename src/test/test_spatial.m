close all;clc;clearvars;

addpath(genpath('../'));

%% Test find fovea
file = '../../tutorials/data/example_raster.vol';
[header, seg, ~, ~] = read_vol(file, 'get_coordinates');
Thickness = compute_thickness(seg, 'TRT', header.scale_z);
TRT = Thickness.TRT;

X = header.X_oct;
Y = header.Y_oct;

methods = {'min','resample_min','smooth_min','flood'};
for m=1:length(methods)
    [x_fovea, y_fovea] = find_fovea(X, Y, TRT, methods{m});
    
    subplot(1,length(methods),m);
    surf(X, Y, TRT, 'EdgeColor', 'none');view(0,90);
    hold on;
    scatter3(x_fovea, y_fovea, max(TRT(:)), 'r', 'filled');
    title(methods{m});
end
%% Test sectorization

% Read vol
file = '../data/raster.vol';
[header, seg, ~, ~] = read_vol(file,'verbose', 'coordinates');

% Compute thickness
Thickness = compute_thickness(seg, 'TRT', header.scale_z);

% Resampling
[X, Y, TRT] = resample_map(header.X_oct, header.Y_oct, Thickness.TRT, ...
    'regular', 'n_point', 100, 'max_d', 2.5);
% surf(X,Y,TRT,'EdgeColor','none'); view(0,90);

% Sectorization
[Z, G] = sectorize_map(X, Y, TRT, 'mean', 'ring', [1.5 2 2.5 3]);
[Z, G] = sectorize_map(X, Y, TRT, 'mean', 'etdrs');

% Visualization

