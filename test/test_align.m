close all;clc;clearvars;

addpath(genpath('../'));

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

% Visualization

