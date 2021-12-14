close all;clc;clearvars;
addpath(genpath('..'));
%% Pit models
% Read file and compute
file = '../data/raster.vol';
[header, seg, ~, ~] = read_vol(file, 'coordinates');
Thickness = compute_thickness(seg, 'TRT', header.scale_z);

% Resampling
[X, Y, TRT] = resample_map(header.X_oct, header.Y_oct, Thickness.TRT, ...
    'star', 'n_angle', 24, 'max_d', 2.5, 'n_point', 100);
[theta, rho] = cart2pol(X, Y);

[Z, fit_c] = fit_pit_model(theta, rho, TRT, 'Breher');

plot(Z');