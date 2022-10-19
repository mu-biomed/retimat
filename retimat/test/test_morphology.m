close all;clc;clearvars;
addpath(genpath('..'));

%% Pit models
% Read file and compute
file = '../data/raster.vol';
[header, seg] = read_vol(file, 'coordinates');
Thickness = compute_thickness(seg, {'TRT', 'GCIP'}, header.scale_z);

% Resampling
[X, Y, TRT] = resample_map(header.X_oct, header.Y_oct, Thickness.TRT, ...
    'star', 'n_angle', 24, 'max_d', 2.5, 'n_point', 100);
[~, ~, GCIP] = resample_map(header.X_oct, header.Y_oct, Thickness.GCIP, ...
    'star', 'n_angle', 24, 'max_d', 2.5, 'n_point', 100);

[theta, rho] = cart2pol(X, Y);

model = 'Yadav';
[TRT_fit, fc_TRT] = fit_pit_model(theta, rho, TRT, model);
[GCIP_fit, fc_GCIP] = fit_pit_model(theta, rho, GCIP, model);

subplot(121);hold on;
plot(TRT', 'k');plot(TRT_fit', 'r');

subplot(122);hold on;
plot(GCIP', 'k');plot(GCIP_fit', 'r');
title('GCIPL');

%% Morphological parameters
clc;
file = '../data/raster.vol';
[header, seg] = read_vol(file, 'coordinates');
Thickness = compute_thickness(seg, {'TRT', 'GCIP'}, header.scale_z);
[X, Y, TRT] = resample_map(header.X_oct, header.Y_oct, Thickness.TRT, ...
    'star', 'n_angle', 24, 'max_d', 2.5, 'n_point', 100);
[theta, rho] = cart2pol(X, Y);

% [TRT_fit, fc_TRT] = fit_pit_model(theta, rho, TRT, 'Scheibe');

X = get_morph_params(rho, TRT, 'all', true);
% X = get_morph_params(rho, TRT, {'cft', 'max_slope'}, true);

%% Pit smoothing
clc;close all;
file = '../data/raster.vol';
[header, seg] = read_vol(file, 'coordinates');
Thickness = compute_thickness(seg, {'TRT', 'GCIP'}, header.scale_z);
[X, Y, Z] = resample_map(header.X_oct, header.Y_oct, Thickness.GCIP, ...
    'star', 'n_angle', 24, 'max_d', 2.5, 'n_point', 100);
[theta, rho] = cart2pol(X, Y);

Z_smooth = smooth_pit(rho, Z, 0.2);

plot(Z','b');hold on;plot(Z_smooth','r')
