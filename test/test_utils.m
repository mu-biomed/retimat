close all;clc;clearvars;
addpath(genpath('../'));

%% Rearange star coordinates
[header, seg] = read_vol('../data/star.vol','coordinates');
X = header.X_oct;
Y = header.Y_oct;

TRT = header.scale_z * abs(seg.ILM - seg.BM);
subplot(121);plot(TRT');

[X_radial, Y_radial, TRT_radial] = rearange_star_coords(X, Y, TRT);
[theta, rho] = cart2pol(X_radial, Y_radial);
subplot(122);plot(rho',TRT_radial');

%% GLCM features
close all;clc;clearvars;

% I = randn(200, 200);
I = imread('cameraman.tif');
GLCM = graycomatrix(I,'NumLevels',250);
X = GLCM_features(GLCM);

%% Lacunarity
close all;clc;clearvars;
I = imread('cameraman.tif');