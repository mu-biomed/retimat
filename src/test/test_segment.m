close all;clc;clearvars;

addpath(genpath('..'));

file = '../data/raster.vol';
[header, ~, bscan,~] = read_vol(file);

%% Retina segmentation
I = bscan(:,:,1);

mask = seg_retina(I, header.scale_z, 100, 'otsu', true);

%% Layers segmentation

I = bscan(:,:,13);

seg = seg_layers(I, header.scale_z, true);
