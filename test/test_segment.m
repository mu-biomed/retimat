close all;clc;clearvars;

addpath(genpath('..'));
file = '../data/raster.vol';
[header, ~, bscan,~] = read_vol(file);

I = bscan(:,:,1);

mask = seg_retina(I, header.scale_z, 100, 'otsu', true);