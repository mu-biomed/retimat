close all;clc;clearvars;

addpath(genpath('..'));
file = '../data/raster.vol';
[header, seg, bscan,~] = read_vol(file);

I = bscan(:,:,1);

mask = seg_retina(I, 'otsu_ascan', true);