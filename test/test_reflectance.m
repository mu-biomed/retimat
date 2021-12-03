close all;clc;clearvars;
addpath(genpath('..'));

%% Reflectance map
file = '../data/raster.vol';
[header, seg, bscan, slo] = read_vol(file, 'coordinates');

MR = reflectance_map(bscan, seg.ILM, seg.BM);
TR = reflectance_map(bscan, seg.ILM, seg.BM, 'total', header.scale_z);

subplot(121);imagesc(MR);
subplot(122);imagesc(TR);
colormap(gray);