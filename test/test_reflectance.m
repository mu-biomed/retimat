close all;clc;clearvars;
addpath(genpath('..'));

%% Reflectance map
file = '../data/raster.vol';
[header, seg, bscan, slo] = read_vol(file, 'coordinates');

% MR = reflectance_map(bscan, seg, 'mean', 'ILM', 'BM');

Tnfl = reflectance_map(bscan, seg, 'total', 'ILM', 'NFL_GCL');
Trpe = reflectance_map(bscan, seg, 'total', 'IDZ_RPE', 'BM');

att = reflectance_map(bscan, seg, 'NFL_RPE_att', header.scale_z, 'ILM', 'BM');

subplot(131);imagesc(Tnfl);
subplot(132);imagesc(Trpe);
subplot(133);imagesc(att);caxis([0 10]);
colormap(gray);