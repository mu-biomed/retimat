close all;clc;clearvars;

addpath(genpath('..'));

file = '../data/raster.vol';
[header, ~, bscan,~] = read_vol(file);

%% Retina segmentation
I = bscan(:,:,1);

mask = seg_retina(I, header.scale_z, 100, 'otsu', true);

%% Layers segmentation
% file = '../data/raster.vol';
% [header, ~, bscan,~] = read_vol(file);

file = '../data/Zeiss_Macula.img';
[bscan, header] = read_img(file);

I = bscan(:,:,1);
layers = {'ilm','isos','elm','bm'};
seg = seg_layers(I, header.scale_z, layers, true);

% figure;
% imagesc(I);hold on;
% for i=1:length(layers)
%     plot(seg.(layers{i}));
% end
% colormap(gray);
% legend(layers);
