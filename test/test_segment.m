close all;clc;clearvars;

addpath(genpath('..'));

% [header, ~, bscan,~] = read_vol(file);

%% Retina segmentation
I = bscan(:,:,1);

mask = seg_retina(I, header.scale_z, 150, 'otsu', true);

%% Layers segmentation
file = '../data/raster.vol';
[header, ~, bscan,~] = read_vol(file);

file = '../data/Zeiss_Macula.img';
[bscan, header] = read_img(file);

for ib=1:10

    I = bscan(:,:,ib);
    I = imresize(I, [400 512]);

    seg = seg_layers(I, [], [], false);
    
    subplot(2,5,ib);
    imagesc(I); hold on; axis off;
    plot(seg.ilm,'r');
    plot(seg.isos,'b');
    plot(seg.bm,'g');
    plot(seg.nfl,'y');
end
colormap(gray);

% layers = {'isos','ilm'};
% seg = seg_layers(I, header.scale_z, layers, false, true);

% figure;
% imagesc(I);hold on;
% for i=1:length(layers)
%     plot(seg.(layers{i}));
% end
% colormap(gray);
% legend(layers);
