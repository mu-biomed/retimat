close all;clc;clearvars;

% file = '/home/drombas/github/retimat/tutorials/data/example_raster.vol';

file = 'C:\Users\dromero\Desktop\GITHUB\retimat\tutorials\data\example_raster.vol';
addpath(genpath('C:\Users\dromero\Desktop\GITHUB\retimat'));

% addpath(genpath('/home/drombas/github/retimat/src'));
% [header, seg, bscan] = read_vol(file,'raw_pixel');
% 
% header.size_x = header.scale_x * (header.n_ascan - 1);
% header.size_y = header.scale_y * (header.n_bscan - 1);
% 
% segment_layers_aura(bscan, header);

%% AURA
[h, seg, bscan] = read_vol(file, 'raw_pixel');

params.resizedata = true;
params.minseg = false;
params.smooth = true;
params.segmethod = 1;

seg2 = segment_layers_aura(bscan, h, params);

figure;
subplot(121);
imagesc(bscan(:,:,12).^0.25);colormap(gray);
hold on;
boundaries = fields(seg);
for i=1:length(boundaries)
    plot(seg.(boundaries{i})(12,:), 'Linewidth', 2);
end
legend(boundaries, 'Interpreter','none');

subplot(122);
imagesc(bscan(:,:,12).^0.25); colormap(gray);
hold on;
boundaries = fields(seg2);
for i=1:length(boundaries)
    plot(seg2.(boundaries{i})(12,:), 'Linewidth', 2);
end
legend(boundaries, 'Interpreter','none');