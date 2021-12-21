close all;clc;clearvars;
addpath(genpath('..'));

%%
file = '../data/raster.vol';
[header, seg, bscan] = read_vol(file);
Thickness = compute_thickness(seg, 'all', header.scale_z);

figure;
i_bscan = 2;
imshow(bscan(:,:,i_bscan));hold on;
layers = fields(seg);
for f=1:length(layers)
    plot(seg.(layers{f})(i_bscan,:), 'Linewidth', 1);
end
legend(layers, 'Interpreter', 'none')

figure;
fname = fields(Thickness);
for f=1:length(fname)
    subplot(3,5,f);
    imagesc(Thickness.(fname{f}));
    title(fname{f});
end