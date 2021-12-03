close all;clc;clearvars;
addpath(genpath('..'));

%% read raster vol
clc;close all;
file = '../data/raster.vol';
[header, seg, bscan, slo] = read_vol(file,'verbose', 'coordinates');

surf(header.X_slo, header.Y_slo, slo, 'EdgeColor', 'none');view(0,90);hold on;
colormap(gray);
scatter3(header.X_oct(:), header.Y_oct(:), repmat(max(slo(:))+1, 1, length(header.X_oct(:))),5,'b','filled');
scatter3(0, 0, max(slo(:))+1,'r','filled');
daspect([1 1 1]);

%% read star vol
file = '../data/star.vol';
% [header, seg, bscan, slo] = read_vol(file, 'visu', 'verbose', 'coordinates');
[header, seg, bscan, slo] = read_vol(file, 'verbose', 'coordinates');


surf(header.X_slo, header.Y_slo, slo, 'EdgeColor', 'none');view(0,90);hold on;
colormap(gray);
scatter3(header.X_oct(:), header.Y_oct(:), repmat(max(slo(:))+1, 1, length(header.X_oct(:))),5,'b','filled');
scatter3(0, 0, max(slo(:))+1,'r','filled');
daspect([1 1 1]);
%% read onh vol
clc;close all;
file = '../data/onh.vol';
[header, seg, bscan, slo] = read_vol(file, 'verbose', 'coordinates');

surf(header.X_slo, header.Y_slo, slo, 'EdgeColor', 'none');view(0,90);hold on;
colormap(gray);
scatter3(0, 0, max(slo(:))+1,'r','filled');
scatter3(header.X_oct, header.Y_oct, repmat(max(slo(:))+1, 1, length(header.X_oct)),5,'b','filled');
scatter3(header.X_oct(1), header.Y_oct(1), max(slo(:))+1,'g','filled');
daspect([1 1 1]);