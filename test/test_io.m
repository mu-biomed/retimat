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

%% read cube img
clc;close all;
% file = '/home/david/Desktop/PNYU061E_Macular Cube 512x128_1-17-2018_10-27-10_OD_sn14547_cube_raw.img';
file = 'C:/Users/dromero/Desktop/PNYU001E_Macular Cube 512x128_9-22-2017_12-17-35_OD_sn13716_cube_raw.img';
% file = 'C:/Users/dromero/Desktop/PNYU001E_HD 5 Line Raster_9-22-2017_12-22-6_OD_sn13722_lineEnhanced.img';
file = 'C:/Users/dromero/Desktop/PNYU001E_Optic Disc Cube 200x200_11-19-2015_14-28-52_OS_sn8004_cube_raw.img';
[bscan, header] = read_img(file,[],true);