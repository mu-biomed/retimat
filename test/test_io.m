close all;clc;clearvars;
addpath(genpath('..'));

%% read raster vol
clc;close all;
file = '../data/raster.vol';
[h, seg, bscan, slo] = read_vol(file,'verbose', 'coordinates');

surf(h.X_slo, h.Y_slo, slo, 'EdgeColor', 'none');view(0,90);hold on;
colormap(gray);
scatter3(h.X_oct(:), h.Y_oct(:), repmat(max(slo(:))+1, 1, length(h.X_oct(:))),5,'b','filled');
scatter3(0, 0, max(slo(:))+1,'r','filled');
daspect([1 1 1]);

%% read star vol
file = '../data/star.vol';
% [h, seg, bscan, slo] = read_vol(file, 'visu', 'verbose', 'coordinates');
[h, seg, bscan, slo] = read_vol(file, 'verbose', 'coordinates');


surf(h.X_slo, h.Y_slo, slo, 'EdgeColor', 'none');view(0,90);hold on;
colormap(gray);
scatter3(h.X_oct(:), h.Y_oct(:), repmat(max(slo(:))+1, 1, length(h.X_oct(:))),5,'b','filled');
scatter3(0, 0, max(slo(:))+1,'r','filled');
daspect([1 1 1]);
%% read onh vol
clc;close all;
file = '../data/onh.vol';
[h, seg, bscan, slo] = read_vol(file, 'verbose', 'coordinates');

surf(h.X_slo, h.Y_slo, slo, 'EdgeColor', 'none');view(0,90);hold on;
colormap(gray);
scatter3(0, 0, max(slo(:))+1,'r','filled');
scatter3(h.X_oct, h.Y_oct, repmat(max(slo(:))+1, 1, length(h.X_oct)),5,'b','filled');
scatter3(h.X_oct(1), h.Y_oct(1), max(slo(:))+1,'g','filled');
daspect([1 1 1]);

%% read cube img
clc;close all;
file = '../data/Zeiss_Macula.img';
% file = 'C:/Users/dromero/Desktop/PNYU001E_Macular Cube 512x128_9-22-2017_12-17-35_OD_sn13716_cube_raw.img';
% file = 'C:/Users/dromero/Desktop/PNYU001E_HD 5 Line Raster_9-22-2017_12-22-6_OD_sn13722_lineEnhanced.img';
% file = 'C:/Users/dromero/Desktop/PNYU001E_Optic Disc Cube 200x200_11-19-2015_14-28-52_OS_sn8004_cube_raw.img';
[h, bscan] = read_img(file,[],true);

%% read bin img
clc;close all;
file = 'PNYU001E_Macular Cube 512x128_4-12-2018_16-18-22_OS_sn15113_lslo.bin';
I = read_bin(file);
imshow(I);

%% read fda
clc;close all;
file = '../data_private/test_3.fda';
[header, seg, bscan, fundus] = read_fda(file, true);

% ib = size(bscan,3)/2;
% imagesc(bscan(:,:,ib)); hold on;colormap(gray)
% boundaries = fields(seg);
% idx = 1:10;
% for i=idx
%     plot(size(bscan,1) - seg.(boundaries{i})(:,ib), 'LineWidth',2);
% end
% legend(boundaries(idx));

en_face = squeeze(mean(bscan,1)).';
trt = seg.RNFL_GCL - seg.IPL_INL;

subplot(131); imshow(fundus);
subplot(132); imagesc(en_face); colormap(gca, 'gray');
subplot(133); imagesc(trt);