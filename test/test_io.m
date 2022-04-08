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

%% Read TABS output
clc;close all;

folder = '../data/1000084_21011_0_0';
[h, seg, bscan, fundus] = read_tabs(folder, true);

% Plotting: fundus
subplot(141);
imagesc(fundus);
hold on;
x = [h.oct_fundus_ul(1) h.oct_fundus_br(1)];
y = [h.oct_fundus_ul(2) h.oct_fundus_br(2)];

scatter(x,y,'g');
plot([x(1) x(2) x(2) x(1) x(1)],[y(1) y(1) y(2) y(2) y(1)],'--g');
axis off;

% Pixel to um
TRT = double(seg.BM - seg.ILM) * h.scale_z * 1e3;

% En face
subplot(142);
en_face = squeeze(mean(bscan,2));
imagesc(en_face);
colormap(gca,gray);
hold on;
scatter(h.fovea_center_px(2), h.fovea_center_px(1),15,'r','filled');
axis off;
title('En face');

% B-scan
subplot(143);
i_bscan = round(h.n_bscan/2);
imagesc(squeeze(bscan(i_bscan,:,:)));hold on;
for i_layer=1:h.n_layer
    layer = h.layers{i_layer};
    plot(seg.(layer)(i_bscan, :));
end
colormap(gca,gray);
axis off;
title(['Bscan:' num2str(i_bscan)]);

% TRT
subplot(144);
surf(h.X_oct,h.Y_oct,TRT,'EdgeColor','none');view(0,90);
hold on;
scatter3(h.x_fovea, h.y_fovea, max(TRT(:)), 15, 'r', 'filled');
axis([-3 3 -3 3]);
title('TRT');
colorbar;