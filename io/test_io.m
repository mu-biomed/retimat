close all;clc;clearvars;

%% read raster vol
clc;close all;
file = '/home/david/Desktop/biocruces-data/data/source/03_OCT_VIRTUAL_GAMES/VIRTUAL_GAMES_year_0_OCT/VOL_macula_raster/cVG06_V_3640.vol';
[header, seg, bscan, slo] = read_vol(file, 'visu','verbose', 'coordinates');

surf(header.X_slo, header.Y_slo, slo, 'EdgeColor', 'none');view(0,90);hold on;
colormap(gray);
scatter3(header.X_oct(:), header.Y_oct(:), repmat(max(slo(:))+1, 1, length(header.X_oct(:))),5,'b','filled');
scatter3(0, 0, max(slo(:))+1,'r','filled');
daspect([1 1 1]);

%% read star vol
file = '/home/david/Desktop/biocruces-data/data/source/03_OCT_VIRTUAL_GAMES/VIRTUAL_GAMES_year_0_OCT/VOL_macula_star/cVG07_c_2167.vol';
[header, seg, bscan, slo] = read_vol(file, 'visu', 'verbose', 'coordinates');


%% read onh vol
clc;close all;
file = '/home/david/Desktop/biocruces-data/data/source/03_OCT_VIRTUAL_GAMES/VIRTUAL_GAMES_year_0_OCT/VOL_onh/cVG18_c_4666_2.vol';
[header, seg, bscan, slo] = read_vol(file, 'visu','verbose', 'coordinates');

surf(header.X_slo, header.Y_slo, slo, 'EdgeColor', 'none');view(0,90);hold on;
colormap(gray);
scatter3(header.x_onh, header.y_onh, max(slo(:))+1,'r','filled');
scatter3(header.X_oct, header.Y_oct, repmat(max(slo(:))+1, 1, length(header.X_oct)),5,'b','filled');
scatter3(header.X_oct(1), header.Y_oct(1), max(slo(:))+1,'g','filled');
daspect([1 1 1]);