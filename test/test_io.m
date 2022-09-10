close all;clc;clearvars;
addpath(genpath('..'));

visu = false;

% Launch test suite
test_vol_raster(visu);
test_vol_star(visu);
test_vol_onh(visu);
test_img_cube(visu);
test_bin_cube(visu);
test_e2e(visu);
test_fda(visu);


function test_vol_raster(visu)
file = '../data/raster.vol';
[h, seg, bscan, slo] = read_vol(file, 'verbose', 'coordinates');

if visu
    surf(h.X_fun, h.Y_fun, slo, 'EdgeColor', 'none');view(0,90);hold on;
    colormap(gray);
    scatter3(h.X_oct(:), h.Y_oct(:), repmat(max(slo(:))+1, 1, length(h.X_oct(:))),5,'b','filled');
    scatter3(0, 0, max(slo(:))+1,'r','filled');
    daspect([1 1 1]);
end
end

function test_vol_star(visu)
file = '../data/star.vol';
[h, seg, bscan, slo] = read_vol(file, 'verbose', 'coordinates');

if visu
    surf(h.X_fun, h.Y_fun, slo, 'EdgeColor', 'none');view(0,90);hold on;
    colormap(gray);
    scatter3(h.X_oct(:), h.Y_oct(:), repmat(max(slo(:))+1, 1, length(h.X_oct(:))),5,'b','filled');
    scatter3(0, 0, max(slo(:))+1,'r','filled');
    daspect([1 1 1]);
end

end

function test_vol_onh(visu)
file = '../data/onh.vol';
[h, seg, bscan, slo] = read_vol(file, 'verbose', 'coordinates');

if visu
    surf(h.X_fun, h.Y_fun, slo, 'EdgeColor', 'none');view(0,90);hold on;
    colormap(gray);
    scatter3(0, 0, max(slo(:))+1,'r','filled');
    scatter3(h.X_oct, h.Y_oct, repmat(max(slo(:))+1, 1, length(h.X_oct)),5,'b','filled');
    scatter3(h.X_oct(1), h.Y_oct(1), max(slo(:))+1,'g','filled');
    daspect([1 1 1]);
end
end

function test_img_cube(visu)
file = '../data/Zeiss_Macula.img';
% file = 'C:/Users/dromero/Desktop/PNYU001E_Macular Cube 512x128_9-22-2017_12-17-35_OD_sn13716_cube_raw.img';
% file = 'C:/Users/dromero/Desktop/PNYU001E_HD 5 Line Raster_9-22-2017_12-22-6_OD_sn13722_lineEnhanced.img';
% file = 'C:/Users/dromero/Desktop/PNYU001E_Optic Disc Cube 200x200_11-19-2015_14-28-52_OS_sn8004_cube_raw.img';
[h, bscan] = read_img(file,[],true);

end

function test_bin_cube(visu)
file = '../data_private/PNYU001E_Macular Cube 512x128_4-12-2018_16-18-22_OS_sn15113_lslo.bin';
I = read_bin(file);
if visu
    imshow(I);
end
end

function test_e2e(visu)
file = '../data_private/oct_1.e2e';
[header, seg, bscan, fundus] = read_e2e(file, 'verbose');
end

function test_fda(visu)
file = '../data_private/test_0.fda';
[header, seg, bscan, fundus] = read_fda(file, 'verbose', 'coordinates');

en_face = squeeze(mean(bscan,1)).';
trt = seg.ILM - seg.BM;

subplot(141); imshow(fundus);
subplot(142); imagesc(en_face); colormap(gca, 'gray');
subplot(143); imagesc(trt);
subplot(144); surf(header.X, header.Y, trt,'EdgeColor','none');view(0,90);
end