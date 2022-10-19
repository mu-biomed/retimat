close all;clc;clearvars;
addpath(genpath('..'));

visu = true;

% Launch test suite
test_vol_raster(visu);
% test_vol_star(visu);
% test_vol_onh(visu);
% test_img_cube(visu);
% test_bin_cube(visu);
test_e2e(visu);
test_fda(visu);
test_iowa_star(visu);
test_iowa_raster(visu);

function test_vol_raster(visu)
file = '../data/raster.vol';
[h, seg, bscan, slo] = read_vol(file, 'verbose', 'get_coordinates');

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
[h, seg, bscan, slo] = read_vol(file, 'verbose', 'get_coordinates');

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
[h, seg, bscan, slo] = read_vol(file, 'verbose', 'get_coordinates');

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
[h, bscan] = read_img(file,[]);

end

function test_bin_cube(visu)
file = '../data_private/PNYU001E_Macular Cube 512x128_4-12-2018_16-18-22_OS_sn15113_lslo.bin';
I = read_bin(file);
if visu
    imshow(I);
end
end

function test_e2e(visu)
file = '../data_private/e2e/oct_1.e2e';
[header, seg, bscan, fundus] = read_e2e(file, 'verbose');
end

function test_fda(visu)
file = '../data_private/fda/test_0.fda';
[header, seg, bscan, fundus] = read_fda(file, 'verbose', 'get_coordinates');

en_face = squeeze(mean(bscan,1)).';
trt = seg.ILM - seg.BM;

subplot(141); imshow(fundus);
subplot(142); imagesc(en_face); colormap(gca, 'gray');
subplot(143); imagesc(trt);
subplot(144); surf(header.X, header.Y, trt,'EdgeColor','none');view(0,90);
end

function test_iowa_star(visu)
%% Star
in_dir = '../data_private/iowa/vol_star';
images = {'HC001AM_H_2341', 'HC001AM_H_2342',...
          'HC002AM_H_2347', 'HC002AM_H_2348',...
          'HC002AF_A_2357', 'HC002AF_A_2358'};
n_image = length(images);
      
layers = {'TRT','RNFL','GCIP','INL'};
n_layer = length(layers);

for i_image=1:n_image
    image = images{i_image};
    
    file_vol  = [in_dir '/' image '/' image '_OCT_Iowa.vol'];
    file_iowa = [in_dir '/' image '/' image '_Surfaces_Iowa.xml'];
        
    [h1,seg1,bscan,slo] = read_vol(file_vol, 'coordinates');
    T1 = compute_thickness(seg1,layers,h1.scale_z);

    % IOWA loading
    [h2,seg2] = read_xml_iowa(file_iowa, 'coordinates');
    T2 = compute_thickness(seg2, layers,h2.scale_z);

    for i_layer=1:n_layer
        layer = layers{i_layer};

        subplot(2,n_layer,i_layer);
        surf(h1.X_oct,h1.Y_oct,T1.(layer),'EdgeColor','none')
        view(0,90);
        axis([-3 3 -3 3]);
        title([layer ' - Spectralis']);

        subplot(2,n_layer,i_layer+n_layer);
        surf(h1.X_oct,h1.Y_oct,T2.(layer),'EdgeColor','none')
        view(0,90);
        axis([-3 3 -3 3]);
        title([layer ' - IOWA']);
    end
end
end

function test_iowa_raster(visu)
%% Raster
in_dir = '../data_private/iowa/vol_raster';

images = {'HC001AF_H_2644', 'HC001AF_H_2645',...
          'HC001AM_H_2339', 'HC001AM_H_2340',...
          'HC002AF_A_2888', 'HC002AF_A_2889'};
n_image = length(images);

layers = {'TRT','RNFL','GCIP','INL'};
n_layer = length(layers);

for i_image=1:n_image
    image = images{i_image};
    
    file_vol  = [in_dir '/' image '/' image '_OCT_Iowa.vol'];
    file_iowa = [in_dir '/' image '/' image '_Surfaces_Iowa.xml'];
        
    [h1,seg1,bscan,slo] = read_vol(file_vol, 'coordinates');
    T1 = compute_thickness(seg1,layers,h1.scale_z);

    % IOWA loading
    [h2,seg2] = read_xml_iowa(file_iowa, 'coordinates');
    T2 = compute_thickness(seg2, layers,h2.scale_z);

    for i_layer=1:n_layer
        layer = layers{i_layer};

        subplot(2,n_layer,i_layer);
        surf(h1.X_oct,h1.Y_oct,T1.(layer),'EdgeColor','none')
        view(0,90);
        axis([-3 3 -3 3]);
        title([layer ' - Spectralis']);

        subplot(2,n_layer,i_layer+n_layer);
        surf(h2.X_oct,h2.Y_oct,T2.(layer),'EdgeColor','none')
        view(0,90);
        axis([-3 3 -3 3]);
        title([layer ' - IOWA']);
    end
end
end