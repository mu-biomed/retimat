close all;clc;clearvars;
addpath(genpath('..'));

global visu
visu = true;

files.vol_raster = '../data_private/vol/raster/HC001AF_H_2644.vol';
files.vol_star   = '../data_private/vol/star/HC001AM_H_2341.vol';
files.vol_wide   = '../data_private/vol/wide/Heidelberg_Macula.vol';
files.vol_onh    = '../data_private/vol/onh/HC001AM_H_2343_0.vol';
files.img_mac    = '../data_private/img/macula/PNYU006E_Macular Cube 512x128_2-19-2016_11-49-59_OD_sn8434_cube_raw.img';
files.bin_mac    = '../data_private/img/macula/PNYU006E_Macular Cube 512x128_2-19-2016_11-49-59_OD_sn8434_lslo.bin';
files.img_onh    = '../data_private/img/onh/PNYU001E_Optic Disc Cube 200x200_11-19-2015_14-22-31_OD_sn7997_cube_raw.img';
files.bin_onh    = '../data_private/img/onh/PNYU001E_Optic Disc Cube 200x200_11-19-2015_14-22-31_OD_sn7997_lslo.bin';
files.e2e_mac    = '../data_private/e2e/oct_1.e2e';
files.fda_mac    = '../data_private/fda/test_0.fda';

files.iowa_star.folder = '../data_private/iowa/vol_star';
files.iowa_star.ids    = {'HC001AM_H_2341', 'HC001AM_H_2342',...
                          'HC002AM_H_2347', 'HC002AM_H_2348',...
                          'HC004AM_H_2357', 'HC004AM_H_2358'};
                
files.iowa_mac.folder = '../data_private/iowa/vol_raster';
files.iowa_mac.ids    = {'HC001AF_H_2644', 'HC001AF_H_2645',...
                         'HC001AM_H_2339', 'HC001AM_H_2340',...
                         'HC002AF_A_2888', 'HC002AF_A_2889'};                
              
files.iowa_onh.folder = '../data_private/iowa/img_onh'; 
files.iowa_onh.ids    = {'PNYU001E_Optic Disc Cube 200x200_11-19-2015_14-22-31_OD_sn7997_cube_raw',...
                         'PNYU001E_Optic Disc Cube 200x200_11-19-2015_14-27-45_OS_sn8002_cube_raw'};

% Launch test suite
test_vol_raster(files);
test_vol_star(files);
test_vol_wide(files);
test_vol_onh(files);
test_img_macular_cube(files);
test_img_onh_cube(files);
test_e2e(files);
test_fda(files);
test_iowa_star(files);
test_iowa_raster(files);
test_iowa_onh(files);

close all;

function test_vol_raster(files)
[header, seg, bscan, fundus] = read_vol(files.vol_raster, 'get_coordinates');
plot_vol_data(header, seg, bscan, fundus);
end

function test_vol_star(files)
[header, seg, bscan, fundus] = read_vol(files.vol_star, 'get_coordinates');
plot_vol_data(header, seg, bscan, fundus)
end

function test_vol_wide(files)
[header, seg, bscan, fundus] = read_vol(files.vol_wide, 'full_header', 'get_coordinates');
plot_vol_data(header, seg, bscan, fundus)
end

function test_vol_onh(files)
[header, seg, bscan, fundus] = read_vol(files.vol_onh, 'verbose', 'get_coordinates');
plot_vol_data(header, seg, bscan, fundus)
end

function test_img_macular_cube(files)
[header, bscan] = read_img(files.img_mac, [], 'get_coordinates');
fundus          = read_bin(files.bin_mac);
plot_img_data(header, bscan, fundus)
end

function test_img_onh_cube(files)
[header, bscan] = read_img(files.img_onh, [], 'get_coordinates');
fundus          = read_bin(files.bin_onh);
plot_img_data(header, bscan, fundus)
end

function test_e2e(files)
[header, seg, bscan, fundus] = read_e2e(files.e2e_mac, 'verbose');
plot_e2e_data(header{1}, seg{1}, bscan{1}, fundus{1});
end

function test_fda(files)
[header, seg, bscan, fundus] = read_fda(files.fda_mac, 'verbose', 'get_coordinates');
plot_fda_data(header, seg, bscan, fundus);
end

function test_iowa_star(files)
%% Star

n_image = length(files.iowa_star.ids);
in_dir = files.iowa_star.folder;

layers = {'TRT','RNFL','GCIPL','INL'};
n_layer = length(layers);

for i_image=1:n_image
    image = files.iowa_star.ids{i_image};
    
    file_vol  = [in_dir '/' image '/' image '_OCT_Iowa.vol'];
    file_iowa = [in_dir '/' image '/' image '_Surfaces_Iowa.xml'];
        
    [h1, seg1] = read_vol(file_vol, 'get_coordinates');
    T1 = compute_thickness(seg1,layers,h1.scale_z);

    % IOWA loading
    [h2, seg2] = read_xml_iowa(file_iowa);
    h2.bscan_pattern = 'star';
    [h2.X_oct, h2.Y_oct] = get_ascan_coordinates(h2);

    T2 = compute_thickness(seg2, layers, h2.scale_z);

    figure;
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
    sgtitle(file_iowa, 'Interpreter', 'none');
end
end

function test_iowa_raster(files)

n_image = length(files.iowa_mac.ids);
in_dir = files.iowa_mac.folder;

layers = {'TRT','RNFL','GCIPL','INL'};
n_layer = length(layers);

for i_image=1:n_image
    image = files.iowa_mac.ids{i_image};
    
    file_vol  = [in_dir '/' image '/' image '_OCT_Iowa.vol'];
    file_iowa = [in_dir '/' image '/' image '_Surfaces_Iowa.xml'];
        
    [h1, seg1] = read_vol(file_vol, 'get_coordinates');
    T1 = compute_thickness(seg1,layers,h1.scale_z);

    % IOWA loading
    [h2,seg2] = read_xml_iowa(file_iowa, 'get_coordinates');
    T2 = compute_thickness(seg2, layers, h2.scale_z);

    figure;
    for i_layer=1:n_layer
        layer = layers{i_layer};

        subplot(1,n_layer,i_layer);
        surf(h1.X_oct,h1.Y_oct,T1.(layer),'EdgeColor','none')
        view(0,90);
        axis([-3 3 -3 3]);
        title([layer ' - IOWA']);
    end
    sgtitle(file_iowa, 'Interpreter', 'none');    
end
end

function test_iowa_onh(files)
n_image = length(files.iowa_onh.ids);
in_dir = files.iowa_onh.folder;

layers = {'TRT','RNFL','GCIPL','INL'};
n_layer = length(layers);

for i_image=1:n_image
    image = files.iowa_onh.ids{i_image};
    
    file_iowa = [in_dir '/' image '/' image '_Surfaces_Iowa.xml'];
        
    % IOWA loading
    [h, seg] = read_xml_iowa(file_iowa, 'get_coordinates');
    T = compute_thickness(seg, layers, h.scale_z);

    figure;
    for i_layer=1:n_layer
        layer = layers{i_layer};

        subplot(2,n_layer,i_layer);
        surf(h.X_oct,h.Y_oct,T.(layer),'EdgeColor','none')
        view(0,90);
        axis([-3 3 -3 3]);
        title([layer ' - Spectralis']);

        subplot(2,n_layer,i_layer+n_layer);
        surf(h.X_oct,h.Y_oct,T.(layer),'EdgeColor','none')
        view(0,90);
        axis([-3 3 -3 3]);
        title([layer ' - IOWA']);
    end
    sgtitle(file_iowa, 'Interpreter', 'none');    
end

end

function plot_vol_data(header, seg, bscan, fundus)
global visu
if ~visu
    return;
end

X_fun = header.X_fun;
Y_fun = header.Y_fun;
Z_fun = max(fundus(:)) + 1;

X_oct = header.X_oct;
Y_oct = header.Y_oct;
Z_oct = repmat(Z_fun, 1, numel(header.X_oct));

figure();

% Fundus
subplot(151); hold on;
surf(X_fun, Y_fun, fundus, 'EdgeColor', 'none'); 
scatter3(X_oct(:), Y_oct(:), Z_oct, 5, 'b', 'filled');
scatter3(0, 0, Z_fun, 'r', 'filled');
if strcmp(header.fixation, 'onh')
    scatter3(X_oct(1), Y_oct(1), Z_fun, 'g', 'filled');
end
view(0,90); 
daspect([1 1 1]);
title('fundus');

% En-face
if ~strcmp(header.fixation, 'onh')
    bscan = double(bscan);
    en_face = squeeze(mean(bscan, 'omitnan')).';
    subplot(152);
    imagesc(en_face);
    daspect([header.scale_y header.scale_x 1]);
    title('en face');
end

% B-scans
layers = fields(seg);
idx_bscan = round(linspace(1, header.n_bscan, 3));
for i=1:length(idx_bscan)
    subplot(1,5,2+i); hold on;
    imshow(bscan(:,:,idx_bscan(i)));
    for i_layer=1:length(layers)
        plot(seg.(layers{i_layer})(idx_bscan(i),:));
    end
    title(sprintf('bscan : %d', idx_bscan(i)));
end

colormap(gray);
end

function plot_img_data(header, bscan, fundus)
global visu
if ~visu
    return
end

subplot(161);
imagesc(fundus);
title('fundus');
daspect([1 1 1]);

subplot(162);
en_face = squeeze(mean(bscan)).';
imagesc(en_face);
daspect([header.scale_y header.scale_x 1]);
title('en face');

subplot(163);
surf(header.X_oct, header.Y_oct, en_face, 'EdgeColor', 'none');
view(0,90);
title('en face (coords)');
daspect([1 1 1]);

idx_bscan = round(linspace(1, header.n_bscan, 3));
for i=1:length(idx_bscan)
    subplot(1,6,3+i);
    imshow(bscan(:,:,idx_bscan(i)));
    title(sprintf('Bscan: %d', idx_bscan(i)));
end
colormap(gray);
end

function plot_e2e_data(header, seg, bscan, fundus)
global visu
if ~visu
    return;
end

figure();

% Fundus
subplot(151); hold on;
imagesc(fundus);
daspect([1 1 1]);
title('fundus');

% En-face
if ~strcmp(header.fixation, 'onh')
    bscan = double(bscan);
    en_face = squeeze(mean(bscan, 'omitnan')).';
    subplot(152);
    imagesc(en_face);
    daspect([header.n_ascan header.n_bscan 1]);
    title('en face');
    caxis([0 prctile(en_face(:),99)]);
end

% B-scans
layers = fields(seg);
idx_bscan = round(linspace(1, double(header.n_bscan), 3));
for i=1:length(idx_bscan)
    subplot(1,5,2+i); hold on;
    imshow(bscan(:,:,idx_bscan(i)));
    for i_layer=1:length(layers)
        plot(seg.(layers{i_layer})(idx_bscan(i),:));
    end
    title(sprintf('bscan : %d', idx_bscan(i)));
end

colormap(gray);
end

function plot_fda_data(header, seg, bscan, fundus)
global visu
if ~visu
    return;
end
en_face = squeeze(mean(bscan,1)).';
trt = seg.ILM - seg.BM;

n_col = 7;

subplot(1,n_col,1);
imshow(fundus);

subplot(1,n_col,2);
imagesc(en_face);
colormap(gca, 'gray');
daspect([header.scale_y header.scale_x 1]);

subplot(1,n_col,3);
imagesc(trt);
daspect([header.scale_y header.scale_x 1]);

subplot(1,n_col,4);
surf(header.X_oct, header.Y_oct, trt,'EdgeColor','none');
view(0,90);
daspect([1 1 1]);

idx_bscan = round(linspace(1, header.n_bscan, 3));
layers = fields(seg);
for i=1:length(idx_bscan)
    subplot(1,n_col,4+i); hold on;
    imagesc(bscan(:,:,idx_bscan(i)));
    set(gca,'YDir','reverse');
    for i_layer=1:length(layers)
        plot(seg.(layers{i_layer})(idx_bscan(i), :));        
    end
    colormap(gca, 'gray');
    title(sprintf('Bscan:%d', idx_bscan(i)));
end
end
