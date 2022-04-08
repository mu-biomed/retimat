close all;clc;clearvars;
addpath(genpath('..'));

%% Reflectance map
file = '../data/raster.vol';
[header, seg, bscan, slo] = read_vol(file, 'coordinates');

% MR = reflectance_map(bscan, seg, 'mean', 'ILM', 'BM');

Tnfl = reflectance_map(bscan, seg, 'total', 'ILM', 'RNFL_GCL');
Trpe = reflectance_map(bscan, seg, 'total', 'IDZ_RPE', 'BM');

att = reflectance_map(bscan, seg, 'NFL_RPE_att', header.scale_z, 'ILM', 'BM');

subplot(131);imagesc(Tnfl);
subplot(132);imagesc(Trpe);
subplot(133);imagesc(att);caxis([0 10]);
colormap(gray);

%% Normalize bscan
bscan_norm = normalize_reflectance(bscan, seg, 'column');
bscan_norm2 = normalize_reflectance(bscan, seg, 'bscan');

n_bscan = 6;
idx_bscan = round(linspace(1,25,n_bscan));

norm_lim = [prctile(bscan_norm(:),1) prctile(bscan_norm(:),99)];
norm_lim2 = [prctile(bscan_norm2(:),1) prctile(bscan_norm2(:),99)];

for i=1:n_bscan
    b = idx_bscan(i);
    
    subplot(3,n_bscan,i);imagesc(bscan(:,:,b));caxis(minmax(bscan(:)'));
    subplot(3,n_bscan,i+n_bscan);imagesc(bscan_norm(:,:,b));caxis(norm_lim);
    subplot(3,n_bscan,i+2*n_bscan);imagesc(bscan_norm2(:,:,b));caxis(norm_lim2);
end

colormap(gray);

%% Attenuation coefficient
clc;clearvars;close all;
file = '../data/raster.vol';
% [header, ~, bscan] = read_vol(file, 'coordinates');
[header, ~, bscan] = read_vol(file, 'coordinates','raw_pixel');

att = compute_attenuation(bscan, header.scale_z);
subplot(121);imagesc(bscan(:,:,13).^0.25);axis off;
subplot(122);imagesc(att(:,:,13));axis off;
colormap(gray);
set(gca,'ColorScale','log')

%% Stakced b-scans
clc;clearvars;close all;
addpath(genpath('..'));
file = '../data/raster.vol';
[header, seg, bscan] = read_vol(file);

stacked = stack_bscans(bscan, seg, {'RNFL','GCIP','RPE'});
subplot(131);imagesc(stacked.RNFL);axis off;title('RNFL');
subplot(132);imagesc(stacked.GCIP);axis off; title('GCIP');
subplot(133);imagesc(stacked.RPE);axis off; title('RPE');
colormap(gray);
