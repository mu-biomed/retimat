close all;clc;clearvars;
addpath(genpath('..'));

file = '../data/raster.vol';
[header, seg, bscan, slo] = read_vol(file, 'coordinates','raw_pixel');
%% Basic reflectance map (all layers)
R = reflectance_map(bscan,'raw','mean');
imagesc(R);

%% Layer reflectance map

Tnfl = reflectance_map(bscan, 'raw', 'total', seg, 'ILM', 'RNFL_GCL');
Tnfl = reflectance_map(bscan, 'normalized', 'total', seg, 'IZ_RPE', 'BM');
% Trpe = reflectance_map(bscan, 'raw', 'total', seg, 'IZ_RPE', 'BM');

att = reflectance_map(bscan, 'attenuation','mean', seg, header.scale_z, 'ILM', 'BM');

subplot(131);imagesc(Tnfl);title('RNFL raw');
subplot(132);imagesc(Trpe);title('RPE raw');
subplot(133);imagesc(att);title('TRT attenuation');
set(gca,'ColorScale','log')

colormap(gray);

%% Normalize bscan
clc;clearvars;close all;
file = '../data/raster.vol';
[header, seg, bscan, slo] = read_vol(file, 'coordinates');

bscan_norm = normalize_reflectance(bscan, seg, 'ascan');
bscan_norm2 = normalize_reflectance(bscan, seg, 'bscan');

n_bscan = 6;
idx_bscan = round(linspace(1,25,n_bscan));

norm_lim = [prctile(bscan_norm(:),1) prctile(bscan_norm(:),99)];
norm_lim2 = [prctile(bscan_norm2(:),1) prctile(bscan_norm2(:),99)];

for i=1:n_bscan
    b = idx_bscan(i);
    
    subplot(3,n_bscan,i);
    imagesc(bscan(:,:,b));
    caxis(minmax(bscan(:)'));
    title('raw');
    
    subplot(3,n_bscan,i+n_bscan);
    imagesc(bscan_norm(:,:,b));
    caxis(norm_lim);
    title('norm-ascan');
    
    subplot(3,n_bscan,i+2*n_bscan);
    imagesc(bscan_norm2(:,:,b));
    caxis(norm_lim2);
    title('norm-bscan');
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

%% Image Quality (Spectralis)
clc;clearvars;close all;
file = '../data/raster.vol';
[header, seg, bscan] = read_vol(file,'full_header');

mTCI = image_quality(bscan,'mTCI','Spectralis');
snr = image_quality(bscan,'snr',seg);

subplot(121);
scatter(header.quality, mTCI, 'filled'); hold on;
xlabel('Spectralis');
ylabel('mTCI');

subplot(122);
scatter(header.quality, snr, 'filled');
xlabel('Spectralis');
ylabel('SNR');

%% Image Quality (Topcon)
clc;clearvars;close all;
folder = '../data_private/example_tabs';
[header, seg, bscan] = read_tabs(folder);
mTCI = image_quality(bscan,'mTCI','3D-OCT-1000');
snr = image_quality(bscan,'snr',seg);
