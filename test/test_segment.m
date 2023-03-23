close all;clc;clearvars;

addpath(genpath('..'));

[header, ~, bscan,~] = read_vol(file);

%% Retina segmentation
I = bscan(:,:,1);

mask = seg_retina(I, header.scale_z, 150, 'otsu', true);

%% Layers segmentation
file = '../data/raster.vol';
[header, ~, bscan,~] = read_vol(file);

% file = '../data/Zeiss_Macula.img';
% [bscan, header] = read_img(file);

for ib=1:20
    I = bscan(:,:,ib);
    I = imresize(I, [496 512]);

    n_filt = [2 4];
    I_dl = conv2(I, [ones(n_filt);-ones(n_filt)], 'same');
    I_ld = conv2(I, [-ones(n_filt);ones(n_filt)], 'same');
    
    I_dl = normalize(I_dl, 'range');
    I_ld = normalize(I_ld, 'range');
    
    mask = true(size(I));
    [seg, D, ~, extra] = segment_layer(I_dl, struct(), mask, 'layer_1');
    for i=1:width(I)
        mask(seg.layer_1(i)-10:seg.layer_1(i)+10,i) = false; 
    end
    [seg, D, mask, extra] = segment_layer(I_dl, seg, mask, 'layer_2');

    if mean(seg.layer_1(:)) > mean(seg.layer_2(:))
        seg.isos = seg.layer_1;
        seg.ilm  = seg.layer_2;
    else
        seg.ilm  = seg.layer_1;
        seg.isos = seg.layer_2;
    end

    for i=1:width(I)
        mask(1:seg.isos(i), i) = false; 
    end
    [seg, D, mask, extra] = segment_layer(I_ld, seg, mask, 'bm');

    mask = false(size(I));
    for i=1:width(I)
        mask(seg.ilm(i)+1:seg.isos(i)-5, i) = true; 
    end
    [seg, D, mask, extra] = segment_layer(I_ld, seg, mask, 'layer_4');
    
    subplot(4,6,ib);
    imagesc(I); hold on; axis off;
    plot(seg.ilm,'r');
    plot(seg.isos,'b');
    plot(seg.bm,'g');
    plot(seg.layer_4,'y');
end
colormap(gray);

% layers = {'isos','ilm'};
% seg = seg_layers(I, header.scale_z, layers, false, true);

% figure;
% imagesc(I);hold on;
% for i=1:length(layers)
%     plot(seg.(layers{i}));
% end
% colormap(gray);
% legend(layers);
