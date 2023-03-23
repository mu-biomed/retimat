function seg = seg_layers(I, scale_z, layers, mask_retina, visu)
%SEG_LAYERS Segment retinal layers from a macular OCT B-scan
%
%   Seg = seg_layers(I, scale_z, visu)
%   Segmentation of retinal layer boundaries based on graphs [1].
%
%   Input arguments:
%  
%   'I'              B-Scan image.
%            
%   'scale_z'        Axial resolution in micrometers.
%            
%   'visu'           If true segmentation results are displayed.
%            
%  
%  
%   Output arguments:
%  
%   'seg'            Struct with the segmented values.          
%  
%
%   
%   Notes
%   -----
%   The function has been tested with Spectralis data only.
%
%
%   References
%   ----------
%   [1] Chiu et al., Automatic segmentation of seven retinal layers in SDOCT
%   images congruent with expert manual segmentation, Optics Express, 2010.
%
%
%   Example
%   ---------      
%   % Example description
%
%     [header,~,bscan] = read_vol(my_vol.vol);
%     seg = seg_layers(bscan(:,:,13), header.scale_z);
%
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2022

[N, M] = size(I);

% Flatten retina
[I_flat, shift, mask] = flatten_retina(I, 'middle', scale_z);
I = I_flat;

% Get retina mask
if mask_retina
    mask_retina = seg_retina(I_flat, scale_z, 50, 'mean', false);
    for i=1:M
        ind = find(mask_retina(:,i)==1);
        start = ind - 10;
        % start <1 -> start=1
        mask_retina(start:ind, i) = 1;
    end
else
    mask_retina = ones(size(I_flat));
end

% Resize masks
% I = imresize(I_flat, 1);
% mask_retina = imresize(mask_retina, 1);
% [N, M] = size(I);

% Compute horizontal intensity gradients
n_filt = [2 4];
I_dl = conv2(I, [ones(n_filt);-ones(n_filt)], 'same');
I_ld = conv2(I, [-ones(n_filt);ones(n_filt)], 'same');

I_dl = normalize(I_dl, 'range');
I_ld = normalize(I_ld, 'range');

% Initialize segmentation
seg = struct;

n_layer = length(layers);

for i_layer=1:n_layer
    layer = layers{i_layer};

    [seg, D, mask, extra] = segment_layer(I_dl, I_ld, seg, mask_retina, layer);
    
    if visu    
        subplot(n_layer,3,1+3*(i_layer-1)); hold on;
        imagesc(I);
        plot(seg.(layer), 'r', 'Linewidth', 1);
        colormap(gca, 'gray');
        set(gca,'YDir','reverse');
        axis off;
        
        subplot(n_layer,3,2+3*(i_layer-1));hold on;
        imagesc(I_dl);
        scatter(1,1,'r')
        scatter(M+2,N,1,'r')
        plot(extra.path_j, extra.path_i,'--r');
        colormap(gca, gray);
        set(gca,'YDir','reverse');
        axis off;
        
        subplot(n_layer,3,3+3*(i_layer-1));hold on;
        imagesc(D);
        scatter(1,1,'r')
        scatter(M+2,N,1,'r')
        plot(extra.path_j, extra.path_i,'--r');
        colors = parula;
        colormap(gca, colors(end:-1:1,:));
        set(gca,'YDir','reverse');
        title(layer);
    end
end

if visu
    figure;
    imagesc(I);hold on;
    for i=1:length(layers)
        plot(seg.(layers{i}),'LineWidth',1);
    end
    colormap(gray);
    legend(layers,'Location','southwest');
    set(gca,'FontSize',14);
    axis off;
end
