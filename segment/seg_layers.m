function seg = seg_layers(I, scale_z, layers, visu)
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

% Get retina mask
mask_retina = seg_retina(I_flat, scale_z, 50, 'mean', false);
for i=1:M
    ind = find(mask_retina(:,i)==1);
    start = ind - 10;
    % start <1 -> start=1
    mask_retina(start:ind, i) = 1;
end

I = I_flat;

% Resize masks
% I = imresize(I_flat, 1);
% mask_retina = imresize(mask_retina, 1);
% [N, M] = size(I);

% Compute horizontal intensity gradients
n_filt = 10;
I_dl = conv2(I, [ones(1,n_filt);-ones(1,n_filt)], 'same');
I_ld = conv2(I, [-ones(1,n_filt);ones(1,n_filt)], 'same');

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
        plot(extra.path_j, extra.path_i,'r');
        colormap(gca, gray);
        set(gca,'YDir','reverse');
        axis off;
        
        subplot(n_layer,3,3+3*(i_layer-1));hold on;
        imagesc(D);
        scatter(1,1,'r')
        scatter(M+2,N,1,'r')
        plot(extra.path_j, extra.path_i,'r');
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

function [seg, D, mask, extra] = segment_layer(I_dl, I_ld, seg, mask_retina, layer)

% Basic info
[N, M] = size(mask_retina);
n_row = sum(mask_retina(:,1));

% Se
switch layer
    case 'ilm'
        % Search mask: half upper retinal mask    
        % Gradient: I_dl
        mask = mask_retina & (cumsum(mask_retina) < n_row/2);
        I = I_dl;
        
    case 'isos'
        % Search mask: half bottom retinal mask
        % Gradient: I_dl
        mask = mask_retina & (cumsum(mask_retina) > n_row/3);
        I = I_dl;
        
    case 'bm'
        % Search mask: down ISOS
        % Gradient: I_ld (I_dl does not work well)
        mask = false(N,M);
        for i=1:M
            mask(seg.isos(i)+5:seg.isos(i)+50,i) = 1; 
        end
        mask = mask_retina .* mask;
        I = I_ld;
        
    case 'elm'
        % Search mask: top ISOS
        % Gradient: I_dl
        mask = false(N,M);
        for i=1:M
            mask(seg.isos(i)-15:seg.isos(i)-4,i) = 1; 
        end
%         mask = [ones(N,1) mask_retina .* mask ones(N,1)];
        I = I_dl;

    otherwise
        error("Unknown layer to segment");
end

% Add padding
mask = [ones(N,1) mask ones(N,1)];  

% Initialize Graph
G = [I(:,1) I I(:,end)];
G(~mask) = -Inf;

% Find shortest path
D = dijkstra_matrix(G);
[path_i, path_j] = get_path(D);

seg.(layer) = path_i(path_j~=M+2 & path_j~=1);
seg.(layer) = flip(seg.(layer));

extra.path_i = path_i;
extra.path_j = path_j;

function D = dijkstra_matrix(G)
    
N = size(G, 1);
M = size(G, 2) - 2;

w_min = 1e-5;
    
V = false(N, M+2); % visited nodes
D = Inf*ones(N, M+2);  % distance to each node
D_un = 1;
idx_un = 1;

D(1,1) = 0;  % start-node distance

current_i = 1;
current_j = 1;

end_i = N;
end_j = M + 2;

stop = false;
while ~stop    
    % Find neighbours
    % Fist column: 4 neighbours (down, right-up, right, right-down)
    % Last column: 2 neighbours (up, down)
    % Regular column: 3 neighbours (right-up, right, right-down)

    % Slower alternative:  idx = sub2ind(size(D), current_i, current_j);
    idx = N * (current_j - 1) + current_i; 
    
    D_un = D_un(idx_un ~= idx);
    idx_un = idx_un(idx_un ~= idx);

    if current_j == 1
        neigh_i = current_i + [-1 0 1 1];
        neigh_j = current_j + [1 1 0 1];            
    elseif current_j == M+2     
        neigh_i = current_i + [-1 1];
        neigh_j = current_j + [0 0];
    else 
        neigh_i = current_i + [-1 0 1];
        neigh_j = current_j + [1 1 1];                       
    end

    in_box = neigh_i>0 & neigh_i<=N & neigh_j>0 & neigh_j<=M+2;    
    n_neigh = sum(in_box);
    if n_neigh==0
        disp('Finished');
        break;
    end
    neigh_i = neigh_i(in_box);
    neigh_j = neigh_j(in_box);

    % Explore neighbors
    for n=1:n_neigh
        i = neigh_i(n);
        j = neigh_j(n);

        % If visited forget about it
        if V(i,j)
            continue;
        end

        % Compute distance
        if isequal([j current_j], [1 1]) | isequal([j current_j], [M+2 M+2])
            dab = w_min;        % First/last columns vertically wmin weight
        else            
            ga = G(current_i, current_j);
            gb = G(i,j);
            dab = 2 - (ga + gb) + w_min;
        end        
        d = D(current_i, current_j) + dab;

        if d < D(i, j)
            D(i, j) = d;

            idx = N * (j - 1) + i; 
%             idx = sub2ind(size(D), i, j);
            if any(idx_un == idx)
                D_ux(idx_un==idx) = d;
            else
                D_un(end+1) = d;
                idx_un(end+1) = idx;
            end
        end        
    end

    % Mark node as visited
    V(current_i, current_j) = true;

    % Stop if we reached the end node
    if current_i==end_i & current_j==end_j
        disp('End node reached');
        break;
    end

    % Choose next node
    [~, next_node] = min(D_un);
    next_node = idx_un(next_node);
    % Slower: [current_i, current_j] = ind2sub(size(D), next_node);
    current_j = ceil(next_node/N);
    current_i = next_node - (current_j-1)*N;    
end    

function [path_i, path_j] = get_path(D)

N = size(D, 1);
M = size(D, 2) - 2; 

% Get path
path_i = nan(1,M+2);
path_j = nan(1,M+2);

current_i = N;
current_j = M+2;

c = 1;
stop = false;
while ~stop
     
    path_i(c) = current_i;
    path_j(c) = current_j;   
    
    if path_j(c) == M+2
        neigh_i = current_i + [-1 0 1 -1 1];
        neigh_j = current_j + [-1 -1 -1 0 0];
    elseif path_j(c) == 1
        neigh_i = current_i + -1;
        neigh_j = current_j + 0;           
    else
        neigh_i = current_i + [-1 0 1];
        neigh_j = current_j + [-1 -1 -1];        
    end
    
    in_box = neigh_i>0 & neigh_i<=N & neigh_j>0 & neigh_j<=M+2;
    n_neigh = sum(in_box);
    if n_neigh==0
        disp('Finished');
        break;
    end
    neigh_i = neigh_i(in_box);
    neigh_j = neigh_j(in_box);
    
    n_neigh = length(neigh_i);
    
    d = nan(1,n_neigh);
    for n=1:n_neigh
        d(n) = D(neigh_i(n), neigh_j(n));
    end
    [~,next_node] = min(d);    
    
    current_i = neigh_i(next_node);
    current_j = neigh_j(next_node);         
        
    c = c + 1;
end

function [I_flat, shift, mask_retina] = flatten_retina(I, method, scale_z)
% Pilot rpe based flattening
M = size(I, 2);

rpe = nan(1,M);
switch method
    case 'middle'
        If = imgaussfilt(I, 1.5);

        mid = nan(1,M);
        for i=1:M
            pdf = If(:,i)/sum(If(:,i));
            cdf = cumsum(pdf);
            [~, mid(i)] = min(abs(cdf - 0.4));
            
            mid_point = round(mid(i));
            [~, rpe(i)] = max([zeros(mid_point-1,1); If(mid_point:end,i)]);
        end
        mask_retina = [];
    case 'mask'
        mask_retina = seg_retina(I, scale_z, 50, 'mean', false);
        n_pixel_rnfl = 1e-3*80/scale_z;
        for i=1:M            
            first_mask = find(mask_retina(:,i), 1);
            mid_point = round(first_mask + n_pixel_rnfl);
            [~, rpe(i)] = max([zeros(mid_point-1,1); I(mid_point:end,i)]);
        end        
    otherwise
        error("Not supported retina flattening method");        
end

% Second order polynomial
p = polyfit(1:M, rpe, 2);
rpe_2 = p(3) + p(2)*(1:M) + p(1)*(1:M).^2;

% Flatten the retina
shift = round(min(rpe_2) - rpe_2);

I_flat = I;
for i=1:M
    I_flat(:,i) = circshift(I_flat(:,i), shift(i));
end
