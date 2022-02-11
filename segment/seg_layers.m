function Seg = seg_layers(I, scale_z, visu)
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
%   'Seg'            Struct with the segmented values.          
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
[I_flat, shift] = flatten_retina(I);

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

% ILM segmentation
[ilm, D, mask, extra] = segment_ilm(I_dl, mask_retina);

if visu    
    % Visualization
    clf;
    subplot(431); hold on;
    imshow(I);
    plot(ilm, 'r', 'Linewidth', 1);
    colormap(gca, 'gray');
    set(gca,'YDir','reverse');
    
    subplot(432);hold on;
    imagesc(I_dl);
    scatter(1,1,'r')
    scatter(M+2,N,1,'r')
    plot(extra.path_j, extra.path_i,'r');
    colormap(gca, gray);
    set(gca,'YDir','reverse');
    
    subplot(433);hold on;
    imagesc(D);
    scatter(1,1,'r')
    scatter(M+2,N,1,'r')
    plot(extra.path_j, extra.path_i,'r');
    colors = parula;
    colormap(gca, colors(end:-1:1,:));
    set(gca,'YDir','reverse');
    title('ILM');
end

% ISOS segmentation
[isos, D, mask, extra] = segment_isos(I_dl, mask_retina);

if visu 
    subplot(434); hold on;
    imshow(I);
    plot(isos, 'g', 'Linewidth', 1);
    
    subplot(435);hold on;
    imagesc(I_dl);
    scatter(1,1,'r')
    scatter(M+2,N,1,'r')
    plot(extra.path_j, extra.path_i,'r');
    colormap(gca, gray);
    set(gca,'YDir','reverse');
    
    subplot(436);hold on;
    imagesc(D);
    scatter(1,1,'r')
    scatter(M+2,N,1,'r')
    plot(extra.path_j, extra.path_i,'r');
    colors = parula;
    colormap(gca, colors(end:-1:1,:));
    set(gca,'YDir','reverse');
end

% BM segmentation
[bm, D, mask, extra] = segment_bm(I_ld, mask_retina, isos);

if visu 
    subplot(437); hold on;
    imshow(I);
    plot(bm, 'g', 'Linewidth', 1);
    
    subplot(438);hold on;
    imagesc(I_dl);
    scatter(1,1,'r')
    scatter(M+2,N,1,'r')
    plot(extra.path_j, extra.path_i,'r');
    colormap(gca, gray);
    set(gca,'YDir','reverse');
    
    subplot(439);hold on;
    imagesc(D);
    scatter(1,1,'r')
    scatter(M+2,N,1,'r')
    plot(extra.path_j, extra.path_i,'r');
    colors = parula;
    colormap(gca, colors(end:-1:1,:));
    set(gca,'YDir','reverse');
end

% ELM segmentation
[elm, D, mask, extra] = segment_elm(I_dl, mask_retina, isos);

if visu 
    subplot(4,3,10); hold on;
    imshow(I);
    plot(elm, 'g', 'Linewidth', 1);
    
    subplot(4,3,11);hold on;
    imagesc(I_dl);
    scatter(1,1,'r')
    scatter(M+2,N,1,'r')
    plot(extra.path_j, extra.path_i,'r');
    colormap(gca, gray);
    set(gca,'YDir','reverse');
    
    subplot(4,3,12);hold on;
    imagesc(D);
    scatter(1,1,'r')
    scatter(M+2,N,1,'r')
    plot(extra.path_j, extra.path_i,'r');
    colors = parula;
    colormap(gca, colors(end:-1:1,:));
    set(gca,'YDir','reverse');
end

end

function [ilm, D, mask, extra] = segment_ilm(I_dl, mask_retina)

[N, M] = size(mask_retina);

n_row = sum(mask_retina(:,1));

mask = mask_retina & (cumsum(mask_retina) < n_row/2);

% Search mask: half upper retinal mask
mask = [ones(N,1) mask ones(N,1)];

% Initialize Graph
G = [I_dl(:,1) I_dl I_dl(:,end)];
G(~mask) = -Inf;

% Find shortest path
D = dijkstra_matrix(G);
[path_i, path_j] = get_path(D);

ilm = path_i(path_j~=M+2 & path_j~=1);
ilm = flip(ilm);

extra.path_i = path_i;
extra.path_j = path_j;
end

function [isos, D, mask, extra] = segment_isos(I_dl, mask_retina)

[N, M] = size(mask_retina);

n_row = sum(mask_retina(:,1));

mask = mask_retina & (cumsum(mask_retina) > n_row/3);

% Search mask: half upper retinal mask
mask = [ones(N,1) mask ones(N,1)];

% Initialize Graph
G = [I_dl(:,1) I_dl I_dl(:,end)];
G(~mask) = -Inf;

% Find shortest path
D = dijkstra_matrix(G);
[path_i, path_j] = get_path(D);

isos = path_i(path_j~=M+2 & path_j~=1);
isos = flip(isos);

extra.path_i = path_i;
extra.path_j = path_j;
end

function [bm, D, mask, extra] = segment_bm(I_ld, mask_retina, isos)
[N, M] = size(mask_retina);

mask_down = false(N,M);
for i=1:M
    mask_down(isos(i)+5:isos(i)+50,i) = 1; 
end
mask = [ones(N,1) mask_retina .* mask_down ones(N,1)];

G = [I_ld(:,1) I_ld I_ld(:,end)];  % I_dl does not work well
G(~mask) = -Inf;
D = dijkstra_matrix(G);
[path_i, path_j] = get_path(D);
bm = flip(path_i(path_j~=M+2 & path_j~=1));

extra.path_i = path_i;
extra.path_j = path_j;
end

function [elm, D, mask, extra] = segment_elm(I_dl, mask_retina, isos)

[N, M] = size(mask_retina);

% ELM segmentation
mask = false(N,M);
for i=1:M
    mask(isos(i)-15:isos(i)-4,i) = 1; 
end
% mask = [ones(N,1) mask_retina .* mask_up ones(N,1)];
G = [I_dl(:,1) I_dl I_dl(:,end)];
G(~mask) = -Inf;
D = dijkstra_matrix(G);
[path_i, path_j] = get_path(D);
elm = flip(path_i(path_j~=M+2 & path_j~=1));

extra.path_i = path_i;
extra.path_j = path_j;
end

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

end

function [I_flat, shift] = flatten_retina(I)
% Pilot rpe based flattening

M = size(I, 2);

If = imgaussfilt(I, 1.5);

rpe = nan(1,M);
mid = nan(1,M);
for i=1:M
    pdf = If(:,i)/sum(If(:,i));
    cdf = cumsum(pdf);
    [~, mid(i)] = min(abs(cdf - 0.4));
    
    mid_point = round(mid(i));
    [~, rpe(i)] = max([zeros(mid_point-1,1); If(mid_point:end,i)]);
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

end