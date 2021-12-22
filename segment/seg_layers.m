function seg_layers
close all;clc;clearvars;
addpath(genpath('..'));

file = '../data/raster.vol';

[header, ~, bscan,~] = read_vol(file);

I = bscan(:,:,1);
[N, M] = size(I);

% Retina flatten
[I_flat, shift] = flatten_retina(I);

% Get retina mask
mask_retina = seg_retina(I_flat, header.scale_z, 50, 'mean', false);
for i=1:M
    ind = find(mask_retina(:,i)==1);
    start = ind - 10;
    % start <1 -> start=1
    mask_retina(start:ind, i) = 1;
end

I = imresize(I_flat, 1);
mask_retina = imresize(mask_retina, 1);
[N,M] = size(I);

% Gradients
I_dl = conv2(I, [1;-1], 'same');
I_ld = conv2(I, [-1;1], 'same');

I_dl = normalize(I_dl, 'range');
I_ld = normalize(I_ld, 'range');

% 1st layer segmentation
mask = [ones(N,1) mask_retina ones(N,1)];

G = [I_dl(:,1) I_dl I_dl(:,end)];
G(~mask) = -Inf;
D = dijkstra_matrix(G);
[path_i, path_j] = get_path(D);

layer = path_i(path_j~=M+2 & path_j~=1);
layer = flip(layer);

If = imgaussfilt(I, 0.5);
th =  graythresh(If);
bright_mask = If > th;
bright_up = zeros(1,M);
for i=1:M
    bright_up(i) = sum(bright_mask(1:layer(i), i));
end

if sum(bright_up)/sum(bright_mask(:)) > 0.025
    Layers.ISOS = layer;    
    next_layer = 'RNFL';
else
    Layers.RNFL = layer;
    next_layer = 'ISOS';
end

clf;
subplot(131); hold on;
imshow(I);
plot(layer, 'r', 'Linewidth', 1);
colormap(gca, 'gray');
set(gca,'YDir','reverse');

subplot(132);hold on;
imagesc(D);
scatter(1,1,'r')
scatter(M+2,N,1,'r')
plot(path_j, path_i,'r');
colors = parula;
colormap(gca, colors(end:-1:1,:));
set(gca,'YDir','reverse');

% 2nd layer segmentation

switch next_layer
    case 'RNFL'
        mask_up = false(N, M);
        for i=1:M
            mask_up(1:Layers.ISOS(i)-5,i) = 1; 
        end
        mask = [ones(N,1) mask_retina .* mask_up ones(N,1)];
        
    case 'ISOS'
        mask_down = false(N, M);
        for i=1:M
            mask_down(Layers.RNFL(i)+5:end,i) = 1; 
        end
        mask = [ones(N,1) mask_retina .* mask_down ones(N,1)];
end

G = [I_dl(:,1) I_dl I_dl(:,end)];
G(~mask) = -Inf;
D = dijkstra_matrix(G);
[path_i, path_j] = get_path(D);

layer = flip(path_i(path_j~=M+2 & path_j~=1));
        
subplot(131); hold on;
plot(layer, 'g', 'Linewidth', 1);
        
subplot(133);hold on;
imagesc(D);
scatter(1,1,'r')
scatter(M+2,N,1,'r')
plot(path_j, path_i,'r');
colors = parula;
colormap(gca, colors(end:-1:1,:));
set(gca,'YDir','reverse');
end

function D = dijkstra_matrix(G)
    
N = size(G, 1);
M = size(G, 2) - 2;

w_min = 1e-5;
    
V = false(N, M+2); % visited nodes
D = Inf*ones(N, M+2);  % distance to each node
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
% D_unvisit = D;
% D_unvisit(V) = Inf;    
% [~, next_node] = min(D_unvisit(:));
[~, next_node] = min(1e3.*V(:) + D(:));
[current_i, current_j] = ind2sub(size(D), next_node);
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