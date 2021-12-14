close all;clc;clearvars;
addpath(genpath('..'));

file = '../data/raster.vol';

[header, seg, bscan,~] = read_vol(file);

I = bscan(:,:,13);
[N, M] = size(I);
%% Pilot RPE

% Gaussian filter
If = imgaussfilt(I, 1.5);
% imshow(If);

rpe = nan(1,M);
mid = nan(1,M);
for i=1:M
    pdf = If(:,i)/sum(If(:,i));
    cdf = cumsum(pdf);
    [~, mid(i)] = min(abs(cdf - 0.4));
    
    mid_point = round(mid(i));
    [~, rpe(i)] = max([zeros(mid_point-1,1); If(mid_point:end,i)]);
%     centroid(i) = (1:N) * pdf;
end

% imshow(I);hold on;
% plot(rpe);
% plot(mid);

% Second order polynomial
p = polyfit(1:M, rpe, 2);
rpe_2 = p(3) + p(2)*(1:M) + p(1)*(1:M).^2;
% plot(rpe_2, 'r');


%% Flatten the retina
shift = round(min(rpe_2) - rpe_2);

for i=1:M
    I_flat(:,i) = circshift(I(:,i), shift(i));
end

I = I_flat;
% figure;imshow(I_flat);

%% Resizing

% I = I(10:230,:);
% I = imresize(I, [size(I,1) 200]);
[N,M] = size(I);

%% Gradients
I_dl = conv2(I, [1;-1], 'same');
I_ld = conv2(I, [-1;1], 'same');

I_dl = (I_dl - min(I_dl(:)))/(max(I_dl(:)) - min(I_dl(:)));
% I_ld = I_ld./max(I_ld(:));

% subplot(121);imshow(I_dl);
% subplot(122);imshow(I_ld);

%% Construct map
colors = parula;

G = [I_dl(:,1) I_dl I_dl(:,end)];
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
    D_unvisit = D;
    D_unvisit(V) = Inf;    
    [~, next_node] = min(D_unvisit(:));
    [current_i, current_j] = ind2sub(size(D), next_node);
end

clf;
subplot(211);imagesc(I);hold on;
colormap(gca, 'gray');

subplot(212);imagesc(D);hold on;
scatter(1,1,'r')
scatter(M+2,N,1,'r')
colormap(gca, colors(end:-1:1,:));

%% Get path
path_i = nan(1,M+2);
path_j = nan(1,M+2);

current_i = end_i;
current_j = end_j;

c = 1;
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
    
%     scatter(current_j, current_i,'r');
    
    c = c + 1;
end

plot(path_j, path_i,'r');
subplot(211);
plot(path_j, path_i,'r');
