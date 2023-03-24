function [seg, D, mask, extra] = segment_layer(I, mask)
% I: image gradient
% mask: to avoid regions

[N, M] = size(I);

% Add padding
mask = [ones(N,1) mask ones(N,1)];  

% Initialize Graph
G = [I(:,1) I I(:,end)];
G(~mask) = -Inf;

% Find shortest path
D = dijkstra_matrix(G);
[path_i, path_j] = get_path(D);

seg = path_i(path_j~=M+2 & path_j~=1);
seg = flip(seg);

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
                D_un(idx_un==idx) = d;
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
