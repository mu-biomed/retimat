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

imshow(I);hold on;
plot(rpe);
plot(mid);

% Second order polynomial
p = polyfit(1:M, rpe, 2);
rpe_2 = p(3) + p(2)*(1:M) + p(1)*(1:M).^2;
plot(rpe_2, 'r');


%% Flatten the retina
shift = round(min(rpe_2) - rpe_2);

for i=1:M
    I_flat(:,i) = circshift(I(:,i), shift(i));
end

figure;imshow(I_flat);

%% Gradient images

I_dl = conv2(I, [1;-1], 'same');
I_ld = conv2(I, [-1;1], 'same');

I_dl = I_dl./max(I_dl(:));
I_ld = I_ld./max(I_ld(:));

subplot(121);imshow(I_dl);
subplot(122);imshow(I_ld);

%% Weighted graph
W = sparse(N*M, N*M);
w_min = 1e-5;

I_dlf = reshape(I_dl,1,[]);  % flatten version (for speed)
for i=1:N
    for j=1:M
        a = N*(i-1) + j;
        
        % Only explore 8 neighbors
        ib = [i-1 i-1 i-1 i i i+1 i+1 i+1];
        jb = [j-1 j j+1 j-1 j+1 j-1 j j+1];
        
        ga = I_dl(i, j);
        
        in_box = ib>0 & jb>0 & ib<=N & jb<=M;
        
        b = N*(ib(in_box)-1) + jb(in_box);
        gb = I_dlf(b);
        W(a, b) = 2 - (ga + gb) + w_min;
    end
    disp(i)
end

% Add extra first/end columns to initialize segmentation
W = [w_min*ones(N*M,1) W w_min*ones(N*M,1)];

%% Weighted graph fast
W = sparse(N*M, N*M);

% Eight different kernels
k_ul = [-1 0 0; 0 1 0; 0 0 0]; % upper-left
k_um = [0 -1 0; 0 1 0; 0 0 0]; % upper-middle
k_ur = [0 0 -1; 0 1 0; 0 0 0]; % upper-right
k_ml = [0 0 0; -1 1 0; 0 0 0]; % middle-left
k_mr = [0 0 0; 0 1 -1; 0 0 0]; % middle-right
k_dl = [0 0 0; 0 1 0; -1 0 0]; % down-left
k_dm = [0 0 0; 0 1 0; 0 -1 0]; % down-middle
k_dr = [0 0 0; 0 1 0; 0 0 -1]; % down-right

W_ul = conv2(I_dl, k_ul,'same');
