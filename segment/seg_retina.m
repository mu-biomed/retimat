close all;clc;clearvars;
addpath(genpath('..'));

file = '../data/raster.vol';

[header, seg, bscan,~] = read_vol(file);

I = bscan(:,:,13);
[N, M] = size(I);

subplot(141);
imagesc(I);colormap(gray);
%% K-Means
mask = nan(N,M);
for n=1:M
    X1 = I(:,n);
    X2 = (1:N)';
    X = [X1 X2];
    X = (X - mean(X))./std(X);
    mask(:,n) = kmeans(X, 3);
%     scatter(X1, X2)
end

subplot(142);
imagesc(mask);
%% OTSU 1
T = graythresh(I);
mask = zeros(size(I));
mask(I<=T) = 1;
subplot(143);imagesc(mask)

%% OTSU 2 levels
thresh = multithresh(I, 2);

mask = nan(size(I));
mask(I <=thresh(1)) = 1;
mask(I > thresh(1) & I <=thresh(2)) = 2;
mask(I >thresh(2)) = 3;

subplot(144);
imagesc(mask)
