
clc;clearvars;close all;
addpath(genpath('..'));
file = '../data/raster.vol';

[~,~,~,I] = read_vol(file);
I = double(I);
[N,M] = size(I);

%% Intensity normalization
close all;

subplot(231);imagesc(I);
subplot(234);histogram(I(:));

subplot(232);
% Background subtraction
nk = 100;
h = 1/nk^2*ones(nk,nk);
back = imfilter(I, h, 'replicate');
I_cor = I- back + mean(back(:));

% Homomorphic filtering
% f = fft2(I);
% 
% If = exp(imgaussfilt(log10(I),1));
% C = mean(I(:))/mean(I(:)./If(:));
% I_cor = If * C;

I_cor = I_cor./max(I_cor(:));
imagesc(adapthisteq(I_cor));
colormap(gray);

subplot(235);
histogram(I_cor(:));

subplot(233);
th = graythresh(I_cor);
th = prctile(I_cor(:),15);
I_th = I_cor < th;
imagesc(I_th);

subplot(236);
se = ones(6,6);
I_s = imclose(I_th, se);
CC = bwconncomp(I_s);
I_s = reshape(I_s,1,[]);
for n=1:CC.NumObjects
    pixels = CC.PixelIdxList{n};
    if length(pixels) < 100
        I_s(pixels) = 0;
    end
end
I_s = reshape(I_s, [N M]);
imagesc(I_s);
