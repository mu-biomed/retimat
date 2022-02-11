close all;clc;clearvars;
addpath(genpath('../'));
file = '../data/raster.vol';

[header, ~, bscan,~] = read_vol(file);


I = bscan(:,:,13);
[N, M] = size(I);
mask_retina = seg_retina(I, 0.0039);

I = adapthisteq(I);

% Average image
% I1 = imgaussfilt(I, 2);
I1 = imfilter(I, ones(3,19));
% Dark pixels in column to 0
I2 = I1;
for n=1:M
%     low_px = I1(:,n) < mean(I1(:,n));
    low_px = I1(:,n) < 0.8*median(I1(mask_retina(:,n),n));
    I2(low_px,n) = 0; 
end
% Second derivative (columnwise)
I3 = diff(I2, 2, 1);
[N,M] = size(I3);
% Mask
mask = I3;
mask(mask>0) = 0;
mask(mask<0) = 1;

% Remove small clusters (columnwise)
mask2 = mask;
for m=1:M
    con = bwconncomp(mask(:,m));
    len = cellfun(@(x) length(x), con.PixelIdxList);
    small = con.PixelIdxList(len < 5);
    for i=1:length(small)
        mask2(small{i},m) = 0; 
    end    
end

% Join close clusters (columnwise)
mask3 = mask2;
for m=1:M
    con = bwconncomp(mask2(:,m));
    first = cellfun(@(x) x(1), con.PixelIdxList);
    last = cellfun(@(x) x(end), con.PixelIdxList);

    dist = first(2:end) - last(1:end-1);
    gap_first = last(1:end-1) + 1;
    gap_last = first(2:end) - 1;
    gap_len = gap_last - gap_first + 1;

    th = 3;
    small_gap_first = gap_first(gap_len < th);
    small_gap_last = gap_last(gap_len < th);
    n_gap = length(small_gap_first);

    for n=1:n_gap
        mask3(small_gap_first(n):small_gap_last(n),m) = 1;
    end
end

n = 2; m = 4;
subplot(n,m,1);imagesc(I);colormap(gray);
subplot(n,m,2);imagesc(I1);colormap(gray);
subplot(n,m,3);imagesc(I2);colormap(gray);
subplot(n,m,4);imagesc(I3);colormap(gray);
subplot(n,m,5);imagesc(mask);colormap(gray);
subplot(n,m,6);imagesc(mask2);colormap(gray);
subplot(n,m,7);imagesc(mask3);colormap(gray);
