close all;clc;clearvars;
addpath(genpath('../'));

file = '../data/raster.vol';

[header, ~, bscan,~] = read_vol(file);

I = bscan(:,:,1);
[N, M] = size(I);
mask_retina = seg_retina(I, 0.0039);

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

% Remove outliers
mask2 = mask;
for m=1:M
    run_length = 0;
    for n=1:N        
        if mask(n,m) == 0
            if run_length > 0 && run_length < 5
                mask2(n-run_length:n-1,m) = 0;                
            end
            run_length = 0;
        else
            run_length = run_length + 1;
        end
    end     
end

n = 2; m = 3;
subplot(n,m,1);imagesc(I);colormap(gray);
subplot(n,m,2);imagesc(I1);colormap(gray);
subplot(n,m,3);imagesc(I2);colormap(gray);
subplot(n,m,4);imagesc(I3);colormap(gray);
subplot(n,m,5);imagesc(mask);colormap(gray);
subplot(n,m,6);imagesc(mask2);colormap(gray);


