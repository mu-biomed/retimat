close all;clc;clearvars;
addpath(genpath('..'));

%% Fractal Dimension
I = imread('../data/D33.gif');
FD = compute_FD(I, true);

%% GLCM features
close all;clc;clearvars;

I = randn(200, 200);
GLCM = graycomatrix(I,'NumLevels',250);
X = compute_glcm_metrics(GLCM);

%% Lacunarity
close all;clc;clearvars;
I = imread('cameraman.tif');
I_noise = I + uint8(randi([-10 10], size(I)));

[Lw, sw] = lacunarity(I);
[Lw_noise, ~] = lacunarity(I_noise);
[Lb, sb] = lacunarity(I, 'box');
[Lb_noise, ~] = lacunarity(I_noise, 'box');

subplot(121);hold on;
plot(log10(sw), Lw, 'LineWidth', 1.5);
plot(log10(sw), Lw_noise, 'LineWidth', 1.5);
xlabel('Log_{10}(s)');
ylabel('L');
grid on;
title('Window method');
legend({'Original','Noise'});

subplot(122);hold on;
plot(log10(sb), Lb, 'LineWidth', 1.5);
plot(log10(sb), Lb_noise, 'LineWidth', 1.5);
xlabel('Log_{10}(s)');
ylabel('L');
grid on;
title('Box method');
legend({'Original','Noise'});
%% LBP feature extraction
close all;clc;clearvars;


%% Wavelet feature extraction
close all;clc;clearvars;

