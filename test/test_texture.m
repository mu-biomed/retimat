close all;clc;clearvars;
addpath(genpath('..'));

%% Fractal Dimension
clc;close all;
% I = imread('../data/D33.gif');
I = imread('cameraman.tif');
I_noise1 = uint8(double(I) + 10*randn(size(I)));
I_noise2 = double(I) + 100*randn(size(I));

FD_dbc = fractal_dimension(I, 'DBC', true);hold on;
FD_irdbc = fractal_dimension(I, 'IR_DBC', true);hold on;
FD_irdbc_noise = fractal_dimension(I_noise1, 'IR_DBC', true);
FD_irdbc_noise2 = fractal_dimension(I_noise2, 'IR_DBC', true);
legend({'DBC', 'IRDBC', 'IRDBC \sigma_{1}', 'IRDBC \sigma_{2}'},...
    'Location','southeast','FontSize',12);
%% GLCM features
close all;clc;clearvars;

I = randn(200, 200);
GLCM = graycomatrix(I,'NumLevels',250);
X = compute_glcm_metrics(GLCM);

%% Lacunarity
close all;clc;clearvars;
I = imread('../data/D33.gif');
% I = ones(257, 256);
% I = imread('cameraman.tif');
I_noise = double(I) + rand(size(I));

[Lw, sw] = lacunarity(I,'box_3d');
[Lw_noise, ~] = lacunarity(I,'window');
% [Lb, sb] = lacunarity(I, 'box');
% [Lb_noise, ~] = lacunarity(I_noise, 'box');

subplot(121);hold on;
plot(log10(sw), Lw, 'LineWidth', 1.5);
plot(log10(sw), Lw_noise, 'LineWidth', 1.5);
xlabel('Log_{10}(s)');
ylabel('L');
grid on;
title('Window method');
legend({'Original','Noise'});

% subplot(122);hold on;
% plot(log10(sb), Lb, 'LineWidth', 1.5);
% plot(log10(sb), Lb_noise, 'LineWidth', 1.5);
% xlabel('Log_{10}(s)');
% ylabel('L');
% grid on;
% title('Box method');
% legend({'Original','Noise'});
%% LBP feature extraction
close all;clc;clearvars;


%% Wavelet feature extraction
close all;clc;clearvars;

