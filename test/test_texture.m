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
close all;clc;
I = imread('../data/D33.gif');
I_noise = double(I) + rand(size(I));

methods = {'local_bin', 'box_3d', 'window'};
r = 2:2:20;

for i=1:3
    [L, r] = lacunarity(I,'box_3d',r);
    [Ln, ~] = lacunarity(I_noise,'box_3d',r);

    subplot(1,3,i);hold on;
    plot(r, L, 'LineWidth', 1.5);
    plot(r, Ln, 'LineWidth', 1.5);
    xlabel('Log(r)');
    ylabel('L');
    grid on;
    title(methods{i});
    legend({'Original','Noise'});

end

%% LBP feature extraction
close all;clc;clearvars;


%% Wavelet feature extraction
close all;clc;clearvars;

