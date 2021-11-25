%% Test GLCM feature extraction
close all;clc;clearvars;

I = randn(200, 200);
GLCM = graycomatrix(I,'NumLevels',250);
X = compute_glcm_metrics(GLCM);

%% Test LBP feature extraction
close all;clc;clearvars;



%% Text Wavelet feature extraction
close all;clc;clearvars;


