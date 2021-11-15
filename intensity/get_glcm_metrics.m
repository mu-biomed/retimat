function X = get_glcm_metrics(M)

% Definitions from:
%
% Zwanenburg et al.The Image Biomarker Standardization Initiative: Standardized
% Quantitative Radiomics for High-Throughput Image-based Phenotyping, 2020,
% Radiology, https://pubs.rsna.org/doi/full/10.1148/radiol.2020191145
% Implementation: https://pyradiomics.readthedocs.io/en/latest/features.html?highlight=imc#radiomics.glcm.RadiomicsGLCM.getImc1FeatureValue


% Tazarjani et al. Retinal OCT Texture Analysis for Differentiating Healthy
% Controls from Multiple Sclerosis (MS) with/without Optic Neuritis, BioMed
% Research International, 2021

% Original paper with definitions:
% Haralick, Textural Features for Image Classification, 1973
% http://haralick.org/journals/TexturalFeaturesHaralickShanmugamDinstein.pdf#page10


% P: GLCM (2D probability density function)

% Get GLCM matrix size and indexes
N = size(M, 1);
[i, j] = meshgrid(1:N, 1:N); % Indexes

% Normalize GLCM matrix as a bidimensional probability distribution
Pij = M./sum(M(:)); 

% Partial probability distributions (rows and columns)
Pi = sum(Pij, 1);
Py = sum(Pij, 2)';
mui = Pi * (1:N)';
muj = Py * (1:N)';
sdx = sqrt(sum(Pij(:) .* (i(:) - mui).^2));
sdy = sqrt(sum(Pij(:) .* (j(:) - muj).^2));
[Px_mat, Py_mat] = meshgrid(Pi, Py); % in matrix form for latter calculations
Pxy = Px_mat .* Py_mat; % Product matrix

% Distribution of i+j (sum of indexes)
P_sum = zeros(1, 2*N-1);
k = 2:2*N;
for ik=1:length(k)
    mask = double((i + j) == k(ik));
    P_sum(ik) = sum(Pij(:).*mask(:));
end

% Distribution of i-j (difference of indexes)
P_dif = zeros(1, N);
d = 0:N-1;
for id=1:length(d)
    mask = double(abs(i - j) == d(id));
    P_dif(id) = sum(Pij(:).*mask(:));
end


% ------------------------ Joint Features ------------------------------------
X.max_prob = max(Pij(:));
X.joint_average = sum(Pij(:) .* i(:) .* j(:));
X.joint_variance = sum(Pij(:) .* (i(:) - X.joint_average).^2);

X.energy = sum(Pij(:).^2);
X.entropy = -sum(Pij(:).*log(Pij(:)));
X.contrast = sum((i(:) - j(:)).^2 .* Pij(:)); 

X.correlation = sum(Pij(:) .* (i(:) - mui)./sdx .* (j(:) - muj)./sdy); 
X.autocorrelation = sum(Pij(:) .* i(:) .* j(:));

X.sum_of_squares = sum(Pij(:) .* (i(:) - mui).^2);
X.cluster_shade = sum(Pij(:) .* (i(:) + j(:) -2*mui).^3);
X.cluster_prominence = sum(Pij(:) .* (i(:) + j(:) -2*mui).^4);
% X.dissimilarity = sum(Pij(:) .* abs(i(:) - j(:))); % equal to difference average

X.inverse_difference = sum(Pij(:) ./ (1 + abs(i(:) - j(:))));
X.inverse_difference_norm = sum(Pij(:) ./ (1 + abs(i(:) - j(:))/N));
X.inverse_difference_moment = sum(Pij(:) ./ (1 + (i(:) - j(:)).^2));% also called homogeneity
X.inverse_difference_moment_norm = sum(Pij(:) ./ (1 + (i(:) - j(:)).^2 ./ N^2));
X.homogeneity = sum(Pij(:) ./ (1 + (i(:) - j(:)).^2));  % Likely to be wrong on Tazarjani,2021 (same dfinition as inverse moment)
X.homogeneity2 = graycoprops(M, 'homogeneity').Homogeneity; % different definition (abs vs squared)

% inverse_variance pending

Hx = -sum(Pi.*log(Pi));  % Hy is equal to Hx (symmetry)
Hxy = X.entropy;

Hxy1 = -sum(Pij(:) .* log(Pxy(:)));
Hxy2 = -sum(Pxy(:) .* log(Pxy(:))); % equal to Hxy1

X.IMC1 = (Hxy - Hxy1)/Hx;
X.IMC2 = sqrt(1 - exp(-2*(Hxy2 - Hxy)));

% Prevent it returning complex values
if ~isreal(X.IMC2)
    warning("Information measure of correlation 2 returned a complex value");
    X.IMC2 = 0;
end


% ------------------------ Sum Features ------------------------------------
X.sum_average = sum(P_sum .* k);
X.sum_variance = sum(P_sum .* (k - X.sum_average).^2);
X.sum_entropy = -sum(P_sum .* log(P_sum));



% ----------------------- Difference features -------------------------------
X.dif_average = sum(P_dif .* d);
X.dif_variance = sum(P_dif .* (d - X.dif_average).^2);
X.dif_entropy = -sum(P_dif .* log(P_dif));
