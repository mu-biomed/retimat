function X = compute_glcm_metrics(GLCM)
%compute_glcm_metrics Compute multiple metrics from a GLCM matrix
%
%   X = compute_glcm_metrics(M)
%   Detail explanation goes here
%
%   Input arguments:
%  
%   'GLCM'           2D GLCM matrix obtained by graycomatrix()
%  
%  
%   Output arguments:
%  
%   'X'           Structure with computed metrics.          
%  
%
%   
%   Notes
%   -----
%   
%
%   References
%   ----------
%   [1] Haralick R. Shanmugam K. and Dinstein I. "Textural Features for Image
%   Classification", IEEE Transactions on Systems, Man and Cybernetics, 1973
%   http://haralick.org/journals/TexturalFeaturesHaralickShanmugamDinstein.pdf#page10
%
%   [2] Zwanenburg et al.The Image Biomarker Standardization Initiative: 
%   Standardized Quantitative Radiomics for High-Throughput Image-based 
%   Phenotyping, 2020,  Radiology, 
%   https://pubs.rsna.org/doi/full/10.1148/radiol.2020191145
%   Implementation: https://pyradiomics.readthedocs.io/en/latest/features.html?highlight=imc#radiomics.glcm.RadiomicsGLCM.getImc1FeatureValue
%
%   [3] Tazarjani et al. Retinal OCT Texture Analysis for Differentiating
%   Healthy Controls from Multiple Sclerosis (MS) with/without Optic Neuritis, 
%   BioMed Research International, 2021
%
%   Example
%   ---------      
%   % Metrics computation
%
%     GLCMS = graycomatrix(I,'NumLevels',9,'G',[])
%     X = get_glcm_metrics(M)
%     
%
%  
%   David Romero-Bascones, Biomedical Engineering Department, Mondragon
%   Unibertsitatea, 2021
%   dromero@mondragon.edu

if ~isnumeric(GLCM)
    error(['Provided GLCM matrix must be numeric but is ' class(GLCM)]); 
end
% Get GLCM matrix size and indexes
N = size(GLCM, 1);
[i, j] = meshgrid(1:N, 1:N);

% Normalize GLCM matrix as a bidimensional probability distribution
Pij = GLCM./sum(GLCM(:)); 

% Partial probability distributions (rows and columns)
Pi = sum(Pij, 1);
Py = sum(Pij, 2)';
mui = Pi * (1:N)';
muj = Py * (1:N)';
sdi = sqrt(sum(Pij(:) .* (i(:) - mui).^2));
sdj = sqrt(sum(Pij(:) .* (j(:) - muj).^2));
[Pi_mat, Pj_mat] = meshgrid(Pi, Py); % in matrix form for latter calculations
P_prod = Pi_mat .* Pj_mat; % Product matrix

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


%------------------------------------------------------------------------------
% Features proposed in Haralick, 1973
%------------------------------------------------------------------------------

% Angular Second Moment or Energy
X.energy = sum(Pij(:).^2); 

% Contrast
X.contrast = d * P_dif'; % Equal to difference mean
X.contrast_2 = sum((i(:) - j(:)).^2 .* Pij(:));  % Different definition in Tazarjani, 2021 and MAtlab itself

% Correlation
X.correlation = sum(Pij(:) .* i(:) .* j(:) - mui*muj)/(sdi*sdj); 
X.correlation_2 = sum(Pij(:) .* (i(:) - mui)./sdi .* (j(:) - muj)./sdj); % Different definition in Tazarjani,2021

% Sum of squares (Variance)
X.sum_of_squares = sum(Pij(:) .* (i(:) - mui).^2);

% Inverse different moment 
X.inverse_difference_moment = sum(Pij(:) ./ (1 + (i(:) - j(:)).^2));
% Called homogeneity in Tazarjani, 2021

% Sum average (Haralick, 1973)
X.sum_average = sum(P_sum .* k);

% Sum variance (Haralick, 1973)
X.sum_variance = sum(P_sum .* (k - X.sum_average).^2);

% Sum entropy (Haralick, 1973)
X.sum_entropy = -sum(P_sum .* log(P_sum));

% Entropy (Haralick, 1973)
X.entropy = -sum(Pij(:).*log(Pij(:))); 

% Difference variance
X.dif_variance = sum(P_dif .* (d - X.dif_average).^2);

% Difference entropy
X.dif_entropy = -sum(P_dif .* log(P_dif));

% Information measures of correlation
Hx = -sum(Pi.*log(Pi));  % Hy is equal to Hx (symmetry)
Hxy = X.entropy;

Hxy1 = -sum(Pij(:) .* log(P_prod(:)));
Hxy2 = -sum(P_prod(:) .* log(P_prod(:))); % equal to Hxy1

X.IMC1 = (Hxy - Hxy1)/Hx;
X.IMC2 = sqrt(1 - exp(-2*(Hxy2 - Hxy)));

if ~isreal(X.IMC2)
    % Prevent it from returning complex values
    warning("Information measure of correlation 2 returned a complex value");
    X.IMC2 = 0;
end

% Maximal Correlation Coefficient
% Complex implementation

% Two different deffinitions found
X.homogeneity = sum(Pij(:) ./ (1 + (i(:) - j(:)).^2));  % Likely to be wrong on Tazarjani,2021 (same dfinition as inverse moment)
X.homogeneity2 = graycoprops(GLCM, 'homogeneity').Homogeneity; % different definition (abs vs squared)


% Maximum probability
X.max_prob = max(Pij(:));

% Joint mean
X.joint_average = sum(Pij(:) .* i(:) .* j(:));

% Joint variance
X.joint_variance = sum(Pij(:) .* (i(:) - X.joint_average).^2);

% Autocorrelation
X.autocorrelation = sum(Pij(:) .* i(:) .* j(:));

% Cluster shade
X.cluster_shade = sum(Pij(:) .* (i(:) + j(:) -2*mui).^3);

% Cluster prominence
X.cluster_prominence = sum(Pij(:) .* (i(:) + j(:) -2*mui).^4);

% Inverse difference
X.inverse_difference = sum(Pij(:) ./ (1 + abs(i(:) - j(:))));

% Normalized inverse difference
X.inverse_difference_norm = sum(Pij(:) ./ (1 + abs(i(:) - j(:))/N));

% Normalized inverse different moment
X.inverse_difference_moment_norm = sum(Pij(:) ./ (1 + (i(:) - j(:)).^2 ./ N^2));


% inverse_variance pending






% ----------------------- Difference features -------------------------------

% Difference average
% Equal to dissimilarity = sum(Pij(:) .* abs(i(:) - j(:))); 
X.dif_average = sum(P_dif .* d);




