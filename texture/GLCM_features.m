function X = GLCM_features(GLCM, features)
%GLCM_FEATURES Compute multiple metrics from a GLCM matrix
%
%   X = GLCM_features(M)
%   Detail explanation goes here
%
%   Input arguments:
%  
%   'GLCM'           2D GLCM matrix obtained by graycomatrix()
%  
%   'features'       Cell array of strings with the lisst of features to be 
%                    computed. By default the full set of features will be
%                    computed.
%
%
%   Output arguments:
%  
%   'X'              Structure with computed metrics.          
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
%   [2] Soh L. "Texture Analysis of SAR Sea Ice Imagery Using Gray Level
%   Co-Ocurrence Matrices", IEEE Transations on Gegoscience and Remote Sensing,
%   1999
%   https://doi.org/10.1109/36.752194

%   [3] Zwanenburg et al.The Image Biomarker Standardization Initiative: 
%   Standardized Quantitative Radiomics for High-Throughput Image-based 
%   Phenotyping, 2020,  Radiology, 
%   https://pubs.rsna.org/doi/full/10.1148/radiol.2020191145
%   Implementation: https://pyradiomics.readthedocs.io/en/latest/features.html?highlight=imc#radiomics.glcm.RadiomicsGLCM.getImc1FeatureValue
%
%   Example
%   ---------      
%   % Metrics computation
%
%     GLCMS = graycomatrix(I,'NumLevels',9,'G',[])
%     X = GLCM_features(M)
%     
%
%  
%   David Romero-Bascones, Biomedical Engineering Department, Mondragon
%   Unibertsitatea, 2021
%   dromero@mondragon.edu

feature_list = {'autocorrelation',...
                'cluster_prominence',...
                'cluster_shade',...
                'contrast',...
                'correlation',...
                'dif_average',...
                'dif_variance',...
                'dif_entropy',...
                'energy', ...
                'entropy',...
                'homogeneity',...
                'IMC1',...
                'IMC2',...
                'inverse_difference',...
                'inverse_difference_norm',...
                'inverse_difference_moment',...
                'inverse_difference_moment_norm',...
                'joint_average',...
                'joint_variance',...
                'max_prob',...
                'sum_of_squares',...
                'sum_average',...
                'sum_variance',...
                'sum_entropy'};
    
if nargin==1
    features = feature_list;
end

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
Pj = sum(Pij, 2)';
mui = Pi * (1:N)';
muj = Pj * (1:N)';
sdi = sqrt(sum(Pij(:) .* (i(:) - mui).^2));
sdj = sqrt(sum(Pij(:) .* (j(:) - muj).^2));

% Compute the distribution of i+j (sum of indexes) if necessary
if any(ismember(features, {'sum_average', 'sum_variance', 'sum_entropy'}))
    P_sum = zeros(1, 2*N-1);
    k = 2:2*N;
    for ik=1:length(k)
        mask = double((i + j) == k(ik));
        P_sum(ik) = sum(Pij(:).*mask(:));
    end
end

% Compute the distribution of i-j (difference of indexes) if necessary
if any(ismember(features, {'dif_average', 'dif_variance', 'dif_entropy'}))
    P_dif = zeros(1, N);
    d = 0:N-1;
    for id=1:length(d)
        mask = double(abs(i - j) == d(id));
        P_dif(id) = sum(Pij(:).*mask(:));
    end
end

% Compute entropies for information measures of correlation
if any(ismember(features, {'IMC1', 'IMC2'}))
    [Pi_mat, Pj_mat] = meshgrid(Pi, Pj); % in matrix form for latter calculations
    P_prod = Pi_mat .* Pj_mat; % Product matrix

    Hx = my_entropy(Pi);  % Hy is equal to Hx (symmetry)
    Hxy = my_entropy(Pij);
    Hxy2 = my_entropy(P_prod); 
    Hxy1 = Hxy2;
%     Hxy1 = -sum(Pij(:) .* log(P_prod(:))); % equal to Hxy2, no need to
%     compute    
end

% Compute features
n_feat = length(features);
for i_feat=1:n_feat
    feature = features{i_feat};
    switch feature
        case 'autocorrelation'
            X.autocorrelation = sum(Pij(:).*i(:).*j(:));
        case 'cluster_prominence'
            X.cluster_prominence = sum(Pij(:) .* (i(:)+j(:)-2*mui).^4);
        case 'cluster_shade'
            X.cluster_shade = sum(Pij(:) .* (i(:)+j(:)-2*mui).^3);
        case 'contrast'
            X.contrast = sum((i(:) - j(:)).^2 .* Pij(:));  % See [1].
        case 'correlation'
            % See [1]
            X.correlation = (sum(Pij(:).* i(:).*j(:)) - mui*muj)/(sdi*sdj); 
        case 'dif_average'
            % Equal to dissimilarity = sum(Pij(:) .* abs(i(:) - j(:))); 
            X.dif_average = sum(P_dif.*d);
        case 'dif_variance'
            X.dif_variance = sum(P_dif .* (d-X.dif_average).^2);  % See [1].
        case 'dif_entropy'
            X.dif_entropy = my_entropy(P_dif);  % See [1].
        case 'energy'            
            X.energy = sum(Pij(:).^2); % Named Angular Second Moment in [1]
        case 'entropy'
            X.entropy = my_entropy(Pij);  % See [1].
        case 'homogeneity'
            X.homogeneity = graycoprops(GLCM, 'homogeneity').Homogeneity;                
        case 'IMC1'
             X.IMC1 = (Hxy - Hxy1)/Hx;            
        case 'IMC2'            
            X.IMC2 = sqrt(1 - exp(-2*(Hxy2 - Hxy)));            
            if ~isreal(X.IMC2)
                % Prevent it from returning complex values
                warning("Information measure of correlation 2 returned a complex value");
                X.IMC2 = 0;
            end
        case 'inverse_difference'
            X.inverse_difference = sum(Pij(:) ./ (1 +abs(i(:)-j(:))));
        case 'inverse_difference_norm'
            X.inverse_difference_norm = sum(Pij(:) ./ (1 + abs(i(:)-j(:))/N));
        case 'inverse_difference_moment'
            X.inverse_difference_moment = sum(Pij(:) ./ (1 + (i(:)-j(:)).^2)); % See [1].            
        case 'inverse_difference_moment_norm'
            X.inverse_difference_moment_norm = sum(Pij(:) ./ (1 + (i(:)-j(:)).^2 ./ N^2));           
        case 'joint_average'
            X.joint_average = sum(Pij(:) .* i(:) .* j(:));
        case 'joint_variance'
            joint_average = sum(Pij(:) .* i(:) .* j(:));
            X.joint_variance = sum(Pij(:) .* (i(:) - joint_average).^2);             
        case 'max_prob'
             X.max_prob = max(Pij(:));
        case 'sum_of_squares'
            % It is actually a variance
            X.sum_of_squares = sum(Pij(:) .* (i(:) - mui).^2); % See [1].
        case 'sum_average'
            X.sum_average = sum(P_sum .* k);  % See [1].
        case 'sum_variance'
            X.sum_variance = sum(P_sum .* (k - X.sum_average).^2);  % See [1].
        case 'sum_entropy'
            X.sum_entropy = my_entropy(P_sum); % See [1].    
        otherwise
            warning(['Unknown feature:' feature]);
    end
end
end

function H = my_entropy(p)
if abs(sum(p(:)) - 1) > 1e-4
    warning('Probabilities do not add to 1. Unable to compute entropy');
    H = nan;
else
    p = p(p~=0); % avoid log(0) = -Inf
    H = -sum(p .* log2(p));
end
end