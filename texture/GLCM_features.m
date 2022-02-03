function X = GLCM_features(GLCM, varargin)
%GLCM_FEATURES Compute multiple metrics from a GLCM matrix
%
%   X = GLCM_features(M)
%   Compute several features from GLCM matrix
%
%   Input arguments (required):
%  
%   'GLCM'           2D GLCM matrix obtained by graycomatrix()
%  
%   Input arguments (optional):
%
%   'features'       Cell array of strings with the lisst of features to be 
%                    computed. By default the full set of features will be
%                    computed.
%
%
%   Output arguments:
%  
%   'X'              Structure with computed features.          
%  
%
%   
%   Notes
%   -----
%   Some features assume GLCM matrix to be symmetric
%   The naming convention of some parameters is messy in the literature. The 
%   following issues are important:
%   - Homogeneity: inverse difference and inverse difference moment.
%   - Dissimilarity: equal to the difference average
%   - Cluster tendency is equal to sum variance 
%   - Intertia is equal to contrast
%   
%   Additional features not yet implemented: maximal correlation coefficient
%
%   References
%   ----------
%   [1] Haralick R. Shanmugam K. and Dinstein I. "Textural Features for Image
%   Classification", IEEE Transactions on Systems, Man and Cybernetics, 1973
%   http://haralick.org/journals/TexturalFeaturesHaralickShanmugamDinstein.pdf#page10
%
%   [2] Soh L. "Texture Analysis of SAR Sea Ice Imagery Using Gray Level
%   Co-Ocurrence Matrices", IEEE Transations on Gegoscience and Remote Sensing,
%   1999. https://doi.org/10.1109/36.752194
%
%   [3] Zwanenburg et al. The Image Biomarker Standardization Initiative: 
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
                'dif_variance',...
                'dif_entropy',...
                'dissimilarity',...
                'energy', ...
                'entropy',...
                'IMC1',...
                'IMC2',...
                'ID',...
                'IDN',...
                'IDM',...
                'IDMN',...
                'joint_average',...
                'joint_variance',...
                'max_prob',...
                'sum_of_squares',...
                'sum_average',...
%                 'sum_variance',...
                'sum_entropy'};
    

[GLCM, features] = parse_inputs(GLCM, varargin, feature_list);

% Get GLCM matrix size and indexes
N = size(GLCM, 1);
[i, j] = meshgrid(1:N, 1:N);

% Normalize GLCM matrix as a bidimensional probability distribution
Pij = GLCM./sum(GLCM(:)); 

% Partial probability distributions (rows and columns)
Pi = sum(Pij, 1);
Pj = sum(Pij, 2)';
mui = Pi * (1:N)'; % equal to sum(Pij(:) .* i(:))
muj = Pj * (1:N)'; % equal to sum(Pij(:) .* j(:))
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

    Hx = entropy_safe(Pi);  
    Hy = entropy_safe(Pj); % Hy is equal to Hx (symmetry)
    Hxy = entropy_safe(Pij);
    Hxy2 = entropy_safe(P_prod); 
    Hxy1 = Hxy2;
%   Hxy1 = -sum(Pij(:) .* log(P_prod(:))); % equal to Hxy2, no need to compute
end

% Compute features
n_feat = length(features);
for i_feat=1:n_feat
    feature = features{i_feat};
    switch feature
        case 'autocorrelation' 
            % See [2]
            X.autocorrelation = sum(Pij(:).*i(:).*j(:));
            
        case 'cluster_prominence'
            % See [2]
            X.cluster_prominence = sum(Pij(:) .* (i(:)+j(:)-mui-muj).^4); 
            
        case 'cluster_shade'
            % See [2]
            X.cluster_shade = sum(Pij(:) .* (i(:)+j(:)-mui-muj).^3);            
            
        case 'contrast'
            % See [1]
            X.contrast = sum((i(:) - j(:)).^2 .* Pij(:)); 
            
        case 'correlation'
            % See [1]
            X.correlation = (sum(Pij(:).* i(:).*j(:)) - mui*muj)/(sdi*sdj); 
        
        case 'dif_variance'
            % See [1]
            X.dif_average = sum(Pij(:) .* abs(i(:) - j(:))); % dissimilarity
            X.dif_variance = sum(P_dif .* (d-X.dif_average).^2);
        
        case 'dif_entropy'
            % See [1]
            X.dif_entropy = entropy_safe(P_dif);  % See [1].
        
        case 'dissimilarity'
            % See [2]. Equal to the difference average: sum(P_dif.*d)
            X.dissimilarity = sum(Pij(:) .* abs(i(:) - j(:))); 
        
        case 'energy'            
            % See [1] (named Angular Second moment)
            X.energy = sum(Pij(:).^2); 
        
        case 'entropy'
            % See [1]
            X.entropy = entropy_safe(Pij); 
          
        case 'IMC1'
            % See [1]
            if isequal([Hx Hy],[0 0])
                X.IMC1 = nan; % prevent division by 0
                warning("Division by 0 encountered. IMC1 not computed.");
            else
                X.IMC1 = (Hxy - Hxy1)/max([Hx Hy]);
            end
            
        case 'IMC2'            
            % See [1]  
            X.IMC2 = sqrt(1 - exp(-2*(Hxy2 - Hxy))); 
            if ~isreal(X.IMC2)
                % Prevent it from returning complex values
                warning("IMC2 returned a complex value");
                X.IMC2 = nan;
            end
            
        case 'ID'
            % Inverse Difference. See [3]
            X.inverse_difference = sum(Pij(:) ./ (1 +abs(i(:)-j(:))));
        
        case 'IDN'
            % Inverse Difference Normalized
            X.inverse_difference_norm = sum(Pij(:) ./ (1 + abs(i(:)-j(:))/N));
        
        case 'IDM'
            % Inverse Different Moment. See [1].
            X.inverse_difference_moment = sum(Pij(:) ./ (1 + (i(:)-j(:)).^2));           
        
        case 'IDMN'
            % Inverse difference Moment Normalized. See [3]
            X.inverse_difference_moment_norm = sum(Pij(:) ./ (1 + (i(:)-j(:)).^2 ./ N^2));           
        
        case 'joint_average'
            % See [3]. Use it only when GLCM is symmetrical
            X.joint_average = sum(Pij(:) .* i(:) .* j(:));
        
        case 'joint_variance'
            joint_average = sum(Pij(:) .* i(:) .* j(:));
            X.joint_variance = sum(Pij(:) .* (i(:) - joint_average).^2);             
        
        case 'max_prob'
            % See [2]
            X.max_prob = max(Pij(:));
        
        case 'sum_of_squares'
            % See [1]. It is actually a variance
            % As per [3] it should be only used with symmetrical GLCM
            X.sum_of_squares = sum(Pij(:) .* (i(:) - mui).^2); 
        
        case 'sum_average'
            % See [1]
            X.sum_average = sum(P_sum .* k); 
        
% Deprecated as it is the same as sum of squares
%         case 'sum_variance'
            % See [1]
%             X.sum_variance = sum(P_sum .* (k - X.sum_average).^2);
        
        case 'sum_entropy'
            % See [1]
            X.sum_entropy = entropy_safe(P_sum);
        
        otherwise
            warning(['Unknown feature:' feature]);
    end
end
end

function H = entropy_safe(p)
% Function to compute entropies in a safe way
if abs(sum(p(:)) - 1) > 1e-4
    warning('Probabilities do not add to 1. Unable to compute entropy');
    H = nan;
else
    p = p(p~=0); % avoid log(0) = -Inf
    H = -sum(p .* log2(p));
end
end

function [GLCM, features] = parse_inputs(GLCM, extra_args, feature_list)

if length(extra_args)> 2
    error("Function expects a maximum of 2 input arguments"); 
end

% Check GLCM
if ~isnumeric(GLCM)
    error(["GLCM matrix is expected to be numeric but is " class(GLCM)]); 
end
has_nan = sum(isnan(GLCM(:))) > 0;
has_neg = sum(GLCM(:) < 0);
has_inf = sum(isinf(GLCM(:))) > 0;
if has_nan | has_neg | has_inf
    error("GLCM cannot have NaN, Inf or negative values");
end

% Check features
if isempty(extra_args)
    features = feature_list;
else
    features = extra_args{1};    
    if ischar(features)
        return;
    elseif iscell(features)
        cell_is_char = cellfun(@(x) class(x), a, 'UniformOutput', false);
        if any(~cell_is_char)
            error("At least one of the features is not of type char"); 
        end
    else
        error(["features is expected to be char or cell array of chars but is ",...
            class(features)]);
    end
end
end