function X = LBP_features(I, n_neighbor, features)
%LBP_features Compute features using Local Binary Pattern 
%
%   X = LBP_features(I, n_neighbor, features)
%   Detail explanation goes here
%
%   Input arguments:
%  
%   'I'              Input grayscale image to be analyzed.
%
%                    Accepted values
%
%                    Default: 
%            
%  
%  
%   Output arguments:
%  
%   'ARG1'           Description of the argument. Type and purpose.          
%  
%
%   
%   Notes
%   -----
%   Important usage informationAnother name for a gray-level co-occurrence matrix is a gray-level
%   spatial dependence matrix.
%
%
%   References
%   ----------
%   [1] 
%
%   Example 1
%   ---------      
%   % Example description
%
%     I = imread('cameraman.tif');
%     X = LBP_features(I,'NumLevels',9,'G',[])
%     
%
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2021


if nargin == 1
    n_neighbor = 8;
end
if nargin <= 2
    features = 'low_order';
end 

LBP = extractLBPFeatures(I, 'NumNeighbors', n_neighbor);
% extra parameters can be entered (radius, rotation-invariance)

switch features
    case 'low_order'
        X.mean = mean(LBP);
        X.median = median(LBP);
        X.std = std(LBP);
        X.iqr = iqr(LBP);
        X.range = max(LBP) - min(LBP);
        X.skewness = skewness(LBP);
        X.kurtosis = kurtosis(LBP);
    case 'histogram'
        X = LBP;
end
