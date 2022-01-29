function [X, LBP_hist] = LBP_features(I, n_neighbor, features)
%LBP_features Compute features using Local Binary Pattern 
%
%   X = LBP_features(I, n_neighbor, features)
%   Detail explanation goes here
%
%   Input arguments:
%  
%   'I'              Input grayscale image to be analyzed.
%
%   'n_neighbor'     Number of neighbors. Default is 8
%
%   'features'       String or a cell array of strings defining the nmerical   
%                    features to be computed. By default all the features will   
%                    be returned. If 'none' then no features are computed and 
%                    only the LBP histogram is returned.
%  
%  
%   Output arguments:
%  
%   'X'              Struct with features computed from LBP histogram.
%
%   'H'              LBP histogram
%
%   
%   Notes
%   -----
%   
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
%     X = LBP_features(I, 8)
%     
%
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2021


if nargin == 1
    n_neighbor = 8;
end
if nargin <= 2
    features = {'mean',...
                'median',...
                'std',...
                'var',...
                'iqr',...
                'range',...
                'skewness',...
                'kurtosis',...
                'entropy'};
end

LBP_hist = extractLBPFeatures(I, 'NumNeighbors', n_neighbor);
% extra parameters can be entered (radius, rotation-invariance)

if ischar(features)
    features = {features};
elseif ~iscell(features)
    error('features argument must be either a string or a cell array of strings');
end

for i=1:length(features)
    switch features{i}
        case 'mean'
            X.mean = mean(LBP_hist);
        case 'median'
            X.median = median(LBP_hist);
        case 'std'
            X.std = std(LBP_hist);
        case 'var'
            X.var = var(LBP_hist);
        case 'iqr'
            X.iqr = iqr(LBP_hist);
        case 'range'
            X.range = range(LBP_hist);
        case 'skewness'
            X.skewness = skewness(LBP_hist);
        case 'kurtosis'
            X.kurtosis = kurtosis(LBP_hist);
        case 'entropy'
            p = LBP_hist./sum(LBP_hist);
            X.entropy = -sum(log2(p(p~=0)).*p(p~=0));
        case 'none'
            X = nan;
        otherwise
            error("Unknown feature.");
    end
end
