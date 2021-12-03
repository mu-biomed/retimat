function [L, s] = lacunarity(I, method)
%lacunarity Compute lacunarity of a gray-scale image
%
%   [L, s] = lacunarity(I, method)
%   Compute image lacunarity for a set of values s.
%
%   Input arguments:
%  
%   'I'              Input gray-scale image. 
%             
%   'method          Method to compute lacunarity.          
%                    Options: ['window', 'box']
%                    Default: 'window'
%
%
%   Output arguments:
%  
%   'L'              1D array with lacunarity values.          
%  
%   's'              1D array with lacunarity values.          
%  
%
%   Notes
%   -----
%   The image is assumed to be gray-scale. The box/window size is fixed to the
%   range s = [1 2 4 8 16 32 64 85 128 256]. This might result in errors or
%   incomplete results when working with images smaller and bigger than 256x256
%   pixels.
%
%
%   References
%   ----------
%   [1] Manikka-Baduge D.C. and Dougherty G. "Texture analysis using lacunarity 
%   and average local variance", Proceedings of SPIE, 2009 
%   https://www.researchgate.net/publication/252776086_Texture_analysis_using_lacunarity_and_average_local_variance
%   
%   [2] Roy A. and Perfect E. 2 Lacunarity Analyses of Multifractal and Natural
%   Grayscale Patterns, Fractals, 2014
%   http://dx.doi.org/10.1142/S0218348X14400039
%
%
%   Example
%   ---------      
%   % Compute lacunarity
%
%     I = imread('cameraman.tif');
%     [L, s] = lacunarity(I)
%     
%
%  
%   David Romero-Bascones dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2021


if nargin == 1
    method = 'window'; 
end

% Adapt input image if necessary
if ndims(I) == 3
    warning("Converting input image to grayscale");
    I = rgb2gray(I);
end

[n_row, n_col] = size(I);

% Define necessary variables
I_max = 256; % gray-scale image
s = [1 2 4 8 16 32 64 85 128 256];
n_scale = length(s);
L = nan(1, n_scale);

for i_scale=1:n_scale
    s_i = s(i_scale);
    
    switch method
        case 'box'
            n_level = floor(I_max/s_i);
            n_box = (n_row - s_i + 1)*(n_col - s_i + 1);
            m = nan(n_level, n_box);
            for i_level=1:n_level
                
                I_lev = (I + 1) - s_i*(i_level - 1);
                I_lev(I_lev > s_i) = s_i;
                I_lev(I_lev < 0) = 0;
                
                m(i_level, :) = reshape(conv2(I_lev', ones(s_i, s_i), 'valid'),[], 1);
            end
            
        case 'window'
                m = reshape(conv2(I, ones(s_i, s_i), 'valid'),[], 1);
        
        otherwise
            error(["Method: " method " not supported"]);
    end
    L(i_scale) = var(m(:))/mean(m(:))^2 + 1;
end
end
