function [L, r] = lacunarity(I, method, visu, r)
%lacunarity Compute lacunarity of a gray-scale image
%
%   [L, r] = lacunarity(I, method)
%   Compute lacunarity of a grayscale/binary image for a set of window sizes r.
%
%   Input arguments:
%  
%   'I'              Input grayscale/binary image. 
%             
%   'method          Method to compute lacunarity.          
%                    Options: ['window', 'box']
%                    Default: 'window'
%
%   'visu'           Optional. If true a plot with the computed values is
%                    shown.
%                    Default = false
%
%   'r'              1D array with window size values. By default 10 values 
%                    from 1 to the minimum between rows and columns are used.
%
%
%   Output arguments:
%  
%   'L'              1D array with lacunarity values.          
%  
%   'r'              1D array with window size values.          
%  
%
%   Notes
%   -----
%   Color images will be converted to grayscale. By default window sizes from 1
%   to min{M, N} are used.
%
%
%   References
%   ----------
%   [1] Manikka-Baduge D.C. and Dougherty G. "Texture analysis using lacunarity 
%   and average local variance", Proceedings of SPIE, 2009 
%   https://www.researchgate.net/publication/252776086_Texture_analysis_using_lacunarity_and_average_local_variance
%   
%   [1] Allain C. and Clitre M. "Charaterizing the lacunarity of random and
%   deterministic fractal sets", Physical Review, 1991.
%   https://doi.org/10.1103/PhysRevA.44.3552
%
%   [2] Roy A. and Perfect E. "Lacunarity Analyses of Multifractal and Natural
%   Grayscale Patterns", Fractals, 2014
%   https://doi.org/10.1142/S0218348X14400039
%
%   [4]"A new approach to estimate lacunarity of texture images"
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
elseif nargin == 2
    visu = false;
end

if ~isnumeric(I) & ~islogical(I)
    error(['Input image must be numeric (grayscale) or logical (binary)', ...
        ' but is: ', class(I)]);
end

if ndims(I) == 3
    warning('Converting input image to grayscale');
    I = rgb2gray(I);
end

if isequal(class(I), 'uint8')
    n_lev = 2^8;
elseif isequal(class(I), 'uint16')
    n_lev = 2^16;
elseif isequal(class(I), 'logical')
    n_lev = 2;
else
    warning(['Input image is of type ', class(I), '. Converting to uint8.']);
    I = uint8(255*(I - min(I(:)))./(max(I(:)) - min(I(:))));
end

[M, N] = size(I);

r_max = min([M N]);  %  maximum window size equals the minimum size
% r_vals = [2.^(0:(log2(r_max) - 1)) r_max];  %  powers of two
r = round(exp(linspace(log(1), log(r_max), 10)));  %  log-uniform points

n_scale = length(r);
L = nan(1, n_scale);

for i_scale=1:n_scale
    r_i = r(i_scale);
    
    switch method
        case 'local_bin'  %  See reference [4] Section 3.1. Slow because of the
            % nested for loops. Always 1 for r=1.
            n_row = M - r_i + 1;  %  number of rows with starting window
            n_col = N - r_i + 1;  %  number of columns with starting window
            m = nan(n_row, n_col);
            for i=1:n_row  
                for j=1:n_row
                    I_box = I(i:i+r_i-1, j:j+r_i-1);  % window
                    m(i,j) = sum(I_box(:) >= mean(I_box(:)));
                end
            end
            
        case 'box_3d'  %  See reference [4] Section 3.2
            n_lev = max(I(:));  % L in the paper
            
            L(i_scale) = lacunarity_3d(I, M, n_lev, 2);
            
        case 'window'  %  See reference [2] Section 2.2
            % Sliding-window through all the image. For each position we
            % compute the sum of the values in each window as mass.
            m = reshape(conv2(I, ones(r_i, r_i), 'valid'),[], 1);
            
%         case 'normalized'
            
        otherwise
            error("Provided method is not supported. Use 'box' or 'window'.");
    end
    L(i_scale) = var(m(:))/mean(m(:))^2 + 1;
end

if visu
    plot(log(r), log(L), '-o', 'LineWidth', 1.5);
    xlabel('Log(s)');
    ylabel('Log(L)');
    grid on;
end
end

function Lac = lacunarity_3d(I, M, L, r)
    S_max = r^3;
    l = L - r + 1;  % number of boxes in Y direction

    vl = zeros(1,l);
    
    for x=0:r-1
        vl = compute_column(I, vl, l, x, 0, r, 1);
    end
    
    vl1 = vl;
    
    n = zeros(1, S_max);
    
    first = true;
    for y=(r-1):(M-1)
        for x=(r-1):(M-1)
            if first
                first = false;
            else
                vl = compute_column(I, vl, x-r, y-r+1, r, -1);
                vl = compute_column(I, vl, x, y-r+1, r, 1);                
            end
            for i=1:(l-1)
                if vl(i)~=0 & vl(i+1)~=S_max
                    n(vl(i)) = n(vl(i)) + 1;  % add boxes with mass vl 
                end
            end
        end
        if y < M -1
            vl1 = compute_row(I, vl1, l, 0, y-r+1, r, -1);
            vl1 = compute_row(I, vl1, l, 0, y+1, r, 1);
            vl = vl1;
        end
        first = true;
    end
    N = 0;
    Z1 = 0;
    Z2 = 0;
    
    for S=0:S_max
        N = N + n(S); 
    end
    for S=0:S_max
        i = S; % ?
        Z1 = Z2 + i * n(S)/N;
        Z2 = Z2 + i^2 * n(S)/N;
    end
    if Z1 == 0
        Lac = 0;
    else
        Lac = Z2/Z1^2; 
    end
end

function vl = compute_column(I, vl, l, x, y, r, t)
    for y1 = y:y+r-1
        [cc, cr] = count_boxes(l, I(x+1,y1+1), r);
        
        vl = vl + t*r;

        i = cc;
        while (i<l) & (cr > 0)
            vl(i) = vl(i) + t * cr;
            cr = cr - 1;
            i = i + 1;
        end
    end
end

function compute_row(I, vl, l, x, y, r, t)
    for x1 = x:x+r-1
        [cc, cr] = count_boxes(l, I(x1+1,y+1), r);

        vl = vl + t * r;
        
        i = cc;
        while (i<l) & (cr > 0)
            vl(i) = vl(i) + t * cr;
            cr = cr - 1;
            i = i + 1;
        end
    end
end

function [cc, cr] = count_boxes(l, h, r)
    % Number of cubes which a given pixel intercepts in function of its gray
    % level intensity;
    % h: pixel intensity (height)
    % cc: crossed cubes (intercepted)
    % l: maximal number of cumbes allowed in the pixel intensity direction
    % cr: cros
    cc = h - r + 1;
    
    if cc <= 0
        cc = 0;
        cr = h;
    else
        if cc < l
            cc = r - 1;
            cr = 0;
        end
    end
end