function [L, r] = lacunarity(I, method, r, visu)
%lacunarity Compute lacunarity of a gray-scale image
%
%   [L, r] = lacunarity(I, method, r, visu)
%   Compute lacunarity of a grayscale/binary image for a set of window sizes r.
%   Lacunarity is a measure of how an image covers the space. Images/masks with
%   bigger gaps will show a higher lacunarity values. Instead of a single value
%   lacunarity is often computed for a set of scales defined by the size of the
%   box used to compute it. It is usually considered as a complementary metric 
%   to fractal dimension.
%
%
%   Input arguments:
%  
%   'I'              Input grayscale/binary image. 
%             
%   'method'         Method to compute lacunarity.          
%                    Options:
%                     'window': Fast and simple. See [2] Section 2.2.
%                     'local_bin': See [3] Section 3.1.
%                     'box_3d': Slow. See [3] Section 3.2.
%                    Default: 'window'
%
%   'r'              Optional. 1D array with window size values. By default 10 
%                    values from 1 to the minimum between rows and columns 
%                    (min{M, N}) are used.
%
%   'visu'           Optional. If true a plot with the computed values is
%                    shown.
%                    Default = false
%
%
%
%   Output arguments:
%  
%   'L'              1D array with lacunarity values.          
%  
%   'r'              1D array with box size values.          
%  
%
%   Notes
%   -----
%   Color images will be converted to grayscale. There are multiple definitions
%   of lacunarity, with different calculation methods. This function covers 3
%   method with different phylosophies. The scale of each method is not
%   directly comparable.
%   There is also a "normalized Lacunarity" proposed in [2] that may be
%   interesting to implement in the future.
%
%
%   References
%   ----------
%   [1] Allain C. and Clitre M. "Charaterizing the lacunarity of random and
%   deterministic fractal sets", Physical Review, 1991.
%   https://doi.org/10.1103/PhysRevA.44.3552
%
%   [2] Roy A. and Perfect E. "Lacunarity Analyses of Multifractal and Natural
%   Grayscale Patterns", Fractals, 2014
%   https://doi.org/10.1142/S0218348X14400039
%
%   [4] Backes A.R. "A new approach to estimate lacunarity of texture images",
%   Pattern Recognition Letters, 2013
%
%
%   Example
%   ---------      
%   % Compute lacunarity
%
%     I = imread('cameraman.tif');
%     [L, s] = lacunarity(I)
%     
%   % Compute Lacunarity choosing a method and custom r values
%     
%     I = imread('cameraman.tif');
%     r = 2:10;
%     [L, ~] = lacunarity(I, 'local_bin', r);
%
%  
%   David Romero-Bascones dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2021


if nargin == 1
    method = 'window'; 
elseif nargin == 2
    visu = false;
elseif nargin == 3
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

if ~isequal(class(I), 'uint8') | isequal(class(I), 'uint16') | isequal(class(I), 'logical')
    warning(['Input image is of type ', class(I), '. Converting to uint8.']);
    I = map2uint8(I);
end

[M, N] = size(I);

if nargin == 2
    r_max = min([M N]);  %  maximum window size equals the minimum size
    % r_vals = [2.^(0:(log2(r_max) - 1)) r_max];  %  powers of two
    r = round(exp(linspace(log(1), log(r_max), 5)));  %  log-uniform points
end

if any(r<=0) | ~isequal(r, floor(r))
    error('r must contain only positive integers');
end

n_scale = length(r);
L = nan(1, n_scale);

for i_scale=1:n_scale
    r_i = r(i_scale);
    
    switch method
        case 'window'  %  See reference [2] Section 2.2
            % Sliding-window through all the image. For each position we
            % compute the sum of the values in each window as mass.
            m = reshape(conv2(I, ones(r_i, r_i), 'valid'),[], 1);

            L(i_scale) = var(m(:))/mean(m(:))^2 + 1;
        
        case 'local_bin'  %  See reference [3] Section 3.1. Slow because of the
            % nested for loops. Always 1 for r=1.
            n_row = M - r_i + 1;  %  number of rows with starting window
            n_col = N - r_i + 1;  %  number of columns with starting window
            m = nan(n_row, n_col);
            for i=1:n_row  
                for j=1:n_col
                    I_box = I(i:i+r_i-1, j:j+r_i-1);  % window
                    m(i,j) = sum(I_box(:) >= mean(I_box(:)));
                end
            end
            
            L(i_scale) = var(m(:))/mean(m(:))^2 + 1;

        case 'box_3d'  %  See reference [3] Section 3.2. Slow method.         
            L(i_scale) = lacunarity_3d(I, M, N, r_i);
            
        otherwise
            error("Provided method is not supported. Use 'window', 'local_bin' or 'box_3d'.");
    end
    
%     disp([num2str(i_scale) '/' num2str(n_scale)]);
end

if visu
    plot(log(r), log(L), '-o', 'LineWidth', 1.5);
    xlabel('Log(s)');
    ylabel('Log(L)');
    grid on;
end
end

function Lac = lacunarity_3d(I, M, N, r)   
    % Intuitively, we build a tower of cubes of rxrxr and we glide the entire
    % tower through the image. The base of the tower is what we call 'box'. 
    % We compute the occupied space/volume for each cube in the tower and store
    % it in vl. (Cubes can overlap).
    %
    % As we glide the tower, there is an overlap with the previous tower so
    % that we do not need to compute vl entirely. We can just remove the
    % occupied space due to the old rows/columns that do not belong to the new
    % tower and add the new/columns. In this way, the occupied space due to the
    % shared columns/rows is kept and the computation is faster.
    %
    % We then only consider those cubes that have partial volume, i.e. that are
    % neither empty (vl=0) or completely filled (vl=S_max=r^3)
    %
    % Parameters:
    % cc: crossed cubes (intercepted). Fully filled boxes with that intensity
    % cr: part of the height that is used to partially fill cubes (note that
    % this part can be also used to fill the previous cubes below)
    % h: pixel intensity (height)
    % l: maximal number of cubes allowed in the pixel intensity direction (z)
    % L: maximum number of gray values
    % n: number of cubes with mass equal to the index (e.g., n(3) is the number
    % of cubes with mass 3)
    % r: cube size
    % S_max: maximum occupied volume in a cube (s^3) 
    % vl: occupied volume of each cube in the tower (z direction)
    
    L = double(max(I(:)));  
    l = L - r + 1; 
    S_max = r^3; 
    I = double(I);
    
    if l <= 0
        warning(['box size (r=', num2str(r), ') is > than the number of', ...
            ' image levels (L=', num2str(L),'). Use a smaller r.']); 
        Lac = nan;
        return;
    end
    
    vl = zeros(1, l);  
    n = zeros(1, S_max);
    
    % Get the occupied space in each cube of the tower located at the first box
    % (first r columns/rows). We do this by gliding through all pixels in the
    % box
    for x=1:r
        vl = compute_column(I, vl, l, x, 1, r, 1);
    end        
    vl1 = vl;  % store the volume of the first box
            
    % Loop through columns
    first = true;
    for y=r:N
        % Loop through rows
        for x=r:M
            if first
                % If it is the first box/tower --> store the volume directly
                first = false;
            else
                % Remove the volume of the column from previous box that does
                % not belong to this new box.
                vl = compute_column(I, vl, l, x-r, y-r+1, r, -1); 

                %  Add the volume from the columns that are new to this box.
                %  The volume of shared columns is kept.
                vl = compute_column(I, vl, l, x, y-r+1, r, 1);                  
            end
            
            for i=1:l
                % Do not count cubes that are empty (vl=0) or completely filled
                % (vl=S_max)
                if vl(i)~=0 & vl(i)~=S_max
                    n(vl(i)) = n(vl(i)) + 1;  % add boxes with mass vl 
                end
            end
        end
        
        % When we have gone through all the rows check if there are more
        % columns remaining (y < N) and if so compute the volume of the first
        % box at the top.
        if y < N
            vl1 = compute_row(I, vl1, l, 1, y-r+1, r, -1);
            vl1 = compute_row(I, vl1, l, 1, y+1, r, 1);
            vl = vl1;            
        end
        first = true;
    end
    
    % Calculate Lacunarity            
    n_cube = sum(n);  % number of partially filled cubes
    if n_cube == 0
        Lac = 0;  % All cubes were fully filled or empty.
    else
        % Probability density distribution of a partially filled cube having
        % certain occupancy. We just divide the number of cubes with each occupancy
        % by the total number of partially filled cubes.
        c = n./sum(n);
    
        Z1 = sum((1:S_max) .* c);  % mean occupied cubes
        Z2 = sum((1:S_max).^2 .* c);  % power
        
        Lac = Z2/Z1^2;
    end
end

function vl = compute_column(I, vl, l, x, y, r, t)
    % t: 1 or -1 to add or remove volume, respectively
    
    for y1=y:y+r-1
        [cc, cr] = count_boxes(l, I(x,y1), r);      

        % Store/remove the volume of each cube that was filled
        vl(1:cc) = vl(1:cc) + t*r;

        % Store the volume of those cubes partially filled
        i = cc + 1;
        while (i<l) & (cr > 0)
            vl(i) = vl(i) + t*cr;  %  add partial volume
            cr = cr - 1;
            i = i + 1;
        end
    end
end

function vl = compute_row(I, vl, l, x, y, r, t)
    % t: 1 or -1 to add or remove volume, respectively

    for x1 = x:x+r-1
        [cc, cr] = count_boxes(l, I(x1,y), r);

        vl(1:cc) = vl(1:cc) + t*r;
        
        i = cc + 1;
        while (i<l) & (cr > 0)
            vl(i) = vl(i) + t*cr;
            cr = cr - 1;
            i = i + 1;
        end
    end
end

function [cc, cr] = count_boxes(l, h, r)
    % Number of cubes which a given pixel intercepts in function of its gray
    % level intensity    
    cc = h - r + 1;
    
    % If not even a single box is filled
    if cc <= 0
        cc = 0;
        cr = h;
    elseif cc < l  % Not all cubes filled (regular case)
        cr = r - 1;  % Typo in the original paper I believe
    else  %  All cubes filled (no height to partially fill cubes cr=0)
        cr = 0;        
    end
end