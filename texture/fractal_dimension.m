function FD = fractal_dimension(I, method, visu)
%fractal_dimension Compute Fractal Dimension of a grayscale image
%
%   FD = fractal_dimension(I)
%   Computes the fractal dimension of a 2D grayscale image 
%
%
%   Input arguments:
%  
%   'I'              Input image as a 2D matrix. Color images are converted to
%                    gray scale values.          
%
%   'method'         Optional. String indicating the method to use for fractal 
%                    dimension computation.
%                    Options: 'DBC', 'IR_DBC' (see references [1] and [2], 
%                    respectively)
%                    Default = 'IR_DBC'
%
%   'visu'           Optional. If true a plot with the computed values is
%                    shown.
%                    Default = false
%
%
%   Output arguments:
%  
%   'FD'             Fractal dimension. Double.          
%  
% 
%   Notes
%   -----
%   This implementation does not incorporate a procedure to compute the fractal
%   dimension of neither color nor binary images. The 'DBC' method has several
%   assumptions on image size so 'IR_DBC' method is prefered.
%
%
%   References
%   ----------
%   [1] Sarkar N. and Chaudhuri B. "An Efficient Differential Box-Counting 
%   Approach to Compute Fractal Dimension of Image", IEEE Transactions On 
%   Systems, Man and Cybernetics, 1991
%
%   [2] Long M. and Peng F. "A Box-Counting Method with Adaptable Box Height
%   for Measuring the Fractal Feature of Images", Radioengineering, 2013
%   
%
%   Example
%   ---------      
%   % Compute fractal dimension
%
%     I = imread('cameraman.tif');
%     FD = fractal_dimension(I);
%     
%
%  
%   David Romero-Bascones, Biomedical Engineering Department, Mondragon
%   Unibertsitatea, 2021
%   dromero@mondragon.edu

if nargin == 1
    method = 'IR_DBC';
end
if nargin <= 2
    visu = false;
end

if ~isnumeric(I)
    error(['Input image must be numeric but is: ', class(I)]);
end

% Convert to gray-scale if it is a color image
if ndims(I) == 3
    warning("Input image has color. It will be converted to gray-scale");
    I = rgb2gray(I);
end

if isequal(class(I), 'uint8')
    n_lev = 2^8;
elseif isequal(class(I), 'uint16')
    n_lev = 2^16;
else
    warning(['Input image is of type ', class(I), '. Converting to uint8.']);
    I = map2uint8(I);
    n_lev = 2^8;
end

[M, N] = size(I); 

switch method
    case 'DBC'  %  See reference [1]        
        if M ~= N
            max_dim = max([M N]);            
            warning(['Number of columns and rows dont match. Resizing ',... 
                'image to ', num2str(max_dim),' x ', num2str(max_dim),'.']);
            I = imresize(I, max_dim*[1 1]);
            M = max_dim;
%             N = max_dim;  %  not used
        end

        % List of s values to analyze
        s_max = floor(M/2);  %  floor to get an integer that maximizes coverage
%         s_vals = [2.^(1:(log2(s_max) - 1)) s_max];  %  powers of two
        s_vals = round(exp(linspace(log(2), log(s_max), 10)));
        
        % Define differential box counting function
        DBC_fun = @(block) max(block.data(:)) - min(block.data(:)) + 1;
                
        Nr = nan(1, length(s_vals));
        r_vals = s_vals/M;
        for i_s=1:length(s_vals)
            s = s_vals(i_s);
            
            % Compute box height (s_z). (s' in the paper. Equal to s only when
            % the number of gray levels equals image size). If the image is too
            % big the computed height might be smaller than 1!.
            s_z = round(n_lev*s/M);
            if s_z < 1
                warning(['Computed box height is < 1.',...
                    ' Using s_z = 1 instead. Results might be biased.',...
                    ' Consider reducing image size to ', num2str(n_lev),...
                    ' x ', num2str(n_lev),'.']);
                s_z = 1;
            end
            
            % From gray level values to box indexes (starting from 1)
            I_box = ceil((double(I) + 1)./s_z);  %  floor((double(I))./s_z); works as well

            % Differential box counting
            % When the blocks do not fit exactly in the image the last block is 
            % processed with the number of pixels that fall in it        
            nr = blockproc(I_box, [s s], DBC_fun);
            Nr(i_s) = sum(nr(:));                        
        end

    case 'IR_DBC'  %  See reference [2]
        
        Q = min([(M*N)^(1/3) M N]);  %  maximum value of r (Q)       
        r_vals = exp(linspace(log(2), log(Q), 10));  %  r-vals to sample (equally divided in log-space)
        
        Nr = nan(1, length(r_vals));
        for i_r=1:length(r_vals)
            r = r_vals(i_r);
            m = floor(M/r);
            n = floor(N/r);            
            p = n_lev/r;
                       
            % Integer ratio differential box counting
            % It does not matter whether I is in range [0,255] or [0,256] as we
            % are only considering the difference.
            % The numel(block.data)/(m*n) is 1 for full blocks (size = m x n) 
            % and smaller when the block was smaller (at the edges).
            IR_DBC_fun = @(block) ((max(block.data(:))-min(block.data(:)))/p + 1) * numel(block.data)/(m*n) ;
                        
            nr = blockproc(double(I), [m n], IR_DBC_fun);
            
            Nr(i_r) = sum(nr(:));
        end

        % Necessary adjustment (probably a notation error in the paper [2])
        % The "r" used for FD calculation in regression is the inverse of what
        % they use for defining block sizes. Caution!        
        r_vals = 1./r_vals;
        
    otherwise
        error("Unsupported method. Use 'DBC' or 'IR_DBC' instead.");
end

% Compute fractal dimension (FD, D in the literature)
x = log(1./r_vals);
y = log(Nr);
p = polyfit(x, y, 1);
FD = p(1);

if visu
    p = plot(x, x*p(1) + p(2),'linewidth',1.5); hold on;
    scatter(x, y,50,'filled','linewidth',1.5,'MarkerFaceColor','white',...
        'MarkerEdgeColor',p.Color,'HandleVisibility','off');
    xlabel('Log(1/r)');
    ylabel('Log(Nr)');
    grid on;
    title(['FD = ' num2str(FD)]);
    set(gca, 'FontSize', 12);
end
end