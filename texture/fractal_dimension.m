function FD = fractal_dimension(I, method, visu)
%fractal_dimension Compute Fractal Dimension of a grayscale image
%
%   FD = fractal_dimension(I)
%   Computes the fractal dimension of a 2D image by using the differential box
%   counting method
%
%   Input arguments:
%  
%   'I'              Input image as a 2D matrix. Color images are converted to
%                    gray scale values.          
%
%   'method'         Optionsl. String indicating the method to use for fractal dimension
%                    computation.
%                    Options: 'DBC', 'IR_DBC' (see references [1] and [2])
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
%   dimension of color images. If the goal is that, check other references. The
%   'DBC' method has several assumptions on image size so 'IR_DBC' is prefered.
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
if nargin<=2
    visu = false;
end

% Convert to gray-scale if it is a color image
if length(size(I)) == 3
    warning("Input image has color. It will be converted to gray-scale");
    I = rgb2gray(I);
end

[M, N] = size(I);  %  image dimension

switch method
    case 'DBC'        
        if M ~= N
            error("The number of columns and rows must be equal");
        end
        if mod(M, 2) ~=0
            error(["Number of rows must be multiple of 2 but is " num2str(M)]);
        end

        % List of s values to analyze
        s_max = M/2;
        s_vals = [2.^(1:(log2(s_max) - 1)) s_max];
        
        % Define differential box counting function
        DBC_fun = @(block) max(block.data(:)) - min(block.data(:)) + 1;
        
        % When the blocks do not fit exactly in the image the last block is processed
        % with the number of pixels that fall on it
        
        Nr = nan(1, length(s_vals));
        r = s_vals/M;
        for i_s=1:length(s_vals)
            s = s_vals(i_s);
            
            % From gray level values to box indexes
            I_box = ceil((double(I) + 1)./s);
            
            % Differential box counting
            nr = blockproc(I_box, [s s], DBC_fun);
            Nr(i_s) = sum(nr(:));
            
            % Slow method
            %     % Number of columns of boxes on each side
            %     n_col_side = M/s;
            %     nr2 = nan(n_col_side, n_col_side);
            %     for i=1:n_col_side
            %         for j=1:n_col_side
            %             mask = false(M, M);
            %             ind_i = (1:s) + s*(i-1); % rows (x) of pixels in column (z)
            %             ind_j = (1:s) + s*(j-1); % columns (y) of pixels in column (z)
            %             mask(ind_i, ind_j) = true;
            %
            %             l = max(I_box(mask));
            %             k = min(I_box(mask));
            %             nr2(i, j) = l - k + 1;
            %
            % %             clf;
            % %             subplot(1,2,1);imagesc(I_box);colorbar;
            % %             subplot(1,2,2);imagesc(mask);
            % %             pause;
            %         end
            %     end
            %     Nr2(i_s) = sum(nr2(:));
            
        end

    case 'IR_DBC'
        n_lev = 256;
        
        Q = min([(M*N)^1/3 M N]);  %  maximum value of r (Q)
        r_vals = linspace(2, Q, 20);  %  Define r values
        
        Nr = nan(1, length(r_vals));
        for i_r=1:length(r_vals)
            r = 2;
            m = floor(M/r);
            n = floor(N/r);            
            p = n_lev/r;
            
            % Integer ratio differential box counting
            IR_DBC_fun = @(block) ((max(block.data(:))-min(block.data(:)))/p + 1) * numel(block.data)/(m*n) ;
            nr = blockproc(I, [m n], IR_DBC_fun);
            Nr(i_r) = sum(nr(:));
        end

    otherwise
        error("Unsupported method. Use 'DBC' or 'IR_DBC' instead.");
end

x = log(1./r);
y = log(Nr);
p = polyfit(x, y, 1);
FD = p(1);

if visu
    plot(x, x*p(1) + p(2),'b','linewidth',1.5); hold on;
    scatter(x, y,50,'filled','linewidth',1.5,'MarkerFaceColor','white','MarkerEdgeColor','blue');
    xlabel('Log(1/r)');
    ylabel('Log(Nr)');
    grid on;
    title(['FD = ' num2str(FD)]);
end
end