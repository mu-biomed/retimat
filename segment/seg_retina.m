function mask = seg_retina(I, method, visu)
%SEG_RETIN Summary of this function goes here
%
%   Usage example OUT = template(IN1)
%   Detail explanation goes here
%
%   Input arguments:
%  
%   'ARG1'           Description of the argument. Type and purpose.          
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
%     I = [1 1 5 6 8 8;2 3 5 7 0 2; 0 2 3 5 6 7];
%     [GLCMS,SI] = graycomatrix(I,'NumLevels',9,'G',[])
%     
%
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2021

if nargin == 2
    visu = false;
end

[N, M] = size(I);


switch method
    case 'k_means'
        mask = zeros(N, M);
        for n=1:M
            X1 = I(:,n);
            X2 = (1:N)';
            X = [X1 X2];
            X = (X - mean(X))./std(X);
            mask(:,n) = kmeans(X, 3);
        end
    case 'otsu'
        mask = zeros(N, M);
        T = graythresh(I);
        mask(I<=T) = 1;
    case 'otsu_rect'         
        % Get scaling
        scale_z = 3.9;
        
        % Define search space
        n_window = 50;
        n_start = 40;
        window = round(linspace(200, 1500, n_window)./scale_z);
        start = round(linspace(1,350, n_start));
        ICV = nan(n_window, n_start); % Inter-class variability
        
        % Brute force search
        for w=1:length(window)
            for s=1:length(start)
                mask = false(size(I));
                endpoint = min([size(I,1) start(s)+window(w)]);
                mask(start(s):endpoint,:) = true;                
                
                ICV(w, s) = intraclass_var(I(mask), I(~mask));                
            end
        end

        % Get minimum
        [~, ind] = min(ICV(:));
        [i, j] = ind2sub(size(ICV),ind);        
        wmin = window(i);
        smin = start(j);
        
        % Build mask
        mask = false(N,M);
        mask(smin:smin+wmin-1,:) = true;
        
        % Plotting
        if visu
            subplot(121); hold on;
            imagesc(I);
            x = [1 1 size(I,2) size(I,2)];
            y = [smin smin+wmin-1 smin+wmin-1 smin];
            p = patch(x,y,1);
            alpha(p, 0.3);
            p.FaceColor = 'green';
            set(gca,'YDir','reverse');
            colormap(gca, 'gray');
            axis off;
            
            subplot(122); hold on;
            imagesc(start, window*scale_z, ICV);
            scatter(smin, wmin*scale_z, 150,'p','filled','MarkerEdgeColor','red','MarkerFaceColor','white','LineWidth',1);
            xlim([start(1) start(end)]);
            ylim(scale_z*[window(1) window(end)]);
            xlabel('start point');
            ylabel('window size');
            set(gca,'FontSize',14);
            colormap(gca, 'turbo');
        end
        
    case 'otsu_ascan'
        scale_z = 3.9;
        
        % Define search space
        n_start = 100;
        window = round(500/scale_z);
        start = round(linspace(1, N-window, n_start));
        ICV = nan(1, n_start); % Inter-class variability
        
        % Brute force search
        mask = false(N, M);
        for i_ascan=1:M
            for s=1:n_start
                                
                s_in = I(start(s) + (0:window-1), i_ascan);
                s_out = I([1:start(s)-1 (start(s)+window):N], i_ascan);
                
%                 ICV(s) = intraclass_var(s_in, s_out);  
                ICV(s) = 1/mean(s_in);
            end
            
            [~, ind] = min(ICV);
            mask(start(ind) + (0:window-1), i_ascan) = true;
        end

        if visu
            imshow(I, 'InitialMag', 'fit'); hold on;            
            green = cat(3, zeros(size(I)), ones(size(I)), zeros(size(I)));
            h = imshow(green);
            set(h, 'AlphaData', mask*0.3)             
        end
        
        
    otherwise
        error("Unsupported method");
end

end

function ICV = intraclass_var(x, y)
                
n_x = numel(x);
n_y = numel(y);
n_all = n_x + n_y;

p_x = n_x/n_all;
p_y = n_y/n_all;

var_x = var(x(:));
var_y = var(y(:));

ICV = p_x*var_x + p_y*var_y;

end
