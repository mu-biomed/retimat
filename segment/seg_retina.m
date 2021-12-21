function mask = seg_retina(I, scale_z, n_win, metric, visu)
%SEG_RETINA Rough segmentation of the region covering the retina
%
%   mask = seg_retina(I, method, visu)
%   Returns a mask covering the retina. 
%
%   Input arguments:
%  
%   'I'              B-Scan image
%
%   'scale_z'        Depth (axial) resolution of the image. Used to limit
%                    leverage known retinal thickness values
%
%   'n_win'          Optional. Number of segmentation windows.
%                    It must be an integer in range [1, n_ascan]
%
%   'metric'         Metric used to segment the retina.
%                    Default: 'mean'
%                    Options: ['mean','otsu']
%
%   'visu'           If true mask is visualized.
%  
%
%   Output arguments:
%  
%   'mask'           Binary mask of the region containing the retina.
%
%   
%   Notes
%   -----
%   Returned mask is intended for rough location and reducing the
%   search space of actual retina layer segmentation process.
%
%
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

if nargin < 5
    visu = false;
end
if nargin < 3
    n_win = 3; 
end
if nargin < 4
    metric = 'mean'; 
end

[N, M] = size(I);

if n_win < 1 | n_win > M
    error(['Number of windows must be in range [1,' num2str(M) ']']); 
end

win_id = discretize(1:M, n_win);
n_pixel = round(0.4/scale_z);  %  number of pixels in window
mask = false(N, M);

% [~, start_col] = unique(win_id);
start_row = 1:N-n_pixel;

if strcmp(metric, 'mean')
    met_fun = @(x, y) 1/mean(x(:));
elseif strcmp(metric, 'otsu')
    met_fun = @(x, y) intraclass_var(x, y);  
else
    error("Unknown metric");    
end

for i_win=1:n_win    
    win_cols = win_id == i_win;
    
    ICV = nan(1, length(start_row));
    for s=1:length(start_row)            
        win_rows = start_row(s) + (0:n_pixel-1);
        s_in = I(win_rows, win_cols);
        s_out = I([1:start_row(s)-1 win_rows(end)+1:N], win_cols);
%         ICV(s) = intraclass_var(s_in, s_out);  
        ICV(s) = met_fun(s_in, s_out);
%         ICV(s) = 1/mean(s_in(:));
    end

    [~, ind] = min(ICV);
    mask(start_row(ind) + (0:n_pixel-1), win_cols) = true;
end


if visu
    imshow(I, 'InitialMag', 'fit'); hold on;            
    green = cat(3, zeros(size(I)), ones(size(I)), zeros(size(I)));
    h = imshow(green);
    set(h, 'AlphaData', mask*0.3)             
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
