function seg = seg_vessels_octa(I, vasculature, verbose)
%SEG_VESSELS_OCTA Segment vasculature from OCTA images
%
%   Usage example seg = seg_vessels_octa(I, depth)
%   Detail explanation goes here
%
%   Input arguments:
%  
%   'I'              OCTA image.
%            
%  
%  
%   Output arguments:
%  
%   'seg'           Description of the argument. Type and purpose.          
%  
%
%   
%   Notes
%   -----
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
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2022

%% Preprocessing
[n, m] = size(I);

I = im2double(I);

I_filt = medfilt2(I, [5, 5]); %% original Maitane [15,15]

% Top Hat
se_size = 75;
se = strel('disk',se_size);
I_th = imtophat(padarray(I_filt,[se_size*2,se_size*2]), se); % fix this
I_th = I_th(se_size*2+1:se_size*2+n, se_size*2+1:se_size*2+m); % fix this
 
switch vasculature
    case 'macro'
        
    % Otsu's threshold
    th  = graythresh(I_th) * 1.7; %% original (1.75)
    seg_1 = I_th > th;

    % Open
    se  = strel('disk', 15); %disc 12 
    seg_2 = imopen(seg_1, se);

    % Reconstruir
    seg_3 = imreconstruct(seg_2, seg_1);

    seg_4 = bwareaopen(seg_3, 5000);


    if verbose
        n_row = 2;
        n_col = 4;
        
        subplot(n_row, n_col, 1); images(I);      title('Original');
        subplot(n_row, n_col, 2); images(I_filt); title('Filtering');
        subplot(n_row, n_col, 3); images(I_th);   title('Top - Hat');
        subplot(n_row, n_col, 4); images(seg_1);  title('Binarized');
        subplot(n_row, n_col, 5); images(seg_2);      title('imopen');
        subplot(n_row, n_col, 5); images(seg_3);      title('imreconstruct');
        subplot(n_row, n_col, 5); images(seg_4);      title('bwareaopen');
    end
       
    seg = seg_4;
    
    case 'micro_macro'
        
        th  = adaptthresh(I_th, 0.4);
        seg = imbinarize(I_th, th);
%         ImBMM=(ImBMM-im2double(FAZ))>0;  % remove FAZ if we have it

    case 'micro'
        se     = strel('disk',7);
        ImBMic = (ImBMM - imdilate(ImB, se)) > 0;
        ImBO   = imopen(ImBMic, se);
        ImBMic = imreconstruct(ImBO,ImBMic);
end