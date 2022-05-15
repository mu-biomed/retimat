function R = reflectance_map(bscan, method, metric, seg, varargin)
%REFLECTANCE_MAP Createa 2D map (en-face image) of each A-Scan reflectance
%
%   R = reflectance_map(bscan, method, metric, seg, varargin)
%   Creates a 2D point map of reflectance. Reflectance may be derived from
%   raw voxel intensities, normalized intensities or attenuation
%   coefficient.
%
%   Input arguments:
%  
%   'bscan'          3D Volume with b-scans images.          
%
%   'method'         Method to compute reflectance.
%                    Default = 'raw'
%                    Options = ['raw', 'normalized', 'attenuation']
%
%   'metric'         Metric to be used.
%                    Default = 'mean'
%                    Options = ['mean', 'total', 'layer_index']
%     
%   'seg'            Struct with boundary segmentation data (in voxel units
%                    measured from the top of each B-Scan). The dimensions
%                    must match the provided volume (bscan). If not
%                    provided, all the voxels in each ascan are used.
%
%   'varargin'       Optional parameters from the list:
%
%                    'scale_z' (double): axial (depth) resolution of the image.
%                    Necessary if the metric is 'total reflectance'.
%
%                    'top' Name of the upper boundary of the layer to be
%                    analyzed. It must correspond to a field in seg.        
%  
%                    'bottom' Name of the bottom boundary of the layer to be
%                    analyzed. It must correspond to a field in seg.     
%
%                    If only one optional argument is provided that is assumed
%                    to be scale_z. When two optionals are provided then those
%                    are considered as top/bottom boundaries. When the three
%                    are provided they are assumed in the order: scale_z, top
%                    and bottom.
%
%
%   Output arguments:
%  
%   'R'              2D matrix with reflectance values.        
%  
%
%   
%   Notes
%   -----
%   When using 'normalized', reflectance is corrected based on the vitreous
%   and the RPE layer. See normalize_reflectance.m for details.
%   
%   Attenuation coefficient computation relies on several assumptions on
%   the optical properties of the tissue. See compute_attenuation.m for
%   details.
%
%   'layer_index' metric is described in [1]. It includes a sort of
%   normalization at a bscan level (the 99 percentile is used for 
%   normalization).
%
%
%   References
%   ----------
%   [1] Varga et al., "Investigating Tissue Optical Properties and Texture
%   Descriptors of the Retina in Patients with Multiple Sclerosis", PLos One,
%   2015
%
%
%
%   Example
%   ---------      
%   % Basic reflectance map calculation
%
%     [header, seg, bscan] = read_vol('my_file.vol');
%     R = reflectance_map(bscan, seg, 'MR', 'ILM', 'RPE')   
%
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2022

if nargin < 2
    method = 'raw';
end
if nargin < 3
    metric = 'mean';
end
if nargin < 4
    seg = []; 
end
if nargin == 5
    scale_z = varargin{1};
elseif nargin == 6
    top = varargin{1};
    bottom = varargin{2};    
    top_seg = seg.(top);
    bottom_seg = seg.(bottom);
elseif nargin == 7
    scale_z = varargin{1};
    top = varargin{2};
    bottom = varargin{3};    
    top_seg = seg.(top);
    bottom_seg = seg.(bottom);
elseif nargin > 7
    error("This function takes a maximum of 7 input arguments");
end

[~, n_ascan, n_bscan] = size(bscan);

switch method
    case 'raw'
        disp(''); % Do nothing        
    case 'normalized'        
        bscan = normalize_reflectance(bscan, seg, 'bscan'); 
    case 'attenuation'
        bscan = compute_attenuation(bscan, scale_z, 'vermeer_2014');
    otherwise
        error("Unknown reflectance method");
end

R = nan(n_bscan, n_ascan);
    
% Loop through all a-scans computing the desired metric.
switch metric
    case 'mean'
        if isempty(seg)
            % Permute to return a [n_bscan n_ascan] matrix
            R = permute(squeeze(nanmean(bscan)),[2 1]);
            return;
        end
        
        for b=1:n_bscan
            for a=1:n_ascan
                roi = round(top_seg(b, a):bottom_seg(b, a));
                if isempty(roi) | isnan(roi)
                    warning('unable to compute A-Scan layer roi');
                    R(b, a) = nan;
                else
                    R(b, a) = mean(bscan(roi, a, b));
                end
            end
        end
    
    case 'total'
        if isempty(seg)
            R = permute(squeeze(nansum(bscan)),[2 1]);
            return;
        end
        
        for b=1:n_bscan
            for a=1:n_ascan
                roi = round(top_seg(b, a):bottom_seg(b, a));
                if isempty(roi) | isnan(roi)
                    warning('unable to compute A-Scan layer roi');
                    R(b, a) = nan;
                else
                    R(b, a) = sum(bscan(roi, a, b));
                end
            end
        end
        
% Deprecated: it is basically the same as total reflectance
%     case 'total_Varga'
%         Thickness = scale_z * abs(top - bottom);        
%         TR = MR .* Thickness / scale_z;        
%         Z = TR;                
        
    case 'layer_index' % See [1]   
        if isempty(seg)
            error("layer_index method requires segmentation as input");
        end
        
        Thickness = get_thickness_map(top, bottom, header);
        MR = reflectance_map(bscan, 'raw', 'mean', seg, top, bottom);
        
        for b=1:n_bscan
            Isa = prctile(reshape(bscan(:,:,b),1,[]), 99);
            R(b,:) = MR(b,:).*Thickness(b,:)/Isa;
        end        
    otherwise
        error("Specified reflectivity metric is not supported");
end
