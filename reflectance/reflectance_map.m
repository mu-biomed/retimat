function Z = reflectance_map(bscan, seg, metric, varargin)
%REFLECTANCE_MAP Createa 2D map (en-face image) of each A-Scan reflectance
%
%   Z = reflectance_map(bscan, )
%   Detail explanation goes here
%
%   Input arguments:
%  
%   'bscan'          3D Volume with b-scans images.          
%     
%   'seg'            Struct with boundary segmentation data (in voxel units
%                    measured from the top of each B-Scan). The dimensions must
%                    match the provided volume (bscan).
%
%
%   'metric'         Reflectivity metric to be used. See references for a
%                    detailed description.
%                    Default = 'mean'
%                    Options = ['mean', 'total', 'layer_index']
%
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
%   'Z'              2D matrix with reflectance values.        
%  
%
%   
%   Notes
%   -----
%   NFL_RPE_att method relies on several assumptions about scattering
%   properties of the layers [1].
%
%
%   References
%   ----------
%   [1] van der Schoot J. et al., "The Effect of Glaucoma on the Optical
%   Attenuation Coefficient of the Retinal Nerve Fiber Layer in Spectral Domain
%   Optical Coherence Tomography Images", IOVS, 2012
%   doi: https://doi.org/10.1167/iovs.11-8436.
%
%   [2] Varga et al., "Investigating Tissue Optical Properties and Texture
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
%     Z = reflectance_map(bscan, seg, 'MR', 'ILM', 'RPE')   
%
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2021


if nargin < 2
    error("This function requires at least 3 input parameters");
elseif nargin == 2
    metric = 'mean';
elseif nargin == 4
    scale_z = varargin{1};
elseif nargin == 5
    top = varargin{1};
    bottom = varargin{2};    
    top_seg = seg.(top);
    bottom_seg = seg.(bottom);
elseif nargin == 6
    scale_z = varargin{1};
    top = varargin{2};
    bottom = varargin{3};    
    top_seg = seg.(top);
    bottom_seg = seg.(bottom);
elseif nargin >6
    error("This function takes a maximum of 5 input arguments");
end

[~, n_ascan, n_bscan] = size(bscan);

switch metric
    case 'mean'
        % Construct mean-reflectance map by looping through all A-Scans
        MR = nan(n_bscan, n_ascan);

        
        for b=1:n_bscan
            for a=1:n_ascan
                roi = round(top_seg(b, a):bottom_seg(b, a));
                if isempty(roi) | isnan(roi)
                    warning('unable to compute A-Scan layer roi');
                    MR(b, a) = nan;
                else
                    MR(b, a) = mean(bscan(roi, a, b));
                end
            end
        end
        Z = MR;
    
    case 'total'
        % Construct total-reflectance map by looping through all A-Scans
        TR = nan(n_bscan, n_ascan);
        for b=1:n_bscan
            for a=1:n_ascan
                roi = round(top_seg(b, a):bottom_seg(b, a));
                if isempty(roi) | isnan(roi)
                    warning('unable to compute A-Scan layer roi');
                    TR(b, a) = nan;
                else
                    TR(b, a) = sum(bscan(roi, a, b));
                end
            end
        end
        Z = TR;
        
    case 'total_Varga'
        Thickness = scale_z * abs(top - bottom);        
        TR = MR .* Thickness / scale_z;        
        Z = TR;                
        
    case 'layer_index'
        LI = nan(n_bscan, n_ascan);
        
        Thickness = get_thickness_map(top, bottom, header);
        
        for b=1:n_bscan
            Isa = prctile(reshape(bscan(:,:,b),1,[]), 99);
            LI(b,:) = MR(b,:).*Thickness(b,:)/Isa;
        end        
        Z = LI;
    case 'NFL_RPE_att'
        % See reference [1] Van der Schoot et al. (IOVS, 2012)
        beta = 2.3; % Empirically computed. Might not work for all.
        
        % NFL thickness in mm
        d = scale_z * abs(seg.ILM - seg.NFL_GCL);
        
        % Total reflectance
        T_nfl = reflectance_map(bscan, seg, 'total', scale_z, 'ILM', 'NFL_GCL');
        T_rpe = reflectance_map(bscan, seg, 'total', scale_z, 'IDZ_RPE', 'BM');
        
        R = T_nfl./T_rpe;
        att = log(R/beta + 1)./(2*d);
        
        Z = att;
    otherwise
        error("Specified reflectivity metric is not supported");
end

