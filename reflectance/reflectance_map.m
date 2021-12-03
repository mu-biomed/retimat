function R = reflectance_map(bscan, top, bottom, metric, varargin)
%REFLECTANCE_MAP Createa 2D map (en-face image) of each A-Scan reflectance
%
%   R = reflectance_map(bscan, )
%   Detail explanation goes here
%
%   Input arguments:
%  
%   'bscan'          3D Volume with b-scans images.          
%  
%   'top'            Upper boundary of the layer to be analyzed (in voxel 
%                    units measured from the top of each B-Scan). The 
%                    dimensions must much the provided volume (bscan).          
%  
%   'bottom'         Bottom boundary of the layer to be analyzed (in voxel 
%                    units measured from the top of each B-Scan). The 
%                    dimensions must much the provided volume (bscan).         
%  
%   'metric'         Reflectivity metric to be used. See references for a
%                    detailed description.
%                    Default = 'mean'
%                    Options = ['mean', 'total', 'layer_index']
%
%   'varargin'       Optional parameters from the list:
%
%                    scale_z (double): axial (depth) resolution of the image.
%                    Necessary if the metric is 'total reflectance'.
%
%   Output arguments:
%  
%   'R'              2D matrix with reflectance values.        
%  
%
%   
%   Notes
%   -----
%   This function assumes that your bscan voxel intensity values have already
%   been normalized.
%
%
%   References
%   ----------
%   [1] Vermeer et al., "RPE-Normalized RNFL Attenuation Coefficient Maps 
%   Derived from Volumetric OCT Imaging for Glaucoma Assessment", IOVS, 2012
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
%     R = reflectance_map(bscan, seg.ILM, seg.BM)   
%
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2021


if nargin < 3
    error("This function requires at least 3 input parameters");
elseif nargin == 3
    metric = 'mean';
elseif nargin == 5
    scale_z = varargin{1};
elseif nargin >5
    error("This function takes a maximum of 5 input arguments");
end

[~, n_ascan, n_bscan] = size(bscan);
% Construct mean-reflectance map by looping through all A-Scans
MR = nan(n_bscan, n_ascan);
for i_bscan=1:n_bscan    
    for i_ascan=1:n_ascan
        roi = round(top(i_bscan, i_ascan):bottom(i_bscan, i_ascan));  
        if isempty(roi) | isnan(roi)
            warning('unable to compute A-Scan layer roi');
            MR(i_bscan, i_ascan) = nan;
        else
            MR(i_bscan, i_ascan) = mean(bscan(roi, i_ascan, i_bscan));
        end
    end
end

switch metric
    case 'mean'
        R = MR;
    
    case 'total'
        Thickness = scale_z * abs(top - bottom);        
        TR = MR .* Thickness / scale_z;        
        R = TR;                
        
    case 'layer_index'
        LI = nan(n_bscan, n_ascan);
        
        Thickness = get_thickness_map(top, bottom, header);
        
        for i_bscan=1:n_bscan
            Isa = prctile(reshape(bscan(:,:,i_bscan),1,[]), 99);
            LI(i_bscan,:) = MR(i_bscan,:).*Thickness(i_bscan,:)/Isa;
        end        
        R = LI;
    otherwise
        error("Specified reflectivity metric is not supported");
end

