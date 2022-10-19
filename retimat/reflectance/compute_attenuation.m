function att = compute_attenuation(bscan, scale_z, method, seg)
%COMPUTE_ATTENUATION Compute voxelwise attenuation coefficient.
%
%   att = compute_attenuation(bscan, header)
%   Compute voxelwise depth-resolved attenuation based on [1] or [3].
%
%   Input arguments:
%  
%   'bscan'          Matrix with B-scans of dimensions n_axial x n_ascan x 
%                    n_bscan. If 2D matrix it will be interpreted as a
%                    single B-scan.
%
%   'scale_z'        Axial resolution in mm.
%
%   'method'         Method used to compute the attenuation coefficient.
%                    Options: 'vermeer_2014', 'vanderschoot_2012'
%                    Default: 'vermeer_2014'
%
%   'seg'            Struct with segmentation data used with the
%                    'vanderschoot_2012' method [3]
%  
%
%   Output arguments:
%  
%   'att'            Attenuation coefficient per voxel.          
%  
%   
%   Notes
%   -----
%   The precise measurement of tissue attenuation coefficient is difficult
%   and still a research topic. This function provides an approximate
%   value based on assumptions on optical tissue properties [1]. You may
%   want to read [2] to check how this method compares to others.
%
%   The method in [3] is currently only supported for RNFL and returns a
%   total attenuation value for each a-scan instead of a 3D matrix.
%
%   For the methods to work it is important to use raw voxel intensity. For
%   instance reading .vol files with the 'raw_pixel' flag.
%
%
%   References
%   ----------
%   [1] Vermeer, Depth-resolved model-based reconstruction of attenuation 
%   coefficients in optical coherence tomography, Biomedical Optics Express
%   2014, https://doi.org/10.1364/BOE.5.000322 
%
%   [2] Chang, Review of methods and applications of attenuation 
%   coefficient measurements with optical coherence tomography, Journal of
%   Biomedical Optics, 2019, https://doi.org/10.1117/1.JBO.24.9.090901
%
%   [3] van der Schoot, The Effect of Glaucoma on the Optical Attenuation
%   Coefficient of the Retinal Nerve Fiber Layer in Spectral Domain Optical
%   Coherence Tomography Images, IOVS, 2012, 
%   https://doi.org/10.1167/iovs.11-8436
%   
%
%
%   Example
%   ---------      
%   % Load a .vol file and compute the attenuation coefficient
%
%     [header,~,bscan] = read_vol('my_file.vol','raw_pixel');
%     att = compute_attenuation(bscan, header.scale_z);
%
%
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2022

if nargin == 2
    method = 'vermeer_2014';
end

if ~ismatrix(bscan) && ndims(bscan) ~= 3
    error('Input b-scan must be either a 2D or a 3D matrix');
end

bscan = double(bscan);

switch method
    case 'vermeer_2014'
        % Trick to compute the sum of pixel intensities bellow each pixel
        Ic = flip(cumsum(flip(bscan),'omitnan'));
        Ic = circshift(Ic,-1,1);
        Ic(end,:,:) = nan; % the last row does not have pixels below --> nan

        % Compute attenuation
        att = 1/(2*scale_z).*log(1 + bscan./Ic);
    case 'vanderschoot_2012'
       % See reference [1] Van der Schoot et al. (IOVS, 2012)
        beta = 2.3; % Empirically computed. Might not work for all.
        
        % RNFL thickness in mm
        d = scale_z * abs(seg.ILM - seg.RNFL_GCL);
        
        % Total reflectance
        T_nfl = reflectance_map(bscan, seg, 'total', scale_z, 'ILM', 'RNFL_GCL');
        T_rpe = reflectance_map(bscan, seg, 'total', scale_z, 'IDZ_RPE', 'BM');
        
        R = T_nfl./T_rpe;
        att = log(R/beta + 1)./(2*d);
        
        % It is already a total metric        
    otherwise
        error("Unsupported attenuation coefficient measurement method");
end
