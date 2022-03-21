function att = compute_attenuation(bscan, scale_z)
%COMPUTE_ATTENUATION Compute voxelwise attenuation coefficient.
%
%   att = compute_attenuation(bscan, header)
%   Compute voxelwise depth-resolved attenuation based on [1].
%
%   Input arguments:
%  
%   'bscan'          Matrix with B-scans of dimensions n_axial x n_ascan x 
%                    n_bscan. If 2D matrix it will be interpreted as a
%                    single B-scan.
%  
%   'scale_z'        Axial resolution in mm.
%  
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
%
%   Example
%   ---------      
%   % Load a .vol file and compute the attenuation coefficient
%
%     [header,~,bscan] = read_vol('my_file.vol');
%     att = compute_attenuation(bscan, header.scale_z);
%
%
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2022

if ~ismatrix(bscan) && ndims(bscan) ~= 3    
    error('Input b-scan must be either a 2D or a 3D matrix');
end

bscan = double(bscan);

% Trick to compute the sum of pixel intensities bellow each pixel
Ic = flip(cumsum(flip(bscan)));
Ic = circshift(Ic,-1,1);
Ic(end,:,:) = nan; % the last row does not have pixels below --> nan

% Compute attenuation
att = 1/(2*scale_z).*log(1 + bscan./Ic);
