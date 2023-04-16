function bscan_norm = normalize_reflectance(bscan, seg, method)
% Normalize the intensity (reflectance) of one or more bscans based on the
% reflectance of the RPE layer.
%
%
% Input arguments
% --------------- 
% * **bscan**:       Matrix with B-scans of dimensions n_axial x n_ascan x N_bscan. If 2D matrix then it will be interpreted as a single B-scan.
%  
% * **seg**:         Struct with the segmentation of the ILM, IZ_RPE and BM.
%  
% * **method**:      Method used to normalize.
% 
%   - 'ascan': to normalize each A-scan separately
%   - 'bscan': to normalize each B-scan as a whole (default).
%
%  
% Output arguments
% ---------------- 
% * **bscan_norm**:  Normalized bscan images.          
%  
%   
% Notes
% -----
% This function relies on an accurate segmentation of the ILM, IZ_RPE and
% BM.
%
% The normalization sets the low and hight reflectance references to be the
% vitreous and the RPE, respectively. Can be thought as the % of deviation
% from those parametesr. Values outside [0,100] are expected.
%
% The performance might be somewhat slow due to the necessity of looping 
% through all the A-scans. There might be computational improvements using 
% poly2mask but that needs to be adapted to handle NaNs properly (usual in
% segmentation). It might be also cleaner to use here reflectance_map to
% avoid duplicated code.
%
%
% References
% ----------
% [1] Sharafeldeen, Precise higher‐order refectivity and morphology models
% for early diagnosis of diabetic retinopathy using OCT images, Scientific
% Reports, 2021. https://doi.org/10.1038/s41598-021-83735-7
%
%
% Example
% -------      
% .. code-block:: matlab
% 
%     [~,seg,bscan] = read_vol('my_file.vol');
%     bscan_norm = normalize_reflectance(bscan, seg);

if ismatrix(bscan)
    [n_axial, n_ascan] = size(bscan);
    n_bscan = 1;
elseif ndims(bscan) == 3
    [n_axial, n_ascan, n_bscan] = size(bscan);    
else
    error('Input b-scan must be either a 2D or a 3D matrix');
end

if any(~isfield(seg, {'ILM','IZ_RPE','BM'}))
    error("Segmentation struct must contain 'ILM','IZ_RPE' and 'BM'");
end

bscan_norm = nan(n_axial, n_ascan, n_bscan);

% Round segmentation
seg.ILM = double(round(seg.ILM));
seg.IZ_RPE = double(round(seg.IZ_RPE));
seg.BM = double(round(seg.BM));

for b=1:n_bscan   
    R_vitreous = nan(1, n_ascan);
    R_rpe = nan(1, n_ascan);    
    
    for a=1:n_ascan
        % Get RPE and vitreous segmentation
        ilm = seg.ILM(b, a);
        bm = seg.BM(b, a);
        iz_rpe = seg.IZ_RPE(b, a);

        if any(isnan([ilm bm iz_rpe]))
            warning("NaN values in segmentation. Results might be inaccurate.");
            continue
        end
        
        % Compute vitreous and RPE reflectance
        R_vitreous(a) = mean(bscan(1:ilm, a, b), 'omitnan');
        R_rpe(a) = mean(bscan(iz_rpe:bm, a, b), 'omitnan');        
    end
    
    % Normalize based on [1] (we use 100 instead of 1000 here).
    switch method
        case 'ascan'
            bscan_norm(:,:,b) = 100*(bscan(:,:,b) - R_vitreous)./(R_rpe - R_vitreous);

        case 'bscan'
            R_rpe = mean(R_rpe, 'omitnan');
            R_vitreous = mean(R_vitreous, 'omitnan');
            
            bscan_norm(:,:,b) = 100*(bscan(:,:,b) - R_vitreous)/(R_rpe - R_vitreous);         
        otherwise
            error("Unknown method. Valid options are: 'ascan' or 'bscan'.");
    end
end

bscan_norm = squeeze(bscan_norm);
