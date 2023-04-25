function seg = segment_layers_aura(bscan, header, varargin)
% Segment retinal layers based on AURA Tools
% IF YOU USE THIS FUNCTION PLEASE CITE THE ORIGINAL WORK [1][2].
%
%
% Input arguments
% --------------- 
% * **bscan**:        3D matrix with B-scans.        
%  
% * **header**:       Struct with metadata.
%
% * **varargin**:     Optional parameters from the list:
%                       
%   - 'segmethod':    Followed by 1, 2 or 3.
%
%   - 'resizedata':   Resize the data.
%
%   - 'smooth_seg':   Smooth final segmentation results.
%
%   - 'minseg':       Minimum segmentation flag.
%
%
% Output arguments
% ---------------- 
% * **seg**:        Structure with retinal layer segmentation.          
%
%
% Notes
% -----
% This function is built by modyfing the main script of AURA tools. To use
% it download OCTLayerSegmentation from AURA tools
% https://www.nitrc.org/projects/aura_tools and add it to the path. If you
% use this function please cite the original authors of AURA [1][2].
%
%
% References
% ----------
% [1] A. Lang, A. Carass, M. Hauser, E. S. Sotirchos, P. A. Calabresi, 
% H. S. Ying, J. L. Prince, "Retinal layer segmentation of macular OCT
% images using boundary classification", Biomed. Opt. Express, vol. 4,
% No. 7, pp. 1133-1152, 2013. http://dx.doi.org/10.1364/BOE.4.001133
%
% [2] P. Bhargava, A. Lang, O. Al-Louzi, A. Carass, S. Saidha, J. Prince,
% P. A. Calabresi. "Cross-platform comparison of retinal layers in multiple
% sclerosis utilizing a novel open-source optical coherence tomography
% automated segmentation algorithm", 2014 Cooperative Meeting of CMSC and
% ACTRIMS, Dallas, TX, International Journal of MS Care, vol. 16, No. S3,
% pp. 11-12, 2014. http://dx.doi.org/10.7224/1537-2073-16.S3.1
%
%
% Examples
% --------         
% Segment a vol file
% ^^^^^^^^^^^^^^^^^^
% .. code-block:: matlab
%
%     [header, ~, bscan] = read_vol('my_oct.vol')
%     seg = segment_layers_aura(bscan, header);

idx = find(strcmp(varargin, 'segmethod'));
if isempty(idx)
    segmethod = 2;
else
    segmethod = varargin{idx+1};
end

resizedata = any(strcmp('resizedata', varargin));
min_seg    = any(strcmp('min_seg', varargin));
smooth_seg = any(strcmp('smooth_seg', varargin));

% Adapt header for AURA
header.SizeX        = header.n_ascan;
header.SizeZ        = header.n_axial;
header.NumBScans    = header.n_bscan;
header.ScaleX       = header.scale_x;
header.ScaleZ       = header.scale_z;
header.Distance     = header.scale_y;
header.ScanPosition = header.eye;

if strcmp(header.scanner, 'cirrus')
    n_ascan = 512;
    model_filename = 'model_rf_edge_nf27_sr100_nv6_nb8_flat1_norm1_dn0__nt60_mt10_cirrus_28-Oct-2013';
elseif strcmp(header.scanner, 'heidelberg')
    n_ascan = 1024;
    if min_seg
        model_filename = 'model_rf_edge_nf27_sr100_nv2_nb4_flat1_norm1_dn0_exp25_nt20_mt5_04-Dec-2013.mat';
    else
        model_filename = 'model_rf_edge_nf27_sr100_nv7_nb8_flat1_norm1_dn0_exp25_nt60_mt10_26-Nov-2013.mat';
    end
else
    error("Unknown scanner");
end

% Check data size
sizeX = header.ScaleX*(double(header.SizeX)-1);
sizeY = header.Distance*(double(header.NumBScans)-1);
max_size = 7;
if sizeX > max_size || sizeY > max_size
    warning('Scan size too large (%2.2f x %2.2f mm), cropping to 6 x 6 mm',sizeX,sizeY)

    fprintf('Cropping volume to 6 x 6 mm...');
    if sizeX > max_size
        npx = round(3/header.ScaleX);
        cpx = round(double(header.SizeX)/2);
        cvals_x = (cpx-npx):(cpx+npx);
        bscan = bscan(:,cvals_x,:);
        header.SizeX = size(bscan,2);
    end

    if sizeY > max_size
        npx = round(3/header.Distance);
        cpx = round(double(header.NumBScans)/2);
        cvals_y = (cpx-npx):(cpx+npx);
        bscan = bscan(:,:,cvals_y);
        header.NumBScans = size(bscan,3);
    end
    fprintf('done!\n');
end

vol_size_crop = size(bscan);
img_vol_crop = bscan;

% Resize A-scan number
if resizedata && size(bscan,2) ~= n_ascan
    fprintf('Resizing data to have %d A-scans...', n_ascan);
    sc = n_ascan / size(bscan, 2);
    bscan = imresize(bscan,[size(bscan,1) n_ascan]);
    bscan(bscan < 0) = 0;
    header.ScaleX = header.ScaleX/sc;
    header.SizeX = n_ascan;
    fprintf('done!\n');
end

% Preprocessing
bscan(isnan(bscan)) = 0;

bscan = normalizeOCTVolume(bscan, 2, header);

if strcmp(header.scanner, 'heidelberg')
    bscan = bscan.^0.25;
end

% Generate retina mask
if strcmp(header.scanner, 'heidelberg')
    [retina_mask, shifts, bds, nbpt] = retinaDetector2_scale(bscan, header);
else
    % Slightly different parameters for cirrus
    p.sigma_lat = 2*16.67;
    p.sigma_ax = 0.5*11.6;
    p.distconst = 96.68;

    % Need to median filter first
    sz = size(bscan);
    dn_k = [3 3 1];
    bscan = permute(bscan,[2 1 3]);
    bscan = medfilt2(bscan(:,:),[dn_k(2) dn_k(1)],'symmetric');
    bscan = reshape(bscan,sz(2),sz(1),sz(3));
    bscan = permute(bscan,[2 1 3]);

    bscan = im2double(bscan);
    [retina_mask, shifts, bds, nbpt] = retinaDetector2_scale(bscan,header,p,false,true);
end
fprintf('done! (%d outlier points)\n',nbpt);

if nbpt > 0.5*size(bscan,2)
    fprintf('Warning: poor fit of retina boundaries detected (%d outlier points). Check for artifacts in the data.\n',nbpt);
end

% Flattening
bscan = retinaFlatten(bscan, shifts, 'linear');
retina_mask = retinaFlatten(retina_mask, shifts, 'nearest');

% Flip if left eye
if ~strncmp(header.ScanPosition, 'OD', 2)
    bscan = flip(bscan, 2);
    retina_mask = flip(retina_mask, 2);
    bds = flip(bds, 1);
end

% Load trained classifier model
load(model_filename, 'trained_model','train_params')
feat_list = train_params.feature_list;

% Compute 3D spatial features
if ~isempty(feat_list.vol)
    feat_vol_3d = calculateFeatures3D(feat_list.vol,bscan,...
                                      retina_mask,header);   
else
    feat_vol_3d = [];
end

% Crop based on retina mask boundaries
sc = sum(sum(retina_mask > 0, 3), 2);
px_buf = round(60 / (header.ScaleZ * 1000));
tv = find(sc > 0, 1, 'first') - px_buf;
bv = find(sc > 0, 1, 'last') + px_buf;
if tv < 1
    tv = 1;
end
if bv > size(bscan,1)
    bv = size(bscan,1);
end

% Crop
if ~isempty(feat_vol_3d)
    feat_vol_3d = feat_vol_3d(tv:bv,:,:,:);
end

% Run the images through the classifier in small groups for efficiency
nGroup = 10; % slices per group
inds = 1:nGroup:size(bscan, 3);

%-- Run each group of data through the classifier
votes_vol = zeros([10 bv-tv+1 size(bscan,2) size(bscan,3)],'uint8');
fprintf('Running data through boundary classifier...');
for l = inds
    fprintf('%3.0f%% ',l / size(votes_vol, 4) * 100);

    % Size of current group
    if (l+nGroup-1) > size(bscan, 3)
        ni = size(bscan,3);
    else 
        ni = (l+nGroup-1);
    end

    % Get image group
    img_grp = bscan(:,:,l:ni);

    % Calculate 2D features
    feat_vec = calculateFeatures2D(feat_list.img,img_grp,...
                                   retina_mask(:,:,l:ni),header);

    % Crop
    feat_vec = feat_vec(tv:bv,:,:,:);
    if ~isempty(feat_list.vol)
        feat_vec = cat(4,feat_vec,feat_vol_3d(:,:,l:ni,:));
    end

    % Use a dilated retina mask to mask pixels to run the
    % classifier on (expand by 15 pixels)
    rm_dil = imdilate(retina_mask(tv:bv,:,l:ni)>0,...
        strel('line',4*round(15/(header.ScaleZ*1000))+1,90));

    feat_vec = permute(feat_vec,[4 1 2 3]);
    feat_vec = feat_vec(:,rm_dil);
    feat_vec = permute(feat_vec,[2 1]);

    % Predict labels 
    [~, votes] = classRF_predict(double(feat_vec), trained_model);

    vv = votes_vol(:,:,:,l:ni);
    vv(:,rm_dil) = permute(votes,[2 1]);
    votes_vol(:,:,:,l:ni) = vv;
end
votes_vol = permute(votes_vol,[2 3 4 1]);
fprintf('...done!\n')

clear trained_model img_vol feat_vol_3d retina_mask vv votes img_grp feat_vec

% Find center of fovea for graph construction
f_cen = foveaFinder(bds(:,:,3)-bds(:,:,1),[header.ScaleX header.Distance]);

fprintf('Calculating final segmentation from votes...')
if segmethod == 2
    bd_pts = votesToSegmentation(votes_vol(:,:,:,2:end),'optimalsurface_bmoe',header);
elseif segmethod == 1
    bd_pts = votesToSegmentation(votes_vol(:,:,:,2:end),'optimalsurface_constrained',header,f_cen);
elseif segmethod == 3
    bd_pts = votesToSegmentation(votes_vol(:,:,:,2:end),'canny',header);
else
    error("Unknown segmentation method")
end

clear votes_vol votes
    
%-- Smoothing
if smooth_seg
    ks = 25; % 25 um smoothing kernel
    bd_pts = imfilter(bd_pts, fspecial('gaussian', [31 1],ks/header.ScaleX/1000),'replicate');
    bd_pts = imfilter(bd_pts, fspecial('gaussian', [1 15],ks/header.Distance/1000),'replicate');
end

if ~strncmp(header.ScanPosition,'OD',2)
    bd_pts = flip(bd_pts, 1);
end
bd_pts = bsxfun(@plus, bd_pts + tv - 1, shifts);

% Resize back
if resizedata
    sc = n_ascan / vol_size_crop(2);
    header.ScaleX = header.ScaleX * sc;
    header.SizeX = vol_size_crop(2);
    bd_pts = imresize(bd_pts, [vol_size_crop(2), vol_size_crop(3)]);
end

bd_pts(bd_pts < 1) = 1;
bd_pts(bd_pts > size(img_vol_crop,1)) = size(img_vol_crop,1);

bd_pts = permute(bd_pts, [2 1 3]);

seg.ILM     = bd_pts(:, :, 1);
seg.NFL_GCL = bd_pts(:, :, 2);
seg.IPL_INL = bd_pts(:, :, 3);
seg.INL_OPL = bd_pts(:, :, 4);
seg.OPL_ONL = bd_pts(:, :, 5);
seg.ELM     = bd_pts(:, :, 6);
seg.MZ_EZ   = bd_pts(:, :, 7);
seg.OSP_IZ  = bd_pts(:, :, 8);
seg.BM      = bd_pts(:, :, 9);
