function seg_aura(filename,params)

scanner_type = 'spectralis';

if strcmp(scanner_type, 'cirrus')
    model_filename = 'model_rf_edge_nf27_sr100_nv6_nb8_flat1_norm1_dn0__nt60_mt10_cirrus_28-Oct-2013';
elseif params.minseg & strcmp(scanner_type,'spectralis')
    model_filename = 'model_rf_edge_nf27_sr100_nv2_nb4_flat1_norm1_dn0_exp25_nt20_mt5_04-Dec-2013.mat';
elseif ~params.minseg & strcmp(scanner_type,'spectralis')
    model_filename = 'model_rf_edge_nf27_sr100_nv7_nb8_flat1_norm1_dn0_exp25_nt60_mt10_26-Nov-2013.mat';
else
    error("Unknown scanner");
end

segmethod = params.segmethod;

if ~exist(params.resultfolder,'dir')
    st = mkdir(params.resultfolder);
    if ~st
        error('Error: could not make output results directory %s',params.outputfolder)
    end
end

fids = [];
if params.logfile
    logfilename = sprintf('log%d%d%d%d%d%d.txt',fix(clock));
    logfilename = fullfile(params.resultfolder,logfilename);
    fids = fopen(logfilename,'w');
    if fids == -1
        warning('Could not open log file')
        fids = [];
    end
end

if params.printtoscreen
    fids = cat(1,fids,1);
end

    
[pathstr, name, ext] = fileparts(filename);
    
if strcmp(scanner_type,'spectralis') && ~strcmp(ext,'.vol')
    mfprintf(fids,'Error: invalid file extension (must be .vol) for spectralis files %s\n\n',filename)
elseif strcmp(scanner_type,'cirrus') && ~strcmp(ext,'.img')
    mfprintf(fids,'Error: invalid file extension (must be .img) for cirrus file %s\n\n',filename)
end
    
pdir = fullfile(params.resultfolder,name);
if ~exist(pdir,'dir')
    try
        mkdir(pdir);
    catch err
        mfprintf(fids,'Error: could not create output directory %s!\nError message: %s\n\n',pdir,err.message);
    end
end
savefilename = fullfile(pdir,[name,'_result']);
if params.skip_completed && exist([savefilename '.mat'],'file')
    mfprintf(fids,'Result already exists, skipping file: %s\n\n',filename);
end       

%-- Read data
mfprintf(fids,'Reading file %s...',filename);
try
    if strcmp(scanner_type,'spectralis')
        [header, ~, ~, img_vol] = openVolFast(filename, 'nodisp');
        img_vol(img_vol > 1) = 0;
    else
        [img_vol,vol_info] = octCirrusReader(filename);
        
        header.PID = vol_info.pid;
        header.SizeX = size(img_vol,2);
        header.NumBScans = size(img_vol,3);
        header.SizeZ = size(img_vol,1);
        header.ScaleX = vol_info.vol_res(2);
        header.Distance = vol_info.vol_res(3);
        header.ScaleZ = vol_info.vol_res(1);
        header.ExamTime = vol_info.scan_date;
        header.ScanPosition = vol_info.eye_side;
        header.ScanType = vol_info.scan_type;
    end
catch err
    mfprintf(fids,'error: could not open file!\n Error message: %s\n\n',err.message);
end    
mfprintf(fids,'done!\n');

[h, ~, bscan] = read_vol(filename, 'raw_pixel');

h.SizeX = h.n_ascan;
h.SizeZ = h.n_axial;
h.NumBScans = h.n_bscan;
h.ScaleX = h.scale_x;
h.ScaleZ = h.scale_z;
h.Distance = h.scale_y;
h.ScanPosition = h.eye;
header = h;

%-- Check data size
sizeX = header.ScaleX*(double(header.SizeX)-1);
sizeY = header.Distance*(double(header.NumBScans)-1);
max_size = 7;
if sizeX > max_size || sizeY > max_size
    warning('Scan size too large (%2.2f x %2.2f mm), cropping to 6 x 6 mm',sizeX,sizeY)

    mfprintf(fids,'Cropping volume to 6 x 6 mm...');
    if sizeX > max_size
        npx = round(3/header.ScaleX);
        cpx = round(double(header.SizeX)/2);
        cvals_x = (cpx-npx):(cpx+npx);
        img_vol = img_vol(:,cvals_x,:);
        header.SizeX = size(img_vol,2);
    end

    if sizeY > max_size
        npx = round(3/header.Distance);
        cpx = round(double(header.NumBScans)/2);
        cvals_y = (cpx-npx):(cpx+npx);
        img_vol = img_vol(:,:,cvals_y);
        header.NumBScans = size(img_vol,3);
    end
    mfprintf(fids,'done!\n');
end

vol_size_crop = size(img_vol);
img_vol_crop = img_vol;

if params.resizedata && strcmp(scanner_type,'spectralis')
    if size(img_vol,2) ~= 1024
        mfprintf(fids,'Resizing data to have 1024 A-scans...');
        sc = 1024/size(img_vol,2);
        img_vol = imresize(img_vol,[size(img_vol,1) 1024]);
        img_vol(img_vol<0) = 0;
        header.ScaleX = header.ScaleX/sc;
        header.SizeX = 1024;
        mfprintf(fids,'done!\n');
    end
elseif params.resizedata && strcmp(scanner_type,'cirrus')
    if size(img_vol,2) ~= 512
        mfprintf(fids,'Resizing data to have 512 A-scans...');
        sc = 512/size(img_vol,2);
        img_vol = imresize(img_vol,[size(img_vol,1) 512]);
        img_vol(img_vol<0) = 0;
        header.ScaleX = header.ScaleX/sc;
        header.SizeX = 512;
        mfprintf(fids,'done!\n');
    end
end

    %% Pre-processing steps

    %-- Normalize intensities
    mfprintf(fids,'Normalizing intensities...');
    img_vol = normalizeOCTVolume(img_vol,2,header);
    mfprintf(fids,'done!\n');
    
    if strcmp(scanner_type,'spectralis')
        img_vol = img_vol.^0.25;
    end

    %-- Generate retina mask
    mfprintf(fids,'Detecting retina boundaries...');
    try
        if strcmp(scanner_type,'spectralis')
            [retina_mask, shifts, bds, nbpt] = retinaDetector2_scale(img_vol,header);
        else
            % Slightly different parameters for cirrus
            p.sigma_lat = 2*16.67;
            p.sigma_ax = 0.5*11.6;
            p.distconst = 96.68;
            
            % Need to median filter first
            sz = size(img_vol);
            dn_k = [3 3 1];
            img_vol = permute(img_vol,[2 1 3]);
            img_vol = medfilt2(img_vol(:,:),[dn_k(2) dn_k(1)],'symmetric');
            img_vol = reshape(img_vol,sz(2),sz(1),sz(3));
            img_vol = permute(img_vol,[2 1 3]);
            
            img_vol = im2double(img_vol);
            [retina_mask, shifts, bds, nbpt] = retinaDetector2_scale(img_vol,header,p,false,true);
        end
        mfprintf(fids,'done! (%d outlier points)\n',nbpt);

        if nbpt > 0.5*size(img_vol,2)
            mfprintf(fids,'Warning: poor fit of retina boundaries detected (%d outlier points). Check for artifacts in the data.\n',nbpt);
        end
    catch err
        mfprintf(fids,' error!\nError message: %s\n\n',err.message);
    end

    %-- Flatten to bottom boundary
    mfprintf(fids,'Flattening data...');
    try
        img_vol = retinaFlatten(img_vol,shifts,'linear');
        retina_mask = retinaFlatten(retina_mask,shifts,'nearest');
        mfprintf(fids,'done!\n');
    catch err
        mfprintf(fids,' error!\nError message: %s\n\n',err.message);
    end

    %-- Flip if left eye
    if ~strncmp(header.ScanPosition,'OD',2)
        img_vol = flipdim(img_vol,2);
        retina_mask = flipdim(retina_mask,2);
        bds = flipdim(bds,1);
    end

    %% Random forest classification

    %-- Load trained classifier model
    mfprintf(fids,'Loading classifier model...');
    load(model_filename)
    mfprintf(fids,'done!\n');
    feat_list = train_params.feature_list;

    % -- Compute 3D spatial features
    if ~isempty(feat_list.vol)
        mfprintf(fids,'Calculating 3D features...');        
        feat_vol_3d = calculateFeatures3D(feat_list.vol,img_vol,...
                                          retina_mask,header);   
        mfprintf(fids,'done!\n');
    else
        feat_vol_3d = [];
    end

    %-- Crop based on retina mask boundaries
    sc = sum(sum(retina_mask>0,3),2);
    px_buf = round(60/(header.ScaleZ*1000));
    tv = find(sc>0,1,'first')-px_buf;
    bv = find(sc>0,1,'last')+px_buf;
    if tv < 1, tv = 1; end
    if bv > size(img_vol,1), bv = size(img_vol,1); end
    % Crop
    if ~isempty(feat_vol_3d)
        feat_vol_3d = feat_vol_3d(tv:bv,:,:,:);
    end

    % Run the images through the classifier in small groups for efficiency
    nGroup = 10; % slices per group
    inds = 1:nGroup:size(img_vol,3);

    %-- Run each group of data through the classifier
    votes_vol = zeros([10 bv-tv+1 size(img_vol,2) size(img_vol,3)],'uint8');
    mfprintf(fids,'Running data through boundary classifier...');
    try
        for l = inds
            mfprintf(fids,'%3.0f%% ',l/size(votes_vol,4)*100);

            % Size of current group
            if (l+nGroup-1) > size(img_vol,3)
                ni = size(img_vol,3);
            else 
                ni = (l+nGroup-1);
            end

            % Get image group
            img_grp = img_vol(:,:,l:ni);

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
            [~, votes] = classRF_predict(double(feat_vec),trained_model);

            vv = votes_vol(:,:,:,l:ni);
            vv(:,rm_dil) = permute(votes,[2 1]);
            votes_vol(:,:,:,l:ni) = vv;
        end
        votes_vol = permute(votes_vol,[2 3 4 1]);
        mfprintf(fids,'...done!\n')
    catch err
        mfprintf(fids,'segmentation error!\nError message: %s\n\n',err.message);
    end

    clear trained_model img_vol feat_vol_3d retina_mask vv votes img_grp feat_vec


    %% Graph-based final segmentation
    
    % Find center of fovea for graph construction
    f_cen = foveaFinder(bds(:,:,3)-bds(:,:,1),[header.ScaleX header.Distance]);

    mfprintf(fids,'Calculating final segmentation from votes...')
    try
        if segmethod == 2
            bd_pts = votesToSegmentation(votes_vol(:,:,:,2:end),'optimalsurface_bmoe',header);
        elseif segmethod == 1
            bd_pts = votesToSegmentation(votes_vol(:,:,:,2:end),'optimalsurface_constrained',header,f_cen);
        elseif segmethod == 3
            bd_pts = votesToSegmentation(votes_vol(:,:,:,2:end),'canny',header);
        else
            error(fids,'error')
        end
    catch err
        mfprintf(fids,'segmentation error!\nError message: %s\n\n',err.message);
    end
    mfprintf(fids,'done!\n');

    clear votes_vol votes
    
    %% Results

    %-- Smoothing
    if params.smooth
        mfprintf(fids,'Smoothing results...')
        ks = 25; % 25 um smoothing kernel
        bd_pts = imfilter(bd_pts,fspecial('gaussian',[31 1],ks/header.ScaleX/1000),'replicate');
        bd_pts = imfilter(bd_pts,fspecial('gaussian',[1 15],ks/header.Distance/1000),'replicate');
        mfprintf(fids,'done!\n');
    end

    if ~strncmp(header.ScanPosition,'OD',2)
        bd_pts = flipdim(bd_pts,1);
    end
    bd_pts = bsxfun(@plus,bd_pts+tv-1,shifts);

    stats = thicknessStatistics(bd_pts,header,params.gridradii,params.displaygrid);

    %-- Resize
    % If we resized the original data, resize results back to original size
    if params.resizedata
        if strcmp(scanner_type,'spectralis')
            sc = 1024/vol_size_crop(2);
        elseif strcmp(scanner_type,'cirrus')
            sc = 512/vol_size_crop(2);
        end
        header.ScaleX = header.ScaleX*sc;
        header.SizeX = vol_size_crop(2);
        bd_pts = imresize(bd_pts,[vol_size_crop(2),vol_size_crop(3)]);
    end
    
    bd_pts(bd_pts < 1) = 1;
    bd_pts(bd_pts > size(img_vol_crop,1)) = size(img_vol_crop,1);

    %-- Save
    mfprintf(fids,'Saving results...')

    % Check if file already exists
    if exist([savefilename '.mat'],'file')
        if ~params.overwrite_results
            ii = 0;
            while(exist([savefilename '.mat'],'file'))
                ii = ii + 1;
                savefilename = fullfile(pdir,[name '_result_' num2str(ii)]);
            end
            mfprintf(fids,'\nWarning: results file already exists, saving to file %s.mat...',savefilename)
        else
            mfprintf(fids,'\nWarning: overwriting results file %s.mat...',savefilename)
        end
    end

    % save mat file
    save(savefilename,'bd_pts','header','params')

    % Save results to csv file
    layernames = {'RNFL','GCL+IPL','INL','OPL','ONL','IS','OS','RPE','Total'};
    sectornames = {'center',...
            'inner superior','inner inferior','inner nasal','inner temporal',...
            'outer superior','outer inferior','outer nasal','outer temporal','macula'};
    
    fid = fopen([savefilename '.csv'],'w');
    fprintf(fid,',%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',sectornames{:});
    fprintf(fid,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n',layernames{1},stats(1,:));
    for i = 2:8
        fprintf(fid,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n',layernames{i},stats(i,:));
    end
    fprintf(fid,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f',layernames{9},sum(stats,1));
    fclose(fid);

    % save XML file
    if params.saveXMLfiles
        % Construct label volume
        % Boundary goes at bottom of layer
        bd_pts_r = round(bd_pts);
        label_vol = zeros(size(img_vol_crop),'uint8');
        for ii = 1:size(label_vol,2)
            for jj = 1:size(label_vol,3)
                for kk = 1:(size(bd_pts,3)-1)
                    label_vol((bd_pts_r(ii,jj,kk)+1):(bd_pts_r(ii,jj,kk+1)),ii,jj) = kk;
                end
            end
        end

        xmlparams.type = 'Unsigned Byte';
        xmlparams.res = [header.ScaleX header.ScaleZ header.Distance];
        xmlparams.sliceSpacing = 1;
        xmlparams.sliceThickness = 0;
        label_vol = permute(label_vol,[2 1 3]);
        writeXml(label_vol,xmlparams,savefilename);

        if ~exist([fullfile(pdir,name) '.xml'],'file') || params.overwrite_results
            if strcmp(scanner_type,'spectralis')
                xmlparams.type = 'Float';
                writeXml(permute(img_vol_crop.^0.25,[2 1 3]),xmlparams,fullfile(pdir,name));
            elseif strcmp(scanner_type,'cirrus')
                xmlparams.type = 'Unsigned Byte';
                writeXml(permute(img_vol_crop,[2 1 3]),xmlparams,fullfile(pdir,name));
            end
        end
    end

    mfprintf(fids,'done! Saved to file %s\n',savefilename)

    mfprintf(fids,'\n')
    
    if params.displayresult
        if strcmp(scanner_type,'spectralis')
            octViewer(img_vol_crop.^0.25,bd_pts);
        elseif strcmp(scanner_type,'cirrus')
            octViewer(img_vol_crop,bd_pts);
        end
    end
    
    clear label_vol


if params.logfile == true
    fclose(fids(fids~=1));
end
