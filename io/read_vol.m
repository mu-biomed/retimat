function [header, segment, bscan, fundus] = read_vol(file, varargin)
%READ_VOL Read a .vol file from Spectralis OCT (Heidelberg Engineering)
%
%   [header, segment, bscan, fundus] = read_vol(file, options)
%
%   This function reads the header, segmentation and image information 
%   contained in a .vol file. 
%
%   Input arguments:
%  
%   'file'           String containing the path to the .vol file to be read.          
%  
%   'varargin'       Optional parameters from the list:
%                       
%                    'visu': Visualize the scanning patter along with B-Scans
%                    and fundus image (slo).
%
%                    'full_header': Retrieve the original header with all the
%                    parameters (By default only a few important parameters are
%                    retrieved).
%
%                    'get_coordinates': retrieve fundus and A-Scan X, Y coordinates
%
%                    'raw_voxel': return raw pixel reflectance instead of
%                    visualization-adapted values.
%
%   Output arguments:
%  
%   'header'         Structure with .vol file header values.          
%  
%   'segment'        Segmenation data stored in the .vol file.
%
%   'bscan'          3D single image with B-Scans.
%
%   'fundus'            2D fundus image.
%   
%
%   Notes
%   -----
%   Spectralis OCT data can be exported into both E2E and vol format. We
%   recommend using the latter as it provides a better access to the header
%   information.
%
%
%   References
%   ----------
%   [1] 
%
%   Examples
%   ---------      
%   % Read all the information in a .vol file
%
%     [header, segment, bscan, fundus] = read_vol('my_oct.vol')
%     
%   % Read only the header (faster) of the .vol file
%
%     header = read_vol('my_oct.vol')
%
%
%   Originally writen by Radim Kolar, Brno University, Czech Republic
%   Modified by Markus Mayer, Pattern Recognition Lab, University of
%   Erlangen-Nuremberg
%
%   Modified by Kris Sheets, Retinal Cell Biology Lab,
%   Neuroscience Center of Excellence, LSU Health Sciences Center, 
%   New Orleans  
%   
%   Modified by Andrew Lang, Image Analysis and Communications Lab, Johns
%   Hopkins University - Modifications to increase efficiency - 5/17/2012 
%
%   Current version modified by David Romero-Bascones, Biomedical Engineering
%   Department, Mondragon Unibertsitatea, 03/2022
%
%   David Romero-Bascones, Biomedical Engineering Department, Mondragon
%   Unibertsitatea, 2022
%   dromero@mondragon.edu


visu            = any(strcmp('visu', varargin));
full_header     = any(strcmp('full_header', varargin));
get_coordinates = any(strcmp('get_coordinates', varargin));
raw_pixel       = any(strcmp('raw_pixel', varargin));

read_seg    = nargout >= 2;
read_bscan  = nargout >= 3; 
read_fundus = nargout == 4;

fid = fopen(file);
 
% Read header
version        = fread(fid, 12, '*int8');
n_ascan        = fread(fid, 1, '*int32');
n_bscan        = fread(fid, 1, '*int32');
n_axial        = fread(fid, 1, '*int32');
scale_x        = fread(fid, 1, '*double');
scale_y        = fread(fid, 1, '*double');
scale_z        = fread(fid, 1, '*double');
size_x_fundus  = fread(fid, 1, '*int32');
size_y_fundus  = fread(fid, 1, '*int32');
scale_x_fundus = fread(fid, 1, '*double');
scale_y_fundus = fread(fid, 1, '*double');
fov_fundus     = fread(fid, 1, '*int32');
scan_focus     = fread(fid, 1, '*double');
eye            = char(fread(fid, 4, '*uchar')');
exam_time      = fread(fid, 1, '*int64');
scan_pattern   = fread(fid, 1, '*int32');
bscan_hdr_size = fread(fid, 1, '*int32');
id             = char(fread(fid, 16, '*uchar')');
reference_id   = char(fread(fid, 16, '*uchar')');
pid            = fread(fid, 1, '*int32');
patient_id     = char(fread(fid, 21, '*uchar')');
padding        = fread(fid, 3, '*int8');
birth_date     = fread(fid, 1, '*double');
vid            = fread(fid, 1, '*int32');
visit_id       = char(fread(fid, 24, '*uchar')');
visit_date     = fread(fid, 1, '*double');
grid_type      = fread(fid, 1, '*int32');
grid_offset    = fread(fid, 1, '*int32');
spare          = fread(fid, 1832, '*int8');
 
if any([n_bscan n_ascan] > 10000) | any([n_bscan n_ascan] <= 0)
    msg = ['Scan dimensions are unnormal: [', num2str(n_bscan), ',', ...
           num2str(n_ascan), ']. File might be corrupt.'];
    warning(msg);
end

% Build the header
header.n_ascan        = double(n_ascan);
header.n_bscan        = double(n_bscan);
header.n_axial        = double(n_axial);
header.scale_x        = scale_x;
header.scale_y        = scale_y;
header.scale_z        = scale_z;
header.size_x_fundus  = double(size_x_fundus);
header.size_y_fundus  = double(size_y_fundus);
header.scale_x_fundus = scale_x_fundus;
header.scale_y_fundus = scale_y_fundus;
header.fov_fundus     = fov_fundus;
header.patient_id     = deblank(patient_id);
header.eye            = deblank(eye);
header.scan_focus     = scan_focus;

switch scan_pattern
    case 2
        header.fixation      = 'onh';
        header.bscan_pattern = 'peripapillary';
    case 3
        header.fixation      = 'macula';
        header.bscan_pattern = 'raster';
    case 5
        header.fixation      = 'macula';
        header.bscan_pattern = 'star';
    otherwise
        header.fixation      = 'unknown';
        header.bscan_pattern = 'unknown';        
end

if ~any([read_fundus read_seg read_bscan full_header get_coordinates])
    return
end

% Read fundus image 
if read_fundus
    fseek(fid, 2048, -1);
    fundus = fread(fid, size_x_fundus * size_y_fundus, '*uint8');
    fundus = reshape(fundus, size_x_fundus, size_y_fundus);
    fundus = fundus';
end
 
% Read BScans, A-Scan coordinates and boundary segmentation
if read_bscan
    fseek(fid, 2048+(size_x_fundus * size_y_fundus), -1);
    bscan=zeros(n_axial, n_ascan, n_bscan, 'single');
end

start_x = zeros(1, n_bscan, 'double');
start_y = zeros(1, n_bscan, 'double');
end_x   = zeros(1, n_bscan, 'double');
end_y   = zeros(1, n_bscan, 'double');
n_seg   = zeros(1, n_bscan, 'int32');
quality = zeros(1, n_bscan, 'single');
shift   = zeros(1, n_bscan, 'int32');
off_seg = zeros(1, n_ascan, 'single');

if read_seg
    ILM      = zeros(n_bscan, n_ascan, 'single');
    BM       = zeros(n_bscan, n_ascan, 'single');
    RNFL_GCL = zeros(n_bscan, n_ascan, 'single');
    GCL_IPL  = zeros(n_bscan, n_ascan, 'single');
    IPL_INL  = zeros(n_bscan, n_ascan, 'single');
    INL_OPL  = zeros(n_bscan, n_ascan, 'single');
    OPL_ONL  = zeros(n_bscan, n_ascan, 'single');
    ELM      = zeros(n_bscan, n_ascan, 'single');
    MZ_EZ    = zeros(n_bscan, n_ascan, 'single');
    OSP_IZ   = zeros(n_bscan, n_ascan, 'single');
    IZ_RPE   = zeros(n_bscan, n_ascan, 'single');
end

for i_bscan = 1:n_bscan
    ii = i_bscan - 1;
   
    % Read Bscan/Ascan positions in SLO image    
    fseek(fid, 16 + 2048 + (size_x_fundus*size_y_fundus)+(ii*(bscan_hdr_size+n_ascan*n_axial*4)), -1);
    start_x(i_bscan) = fread(fid, 1, '*double'); 
    start_y(i_bscan) = fread(fid, 1, '*double');  
    end_x(i_bscan)   = fread(fid, 1, '*double');
    end_y(i_bscan)   = fread(fid, 1, '*double');  
                     
    % Read Bscan info 
    n_seg(i_bscan)   = fread(fid, 1, '*int32'); % number of segmentated boundaries
    off_seg(i_bscan) = fread(fid, 1, '*int32'); % just in case it is useful
    quality(i_bscan) = fread(fid, 1, '*float32');
    shift(i_bscan)   = fread(fid, 1, '*int32');
        
    % Read Bscan (images)    
    if read_bscan
        fseek(fid, bscan_hdr_size + 2048 + (size_x_fundus*size_y_fundus) + (ii*(bscan_hdr_size+n_ascan*n_axial*4)), -1);
        oct = fread(fid, n_ascan*n_axial, '*float32');
        oct = reshape(oct, n_ascan, n_axial);
        
        % As per Van der Schoot et al. (IOVS, 2012) normalizes the range to
        % [0,1] for visualization. To compute reflectance we need however
        % raw intensities.
        if ~raw_pixel
            oct = oct.^0.25;  
        end
        
        bscan(:,:,i_bscan) = oct';
        bscan(bscan > 1e3) = nan;  % remove outliers at the edges
    end
         
    % Read segmentation     
    fseek(fid, 256 + 2048 + (size_x_fundus*size_y_fundus) + (ii*(bscan_hdr_size+n_ascan*n_axial*4)), -1);
    if read_seg
        seg = (fread(fid, n_seg(i_bscan)*n_ascan, '*float' ))';
        
        ILM(i_bscan,:)      = seg(1:n_ascan);
        BM(i_bscan,:)       = seg((1:n_ascan) + n_ascan);
        RNFL_GCL(i_bscan,:) = seg((1:n_ascan) + 2*n_ascan);
        
        if n_seg(i_bscan) > 3
            GCL_IPL(i_bscan,:) = seg((1:n_ascan) + 3*n_ascan);
            IPL_INL(i_bscan,:) = seg((1:n_ascan) + 4*n_ascan);
            INL_OPL(i_bscan,:) = seg((1:n_ascan) + 5*n_ascan);
            OPL_ONL(i_bscan,:) = seg((1:n_ascan) + 6*n_ascan);
            
            ELM(i_bscan,:)     = seg((1:n_ascan) + 8*n_ascan);
            
            MZ_EZ(i_bscan,:)   = seg((1:n_ascan) + 14*n_ascan);
            OSP_IZ(i_bscan,:)  = seg((1:n_ascan) + 15*n_ascan);
            IZ_RPE(i_bscan,:)  = seg((1:n_ascan) + 16*n_ascan);
        end
    end
end

fclose(fid);

% Return the entire header only if required
if full_header
    header.version        = version;
    header.exam_time      = datestr(double(exam_time(1)/(1e7*60*60*24)+584755+(2/24)));
    header.scan_pattern   = scan_pattern;
    header.bscan_hdr_size = bscan_hdr_size;
    header.id             = deblank(id);
    header.reference_id   = deblank(reference_id);
    header.pid            = pid;
    header.padding        = padding;
    header.birth_date     = datestr(birth_date+693960);
    header.vid            = vid;
    header.visit_id       = deblank(visit_id);
    header.visit_date     = datestr(double(visit_date+693960));
    header.grid_type      = grid_type;
    header.grid_offset    = grid_offset;
    header.spare          = spare;

    header.quality        = quality;
    header.off_seg        = off_seg;    
    header.start_x        = start_x;
    header.start_y        = start_y;
    header.end_x          = end_x;
    header.end_y          = end_y;
    header.n_seg          = n_seg;
    header.shift          = shift;
    header.off_seg        = off_seg;
end

% Compute A-Scan coordinates following the convention:
%  X: left to right
%  Y: inferior to superior
if get_coordinates    
    known_scan_types = {'raster','star','peripapillary'};
    if any(strcmp(header.bscan_pattern, known_scan_types))
        [X_fun, Y_fun, X_oct, Y_oct] = compute_coordinates(header, start_x, start_y,end_x, end_y);
        header.X_fun = X_fun;
        header.Y_fun = Y_fun;
        header.X_oct = X_oct;
        header.Y_oct = Y_oct;
    else
        warning(['Unable to compute coordinates for bscan pattern type ',...
                 header.bscan_pattern]);
    end
end

% Save segmentation data
if read_seg   
    segment.ILM      = ILM;
    segment.BM       = BM;
    segment.RNFL_GCL = RNFL_GCL;
    segment.GCL_IPL  = GCL_IPL;
    segment.IPL_INL  = IPL_INL;
    segment.INL_OPL  = INL_OPL;
    segment.OPL_ONL  = OPL_ONL;
    segment.ELM      = ELM;
    segment.MZ_EZ    = MZ_EZ;
    segment.OSP_IZ   = OSP_IZ;
    segment.IZ_RPE   = IZ_RPE;
    
    % Remove outliers if present
    layers = fields(segment);
    for i=1:11
        segment.(layers{i})(abs(segment.(layers{i})) > n_axial) = nan;
    end

end

% Visualize the acquisition patter if asked
if visu
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 0 scrsz(3) scrsz(4)-70]);
    subplot(1,2,1);
    imshow(fundus);

    for i_bscan = 1:n_bscan       
        start_x_px = round(start_x(i_bscan)/scale_x_fundus); %StartX in pixels
        start_y_px = round(start_y(i_bscan)/scale_y_fundus); %StartY in pixels
        end_x_px   = round(end_x(i_bscan)/scale_x_fundus); %EndX in pixels
        end_y_px   = round(end_y(i_bscan)/scale_y_fundus); %EndY in pixels
        
        subplot(1,2,1); 
        hold on;
        line([start_x_px end_x_px],[start_y_px end_y_px],'color','b');
                
        subplot(1,2,2);
        imshow(bscan(:,:,i_bscan),[0 1]);
        drawnow
    end
end

function [X_fun, Y_fun, X_oct, Y_oct] = compute_coordinates(header, start_x, ...
    start_y, end_x, end_y)
% Compute X,Y coordinates of fundus image and each A-Scan.
% Axis convention
% - X: left to right
% - Y: inferior to superior

size_x_fundus  = double(header.size_x_fundus);
size_y_fundus  = double(header.size_y_fundus);
scale_x_fundus = header.scale_x_fundus;
scale_y_fundus = header.scale_y_fundus;
n_ascan        = header.n_ascan;
n_bscan        = header.n_bscan;

% 1. Define fundus coordinate grid. Initial origin is at upper-left corner of
% the image and B-Scan delimiting values are in mm values
x_range_fundus    = linspace(0, scale_x_fundus * (size_x_fundus - 1), size_x_fundus);
y_range_fundus    = linspace(0, scale_y_fundus * (size_y_fundus - 1), size_y_fundus);   
[X_fun, Y_fun] = meshgrid(x_range_fundus, y_range_fundus);

% 2. Revert the y axis (we want y axis in inferior-superior direction)
Y_fun   = -Y_fun;
start_y = -start_y;
end_y   = -end_y;

% -----------------------------------------------------------------------------
% 3. Set coordinate origin at the center of the fundus image
%------------------------------------------------------------------------------
% Distance to fundus image center
x_offset = scale_x_fundus * (size_x_fundus - 1)/2;
y_offset = -scale_y_fundus * (size_y_fundus - 1)/2;

% Translate all coordinates (fundus + start/end A-scans)
X_fun   = X_fun - x_offset;
Y_fun   = Y_fun - y_offset;
start_x = start_x - x_offset;
start_y = start_y - y_offset;
end_x   = end_x - x_offset;
end_y   = end_y - y_offset;

% -----------------------------------------------------------------------------
% 4. Compute A-Scan coordinates (in fundus image space)
% -----------------------------------------------------------------------------
     
switch header.bscan_pattern
    case 'raster'          
        % Compute A-Scan coordinates
        X_oct = nan(n_bscan, n_ascan);
        Y_oct = nan(n_bscan, n_ascan);
        for i_bscan = 1:n_bscan
            X_oct(i_bscan, :) = linspace(start_x(i_bscan), end_x(i_bscan), n_ascan);
            Y_oct(i_bscan, :) = linspace(start_y(i_bscan), end_y(i_bscan), n_ascan);
        end

        mid_bscan = (n_bscan + 1)/2; % Central B-Scan n_bscan = 25 -> mid_bscan=13        
        mid_ascan = (n_ascan + 1)/2; % Central A-Scan

        % If n_bscan or n_ascan there is no central so use two adjacent scans
        if mod(n_bscan, 2) == 0
            mid_bscan = [floor(mid_bscan) ceil(mid_bscan)];
        end
        
        if mod(n_ascan, 2) == 0
            mid_ascan = [floor(mid_ascan) ceil(mid_ascan)];
        end
        
        % Get Y coordinate of central B-Scan (average for small differences)
        x_oct_center = mean(X_oct(mid_bscan, mid_ascan), 'all');        
        y_oct_center = mean(Y_oct(mid_bscan, mid_ascan), 'all');        
        
    case 'star'             
        % Compute A-Scan coordinates
        X_oct = nan(n_bscan, n_ascan);
        Y_oct = nan(n_bscan, n_ascan);
        for i_bscan = 1:n_bscan
            X_oct(i_bscan, :) = linspace(start_x(i_bscan), end_x(i_bscan), n_ascan);
            Y_oct(i_bscan, :) = linspace(start_y(i_bscan), end_y(i_bscan), n_ascan);
        end
        
        % Get centre Xc based on the first (vertical) B-Scan
        x_oct_center = mean(X_oct(1, :)); % first vertical B-Scan
        
        % Get centre Yc based on the horizontal B-Scan
        hor_bscan = n_bscan/2 + 1;
        y_oct_center = mean(Y_oct(hor_bscan,:));
        
        % Alternative (maybe better)
%         mid_ascan = (n_ascan + 1)/2; % Central A-Scan
%         if mod(n_ascan, 2) == 0
%             mid_ascan = [floor(mid_ascan) ceil(mid_ascan)];
%         end
%         x_oct_center = mean(X_oct(:, mid_ascan),'all'); % Alternative
%         y_oct_center = mean(Y_oct(:, mid_ascan),'all'); % Alternative

    case 'peripapillary'
        % ONH center
        x_onh = end_x(1);
        y_onh = end_y(1);
                
        % Compute A-Scan coordinates        
        radius = sqrt((x_onh - start_x)^2 + (y_onh - start_y)^2);  % onh radius
        
        % Generate a circle (clockwise or anticlockwise depending on eye type)
        if strcmp(header.eye,'OS')
            theta = linspace(0, 2*pi, n_ascan);
        elseif strcmp(header.eye,'OD')
            theta = linspace(pi, -pi, n_ascan);
        end
        [X_oct, Y_oct] = pol2cart(theta, repmat(radius, 1, n_ascan));
        
        % Translate the circle to the onh center
        X_oct = X_oct + x_onh;
        Y_oct = Y_oct + y_onh;
        
        x_oct_center = x_onh;
        y_oct_center = y_onh;
end
 
% -----------------------------------------------------------------------------
% 5. Set coordinate origin at the acquisition center
% -----------------------------------------------------------------------------            

% Translate all coordinates
X_fun = X_fun - x_oct_center;
Y_fun = Y_fun - y_oct_center;

X_oct = X_oct - x_oct_center;
Y_oct = Y_oct - y_oct_center;

% Set to zero horizontal and vertical
X_oct(abs(X_oct) < 1e-6) = 0;
Y_oct(abs(Y_oct) < 1e-6) = 0;
