function [header, segment, bscan, slo] = read_vol(file, varargin)
%read_vol Read .vol file exported from Spectralis OCT (Heidelberg Engineering)
%
%   [header, segment, bscan, slo] = read_vol(file, options)
%
%   This function reads the header, segmentation and image information 
%   contained in the .vol files. 
%
%   Input arguments:
%  
%   'file'           String containing the path to the .vol file to be read.          
%  
%   'varargin'       Optional parameters from the list:
%                       
%                    'visu': Visualize the scanning patter along with B-Scans
%                    and slo image.
%                       
%                    'verbose': Display header info during read.
%
%                    'full_header': Retrieve the original header with all the
%                    parameters (By default only a few important parameters are
%                    retrieved).
%
%                    'coordinates': retrieve fundus and A-Scan X, Y coordinates
%
%   Output arguments:
%  
%   'header'         Structure with .vol file header values.          
%  
%   'segment'        Segmenation data stored in the .vol file.
%
%   'bscan'          3D single image with B-Scans.
%
%   'slo'            2D fundus image.
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
%     file = 'my_oct.vol';
%     [header, segment, bscan, slo] = read_vol(file)
%     
%
%   % Read only the header (faster) of the .vol file
%     file = 'my_oct.vol';
%     header = read_vol(file)
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
%   Department, Mondragon Unibertsitatea, 2021
%
%   David Romero-Bascones, Biomedical Engineering Department, Mondragon
%   Unibertsitatea, 2021
%   dromero@mondragon.edu


% Configure visualization options

visu = any(strcmp('visu', varargin));
verbose = any(strcmp('verbose', varargin));
full_header = any(strcmp('full_header', varargin));
coordinates = any(strcmp('coordinates', varargin));
read_seg = nargout >= 2;
read_bscan = nargout >= 3; 
read_slo = nargout == 4;

% Open the file
fid = fopen(file);
 
% Read header
version = fread(fid, 12, '*int8');
n_ascan = fread(fid, 1, '*int32');
n_bscan = fread(fid, 1, '*int32');
n_axial = fread(fid, 1, '*int32');
scale_x = fread(fid, 1, '*double');
scale_y = fread(fid, 1, '*double');
scale_z = fread(fid, 1, '*double');
size_x_slo = fread(fid, 1, '*int32');
size_y_slo = fread(fid, 1, '*int32');
scale_x_slo = fread(fid, 1, '*double');
scale_y_slo = fread(fid, 1, '*double');
fov_slo = fread(fid, 1, '*int32');
scan_focus = fread(fid, 1, '*double');
eye = char(fread(fid, 4, '*uchar')');
exam_time = fread(fid, 1, '*int64');
scan_pattern = fread(fid, 1, '*int32');
bscan_hdr_size = fread(fid, 1, '*int32');
id = char(fread(fid, 16, '*uchar')');
reference_id = char(fread(fid, 16, '*uchar')');
pid = fread(fid, 1, '*int32');
patient_id = char(fread(fid, 21, '*uchar')');
padding = fread(fid, 3, '*int8');
dob = fread(fid, 1, '*double');
vid = fread(fid, 1, '*int32');
visit_id = char(fread(fid, 24, '*uchar')');
visit_date = fread(fid, 1, '*double');
grid_type = fread(fid, 1, '*int32');
grid_offset = fread(fid, 1, '*int32');
spare = fread(fid, 1832, '*int8');
 
if verbose
    disp('---------------------------------------------');
    disp(['           Version: ' char(version')]);
    disp(['         NumAScans: ' num2str(n_ascan)]);
    disp(['         NumBScans: ' num2str(n_bscan)]);
    disp(['             SizeZ: ' num2str(n_axial)]);
    disp(['            ScaleX: ' num2str(scale_x) ' mm']);
    disp(['            ScaleY: ' num2str(scale_y) ' mm']);
    disp(['            ScaleZ: ' num2str(scale_z) ' mm']);
    disp(['          SizeXSlo: ' num2str(size_x_slo)]);
    disp(['          SizeYSlo: ' num2str(size_y_slo)]);
    disp(['         ScaleXSlo: ' num2str(scale_x_slo) ' mm']);
    disp(['         ScaleYSlo: ' num2str(scale_y_slo) ' mm']);
    disp(['FieldSizeSlo (FOV): ' num2str(fov_slo) 'deg']);
    disp(['         ScanFocus: ' num2str(scan_focus) ' dpt']);
    disp(['               Eye: ' char(eye)]);
    disp(['          ExamTime: ' datestr(double(exam_time(1)/(1e7*60*60*24)+584755+(2/24)))]);
    disp(['       ScanPattern: ' num2str(scan_pattern)]);
    disp(['      BScanHdrSize: ' num2str(bscan_hdr_size) ' bytes']);
    disp(['                ID: ' char(id)]);
    disp(['       ReferenceID: ' char(reference_id)]);
    disp(['               PID: ' num2str(pid)]);
    disp(['         PatientID: ' char(patient_id)]);
    disp(['               DOB: ' datestr(dob+693960)]);
    disp(['               VID: ' num2str(vid)]);
    disp(['           VisitID: ' char(visit_id)]);
    disp(['         VisitDate: ' datestr(double(visit_date+693960))]);
    disp(['          GridType: ' num2str(grid_type)]);
    disp(['        GridOffset: ' num2str(grid_offset)]);
    disp('---------------------------------------------');
end

% Read fundus image (slo)
status = fseek(fid, 2048, -1 );

if read_slo
    slo = fread(fid, size_x_slo*size_y_slo, '*uint8');
    slo = reshape(slo, size_x_slo, size_y_slo);
    slo = slo';
end
 
% Read BScans, A-Scan coordinates and boundary segmentation
status = fseek(fid, 2048+(size_x_slo*size_y_slo), -1);
if read_bscan
    bscan=zeros(n_axial, n_ascan, n_bscan, 'single');
end

start_x = zeros(1, n_bscan, 'double');
start_y = zeros(1, n_bscan, 'double');
end_x = zeros(1, n_bscan, 'double');
end_y = zeros(1, n_bscan, 'double');
n_seg = zeros(1, n_bscan, 'int32');
quality = zeros(1, n_bscan, 'single');
shift = zeros(1, n_bscan, 'int32');
off_seg = zeros(1, n_ascan, 'single');

if read_seg
    ILM = zeros(n_bscan, n_ascan, 'single');
    BM = zeros(n_bscan, n_ascan, 'single');
    NFL_GCL = zeros(n_bscan, n_ascan, 'single');
    GCL_IPL = zeros(n_bscan, n_ascan, 'single');
    IPL_INL = zeros(n_bscan, n_ascan, 'single');
    INL_OPL = zeros(n_bscan, n_ascan, 'single');
    OPL_ONL = zeros(n_bscan, n_ascan, 'single');
    ELM = zeros(n_bscan, n_ascan, 'single');
    MZ_EZ = zeros(n_bscan, n_ascan, 'single');
    PHROS_IDZ = zeros(n_bscan, n_ascan, 'single');
    IDZ_RPE = zeros(n_bscan, n_ascan, 'single');
end

for i_bscan = 1:n_bscan
    ii = i_bscan - 1;
   
    % Read Bscan/Ascan positions in SLO image    
    status = fseek(fid, 16 + 2048 + (size_x_slo*size_y_slo)+(ii*(bscan_hdr_size+n_ascan*n_axial*4)), -1);
    start_x(i_bscan) = fread(fid, 1, '*double'); 
    start_y(i_bscan) = fread(fid, 1, '*double');  
    end_x(i_bscan) = fread(fid, 1, '*double');
    end_y(i_bscan) = fread(fid, 1, '*double');  
                     
    % Read Bscan info 
    n_seg(i_bscan) = fread(fid, 1, '*int32'); % number of segmentated boundaries
    off_seg(i_bscan) = fread(fid, 1, '*int32'); % just in case it is useful
    quality(i_bscan) = fread(fid, 1, '*float32');
    shift(i_bscan) = fread(fid, 1, '*int32');
        
    % Read Bscan (images)    
    if read_bscan
        status = fseek( fid, bscan_hdr_size + 2048 + (size_x_slo*size_y_slo) + (ii*(bscan_hdr_size+n_ascan*n_axial*4)), -1);
        oct = fread(fid, n_ascan*n_axial, '*float32');
        oct = reshape(oct, n_ascan, n_axial);
        oct = oct.^0.25;  % necessary to enhance the contrast but not sure why is that value concretely.
        bscan(:,:,i_bscan)=oct';
    end
         

    % Read segmentation     
    status = fseek(fid, 256 + 2048 + (size_x_slo*size_y_slo) + (ii*(bscan_hdr_size+n_ascan*n_axial*4)), -1 );
    if read_seg
        seg = (fread(fid, n_seg(i_bscan)*n_ascan, '*float' ))';
        
        ILM(i_bscan,:) = seg(1:n_ascan);
        BM(i_bscan,:) = seg((1:n_ascan) + n_ascan);
        NFL_GCL(i_bscan,:) = seg((1:n_ascan) + 2*n_ascan);
        
        if n_seg(i_bscan) > 3
            GCL_IPL(i_bscan,:) = seg((1:n_ascan) + 3*n_ascan);
            IPL_INL(i_bscan,:) = seg((1:n_ascan) + 4*n_ascan);
            INL_OPL(i_bscan,:) = seg((1:n_ascan) + 5*n_ascan);
            OPL_ONL(i_bscan,:) = seg((1:n_ascan) + 6*n_ascan);
            
            ELM(i_bscan,:) = seg((1:n_ascan) + 8*n_ascan);
            
            MZ_EZ(i_bscan,:) = seg((1:n_ascan) + 14*n_ascan);
            PHROS_IDZ(i_bscan,:) = seg((1:n_ascan) + 15*n_ascan);
            IDZ_RPE(i_bscan,:) = seg((1:n_ascan) + 16*n_ascan);
        end
    end
end

fclose(fid);

% Build the header

header.n_ascan = double(n_ascan);
header.n_bscan = double(n_bscan);
header.n_axial = double(n_axial);
header.scale_x = scale_x;
header.scale_y = scale_y;
header.scale_z = scale_z;

header.size_x_slo = double(size_x_slo);
header.size_y_slo = double(size_y_slo);
header.scale_x_slo = scale_x_slo;
header.scale_y_slo = scale_y_slo;
header.fov_slo = fov_slo;

header.patient_id = patient_id;
header.eye = eye(1:2);
header.scan_focus = scan_focus;

switch scan_pattern
    case 2
        header.scan_type = 'peripapillary';
    case 3
        header.scan_type = 'macula_raster';
    case 5
        header.scan_type = 'macula_star';
    otherwise
        header.scan_type = 'unknown';
end

% Return the entire header only if required
if full_header
    header.version = version;
    header.exam_time = exam_time;
    header.scan_pattern = scan_pattern;
    header.bscan_hdr_size = bscan_hdr_size;
    header.id = id;
    header.reference_id = reference_id;
    header.pid = pid;
    header.padding = padding;
    header.dob = dob;
    header.vid = vid;
    header.visit_id = visit_id;
    header.visit_date = visit_date;
    header.grid_type = grid_type;
    header.grid_offset = grid_offset;
    header.spare = spare;

    header.quality = quality;
    header.off_seg = off_seg;    
    header.start_x = start_x;
    header.start_y = start_y;
    header.end_x = end_x;
    header.end_y = end_y;
    header.n_seg = n_seg;
    header.shit = shift;
    header.off_seg = off_seg;
end

% Compute A-Scan coordinates following the convention:
%  X: left to right
%  Y: inferior to superior
if coordinates    
    known_scan_types = {'macula_raster','macula_star','peripapillary'};
    if any(strcmp(header.scan_type, known_scan_types))
        [X_slo, Y_slo, X_oct, Y_oct] = get_coordinates(header, start_x, start_y, ...
        end_x, end_y);
        header.X_slo = X_slo;
        header.Y_slo = Y_slo;
        header.X_oct = X_oct;
        header.Y_oct = Y_oct;
    else
        warning(['Unable to compute coordinates for scan_type: ' header.scan_type]);
    end
end

% Save segmentation data
if read_seg   
    segment.ILM = ILM;
    segment.BM = BM;
    segment.NFL_GCL = NFL_GCL;
    segment.GCL_IPL = GCL_IPL;
    segment.IPL_INL = IPL_INL;
    segment.INL_OPL = INL_OPL;
    segment.OPL_ONL = OPL_ONL;
    segment.ELM = ELM;
    segment.MZ_EZ = MZ_EZ;
    segment.PHROS_IDZ = PHROS_IDZ;
    segment.IDZ_RPE = IDZ_RPE;
    
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
    imshow(slo);

    for i_bscan = 1:n_bscan       
        start_x_px = round(start_x(i_bscan)/scale_x_slo); %StartX in pixels
        start_y_px = round(start_y(i_bscan)/scale_y_slo); %StartY in pixels
        end_x_px = round(end_x(i_bscan)/scale_x_slo); %EndX in pixels
        end_y_px = round(end_y(i_bscan)/scale_y_slo); %EndY in pixels
        
        subplot(1,2,1); 
        hold on;
        line([start_x_px end_x_px],[start_y_px end_y_px],'color','b');
                
        subplot(1,2,2);
        imshow(bscan(:,:,i_bscan),[0 1]);
        drawnow
    end
end
end

function [X_slo, Y_slo, X_oct, Y_oct] = get_coordinates(header, start_x, ...
    start_y, end_x, end_y)
% Compute X,Y coordinates of fundus image and each A-Scan.
% Axis convention
% - X: left to right
% - Y: inferior to superior

size_x_slo = double(header.size_x_slo);
size_y_slo = double(header.size_y_slo);
scale_x_slo = header.scale_x_slo;
scale_y_slo = header.scale_y_slo;
n_ascan = header.n_ascan;
n_bscan = header.n_bscan;

% 1. Define fundus coordinate grid. Initial origin is at upper-left corner of
% the image and B-Scan delimiting values are in mm values
x_range_slo = linspace(0, scale_x_slo * (size_x_slo - 1), size_x_slo);
y_range_slo = linspace(0, scale_y_slo * (size_y_slo - 1), size_y_slo);   
[X_slo, Y_slo] = meshgrid(x_range_slo, y_range_slo);

% 2. Revert the y axis (we want y axis in inferior-superior direction)
Y_slo = -Y_slo;
start_y = -start_y;
end_y = -end_y;

% -----------------------------------------------------------------------------
% 3. Set coordinate origin at the center of the fundus image
%------------------------------------------------------------------------------
% Distance to fundus image center
x_offset = scale_x_slo * (size_x_slo - 1)/2;
y_offset = -scale_y_slo * (size_y_slo - 1)/2;

% Translate all coordinates (fundus + start/end A-scans)
X_slo = X_slo - x_offset;
Y_slo = Y_slo - y_offset;
start_x = start_x - x_offset;
start_y = start_y - y_offset;
end_x = end_x - x_offset;
end_y = end_y - y_offset;

% -----------------------------------------------------------------------------
% 4. Compute A-Scan coordinates (in fundus image space)
% -----------------------------------------------------------------------------
     
switch header.scan_type
    case 'macula_raster'          
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
        x_oct_center = mean(X_oct(mid_bscan, mid_ascan),'all');        
        y_oct_center = mean(Y_oct(mid_bscan, mid_ascan),'all');        
        
    case 'macula_star'             
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
X_slo = X_slo - x_oct_center;
Y_slo = Y_slo - y_oct_center;

X_oct = X_oct - x_oct_center;
Y_oct = Y_oct - y_oct_center;

% Set to zero horizontal and vertical
X_oct(abs(X_oct) < 1e-6) = 0;
Y_oct(abs(Y_oct) < 1e-6) = 0;
       
end
