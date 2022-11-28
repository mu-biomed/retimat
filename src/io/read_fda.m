function [header, segment, bscan, fundus] = read_fda(file, varargin)
%READ_FDA Read Topcon OCT files (fda)
%
%   [header, segment, bscan, fundus] = read_fda(file)
%
%   Reads the header, segmentation and images (bscan + fundus) contained in  
%   an .fda Topcon file. 
%
%   Input arguments:
%  
%   'file'           Path of the .fda file to read.
%
%   'varargin'       Optional flags from the list:
%
%                    'verbose': If provided, reading info is displayed.
%  
%                    'get_coordinates': If provided A-scan coordinates are 
%                    returned. This only works if the scanning pattern is
%                    horizontal raster (pending to implement others).
%
%
%   Output arguments:
%  
%   'header'         Structure with .fda file header values.          
%  
%   'segment'        Segmenation data stored in the .fda file.
%
%   'bscan'          3D single image with B-Scans.
%
%   'fundus'         2D fundus image.
%
%   
%   Notes
%   -----
%   This code is heavily based on [1,2], which were developed by reverse 
%   engineering fda files. Therefore, data read using this function may be
%   incomplete/incorrect.
%
%   The reading process is as follows:
%   1. Identify data chunks present in file (not all are always present)
%   2. Read specific chunks with the relevant information
%
%
%   References
%   ----------
%   [1] uocte documentation
%   https://bitbucket.org/uocte/uocte/wiki/Topcon%20File%20Format
%
%   [2] OCT-Converter
%   https://github.com/marksgraham/OCT-Converter
%
%
%   Example
%   ---------      
%   % Read fda file
%
%     [header, segment, bscan, fundus] = read_fda(file)
%     
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2022

verbose         = any(strcmp('verbose', varargin));
get_coordinates = any(strcmp('get_coordinates', varargin));

if ~isfile(file)
    error('Unable to find the file. Check the path.');
end

fid = fopen(file);

%% Preliminary values (not useful)
type     = fread(fid, 4, '*char')';  % FOCT
fixation = fread(fid, 3, '*char')';  % FDA (not really fixation)
unknown  = fread(fid, 2, '*uint32');  

%% Discover all chunks (name-size-position)
chunks = struct;
i = 1;

more_chunks = true;
while more_chunks
    chunk_name_size = fread(fid,1,'*uint8');
    
    if isempty(chunk_name_size) | chunk_name_size == 0
        break;
    end
    
    chunks(i).name = char(fread(fid, chunk_name_size, '*uint8')');    
    chunks(i).size = fread(fid, 1, '*uint32');  % Bytes!
    chunks(i).pos = ftell(fid);       
    fseek(fid, chunks(i).pos + chunks(i).size, 'bof');
    
    i = i + 1;
end
chunks = struct2table(chunks);

if verbose
    disp([num2str(size(chunks,1)) ' chunks found in file']);
    disp(chunks.name);
end

%% Data reading (using corresponding chunks)
header = read_header(fid, chunks);
if get_coordinates
    [X, Y] = get_ascan_coordinates(header);
    header.X_oct = X;
    header.Y_oct = Y;
end
header = reorder_header(header);

if nargout == 1
    return
end

segment = read_segmentation(fid, chunks, header);
if nargout == 2
    return
end

bscan = read_bscan(fid, chunks);
if nargout == 3
    return
end

fundus = read_fundus(fid, chunks);
fclose(fid);

function header = read_header(fid, chunks)
% We construct a header with data from different chunks
% Important parameters are eye, fixation , pattern and dimensions

% 1. Eye + scan date
idx = find(strcmp(chunks.name, '@CAPTURE_INFO_02'));
if isempty(idx)
    warning("@CAPTURE_INFO_02 chunk was not found, header will be incomplete");    
else
    fseek(fid, chunks.pos(idx), 'bof');
    data = read_chunk(fid, '@CAPTURE_INFO_02');

    header.eye       = data.eye;
    header.scan_date = datetime(data.year, data.month, data.day, ...
                                data.hour,data.minute,data.second);
end

% 2. Version (not very useful)                        
idx = find(strcmp(chunks.name, '@FDA_FILE_INFO'));
if ~isempty(idx)
    fseek(fid, chunks.pos(idx), 'bof');
    data = read_chunk(fid, '@FDA_FILE_INFO');
    header.version = deblank(data.version);
end

% 3. Some dimensions
idx = find(strcmp(chunks.name, '@PARAM_SCAN_04'));
if isempty(idx)
    warning("@PARAM_SCAN_04 chunk was not found. header will be incomplete");    
else
    fseek(fid, chunks.pos(idx), 'bof');
    data = read_chunk(fid, '@PARAM_SCAN_04');

    header.size_x   = data.size_x;
    header.size_y   = data.size_y;
    header.scale_z  = data.scale_z * 1e-3;  % in mm
    header.n_bscan  = double(data.n_bscan);
    header.fixation = data.fixation;
    
    header.scale_y  = header.size_y / (header.n_bscan - 1);
end

% 4. More dimensions (from B-Scan image chunk itself)
idx = find(strcmp(chunks.name, '@IMG_JPEG'));
if isempty(idx)
    warning("@IMG_JPEG chunk was not found. header will be incomplete");    
else
    fseek(fid, chunks.pos(idx), 'bof');

    scan_pattern  = fread(fid, 1, '*uint8');
    switch scan_pattern
        case 0
            header.bscan_pattern = 'line';                            
        case 2
            header.bscan_pattern = 'raster'; % horizontal                
        case 3
            header.bscan_pattern = 'star';
        case 6
            header.bscan_pattern = 'raster_v'; % vertical
        case 7
            header.bscan_pattern = '7_line';
        case 11
            header.bscan_pattern = 'two_5_line';
        otherwise
            header.bscan_pattern = 'unknown';
    end
    
    unknown         = fread(fid, 2, '*uint32');        
    header.n_ascan  = double(fread(fid, 1, '*uint32'));  % width
    header.n_axial  = double(fread(fid, 1, '*uint32'));  % height
    header.n_bscan  = double(fread(fid, 1, '*uint32'));  % num_slices
    
    header.scale_x  = header.size_x / (header.n_ascan -1);
    header.size_z   = header.scale_z * header.n_axial; 
end     
      
function header_ord = reorder_header(header)
% Reorder header fields alphabetically to read them easier
vars = fields(header);
vars = sort(vars);

for i=1:length(vars)
    header_ord.(vars{i}) = header.(vars{i});
end

function segment = read_segmentation(fid, chunks, header)
% Stored in multiple @CONTOUR_INFO chunks (one per segmented boundary)
% 
% Boundary names are defined to match the APOSTEL 2.0 convention:
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8279566/figure/F2/
%
% Segmentation is measured in pixels from the bottom of the image
% (as opposed to Spectralis). For plotting use n_axial - segment.layer

boundary_name = struct('RETINA_1',       'ILM',...
                       'RETINA_2',       'IZ_RPE',...
                       'RETINA_3',       'MZ_EZ',...
                       'RETINA_4',       'BM',...
                       'MULTILAYERS_1',  'ILM',...
                       'MULTILAYERS_2',  'RNFL_GCL',...
                       'MULTILAYERS_3',  'GCL_IPL',...
                       'MULTILAYERS_4',  'IPL_INL',...
                       'MULTILAYERS_5',  'MZ_EZ',...
                       'MULTILAYERS_6',  'IZ_RPE',...
                       'MULTILAYERS_7',  'BM',...
                       'MULTILAYERS_8',  'INL_OPL',...
                       'MULTILAYERS_9',  'ELM',...
                       'MULTILAYERS_10', 'CSI');
               
idx = find(strcmp(chunks.name, '@CONTOUR_INFO'));
n_seg_chunk = length(idx);
if n_seg_chunk == 0
    segment = [];
    return;
end

segment = struct;
for i=1:n_seg_chunk
    fseek(fid, chunks.pos(idx(i)), 'bof');
    data = read_chunk(fid, '@CONTOUR_INFO');   
    
    z = double(data.segment);   % uint16 -> double
    z(z > 10000) = nan;     % invalid values -> NaN
    z = permute(z, [2 1]);  % arange it into [bscan x ascan]
    z = flip(z, 1);         % flip top-bottom
    
    z = header.n_axial - z;
    segment.(boundary_name.(data.id)) = z;
end
        
function bscan = read_bscan(fid, chunks)

idx = find(strcmp(chunks.name, '@IMG_JPEG'));
fseek(fid, chunks.pos(idx), 'bof');
data = read_chunk(fid, '@IMG_JPEG');
bscan = data.bscan;

function fundus = read_fundus(fid, chunks)

idx = find(strcmp(chunks.name, '@IMG_FUNDUS'));
fseek(fid, chunks.pos(idx), 'bof');
data = read_chunk(fid, '@IMG_FUNDUS');

fundus = data.fundus;

function data = read_chunk(fid, chunk_name)

% Known chunks:
% '@CAPTURE_INFO_02': eye, date
% '@CONTOUR_INFO': segmentation data
% '@FDA_FILE_INFO'
% '@HW_INFO_03'
% '@IMG_FUNDUS'
% '@IMG_JPEG'
% '@IMG_TRC_02'
% '@PARAM_SCAN_04'
% '@PATIENT_INFO_02'
% '@PATIENT_INFO_02'
%
%
% Unknown chunks:
% '@ANGIO_TRACKING_INFO'
% '@ANTERIOR_CALIB_INFO'
% '@ALIGN_INFO'
% '@CONTOUR_MASK_INFO'
% '@EFFECTIVE_SCAN_RANGE'
% '@FAST_Q2_INFO'
% '@FDA_DISC_SEGMENTATION'
% '@FOCUS_MOTOR_INFO'
% '@GLA_LITTMANN_01'
% '@IMG_EN_FACE'
% '@IMG_MOT_COMP_03'
% '@IMG_PROJECTION'
% '@MAIN_MODULE_INFO'
% '@PARAM_EN_FACE'
% '@PARAM_OBS_02'
% '@PARAM_TRC'
% '@PARAM_TRC_02'        
% '@PATIENTEXT_INFO'
% '@REF_IMG_SCAN'
% '@REGIST_INFO'
% '@REPORT_INFO'
% '@RESULT_CORNEA_CURVE'
% '@SCAN_POS_COMP_DATA'
% '@THUMBNAIL'
% '@TOPQEXT_INFO'
           

data = struct();
switch chunk_name

    case '@CAPTURE_INFO_02'
        eye = fread(fid, 1, '*uint8');
        if eye == 0
            data.eye = 'OD';
        elseif eye == 1
            data.eye = 'OS';
        else
            data.eye = 'unknown';
            warning('Unknown eye value');
        end
        
        data.unknown = fread(fid, 1, '*uint8');  % always 2 apparently        
        data.zero    = fread(fid, 52, '*uint16');        
        data.year    = fread(fid, 1, '*uint16');
        data.month   = fread(fid, 1, '*uint16');
        data.day     = fread(fid, 1, '*uint16');
        data.hour    = fread(fid, 1, '*uint16');
        data.minute  = fread(fid, 1, '*uint16');
        data.second  = fread(fid, 1, '*uint16');
        
    case '@CONTOUR_INFO'
        % Segmentation
        data.id      = fread(fid, 20, '*char')';        
        data.type    = fread(fid, 1, '*uint16');
        data.n_ascan = fread(fid, 1, '*uint32');  % width
        data.n_bscan = fread(fid, 1, '*uint32');  % height
        data.size    = fread(fid, 1, '*uint32');  % useful?
        
        data.id(data.id == 0) = ' ';  % necessary to trim it later
        data.id = strtrim(data.id);
        
        n_voxel = data.n_ascan * data.n_bscan;
        if data.type == 0
            segment = fread(fid, n_voxel, '*uint16');
        else
            segment = fread(fid, n_voxel, '*float64');
        end
                
        data.segment = reshape(segment, data.n_ascan, data.n_bscan);      
        data.seg_version = fread(fid, 32, '*char')';
   
    case '@FDA_FILE_INFO'
        data.info1   = fread(fid, 1, '*uint32');
        data.info2   = fread(fid, 1, '*uint32');
        data.version = fread(fid, 32, '*char')';        
        
    case '@HW_INFO_03'
        data.model_name      = char(fread(fid, 16, '*uint8')');        
        data.serial_number   = fread(fid, 16, '*uint8');        
        data.zeros           = fread(fid, 32, '*uint8');        
        data.version         = fread(fid, 16, '*uint8');        
        data.build_year      = fread(fid, 1, '*uint16');        
        data.build_month     = fread(fid, 1, '*uint16');        
        data.build_day       = fread(fid, 1, '*uint16');        
        data.build_hour      = fread(fid, 1, '*uint16');        
        data.build_minute    = fread(fid, 1, '*uint16');        
        data.build_second    = fread(fid, 1, '*uint16');        
        data.zeros2          = fread(fid, 8, '*uint8');        
        data.version_numbers = fread(fid, 5*16, '*uint8');        
        
    case '@IMG_FUNDUS'
        data.width          = fread(fid, 1, '*uint32');
        data.height         = fread(fid, 1, '*uint32');
        data.bits_per_pixel = fread(fid, 1, '*uint32');
        data.num_slices     = fread(fid, 1, '*uint32');
        data.val            = fread(fid, 1, '*uint32');        
        data.size           = fread(fid, 1, '*uint32');        
        data.fundus_bytes   = fread(fid, data.size, '*uint8');

        % Ugly trick to decode JPEG
        temp_file = 'temp.jpg';
        fid2 = fopen(temp_file, 'wb');
        fwrite(fid2, data.fundus_bytes, 'uint8');             
        fclose(fid2);        
        fundus = imread(temp_file);          
        delete(temp_file);
        
        % BGR to RGB 
        data.fundus = flip(fundus, 3);        

    case '@IMG_JPEG'
        % Read preliminary data
        scan_pattern  = fread(fid, 1, '*uint8');
        
        switch scan_pattern
            case 2
                data.scan_pattern = 'star';                
            case 3
                data.scan_pattern = 'raster';
            otherwise
                data.scan_pattern = 'unknown';
        end
        
        data.unknown1 = fread(fid, 1, '*uint32');
        data.unknown2 = fread(fid, 1, '*uint32');
        data.n_ascan  = fread(fid, 1, '*uint32');  % width
        data.n_axial  = fread(fid, 1, '*uint32');  % height
        data.n_bscan  = fread(fid, 1, '*uint32');  % num_slices
        data.constant = fread(fid, 1, '*uint32');

        % Read B-scans        
        temp_file = 'temp.jpg';

        data.bscan = nan(data.n_axial, data.n_ascan, data.n_bscan);
        for i=1:data.n_bscan
            slice_size = fread(fid, 1, '*uint32');
            bscan_bytes = fread(fid, slice_size, '*uint8');

            % Ugly trick to read JPEG compressed info. It would be more
            % efficient to decode it without saving it.
            % Possible option
            % jImg = javax.imageio.ImageIO.read(java.io.ByteArrayInputStream(bscan_bytes));
            fid2 = fopen(temp_file, 'wb');
            fwrite(fid2, bscan_bytes, 'uint8'); 
            fclose(fid2);
            data.bscan(:,:,i) = imread(temp_file);              
        end
        delete(temp_file);       
        
    case '@IMG_MOT_COMP_03'
        data.unknown        = fread(fid, 1, '*uchar');
        data.n_ascan        = fread(fid, 1, '*uint32');
        data.n_axial        = fread(fid, 1, '*uint32');
        data.bits_per_pixel = fread(fid, 1, '*uint32');
        data.n_slices       = fread(fid, 1, '*uint32');
        data.unknown        = fread(fid, 1, '*uchar');
        data.size           = fread(fid, 1, '*uchar');
        
    case '@IMG_TRC_02'
        % Grayscale fundus image (lower quality)
        % It is stored twice
        % Probably not very useful
        width = fread(fid, 1, '*uint32');
        height = fread(fid, 1, '*uint32');
        bits_per_pixel = fread(fid, 1, '*uint32');
        num_slices = fread(fid, 1, '*uint32'); % 2 always
        val = fread(fid, 1, '*uint8');
        
        size = fread(fid, 1, '*uint32');        
        fundus_bytes = fread(fid, size, '*uint8');

        % Ugly trick to decode jp2
        temp_file = 'temp.jpg';
        fid2 = fopen(temp_file, 'wb');
        fwrite(fid2, fundus_bytes, 'uint8');             
        fclose(fid2);        
        fundus = imread(temp_file);          
        delete(temp_file);
        
        % No need to have 3 dimensions (it is grayscale not RGB)
        fundus = fundus(:,:,1);       
                
    case '@PARAM_SCAN_04'
        
        data.unknown  = fread(fid, 3, '*uint8');  % first is always three
        
        fixation      = fread(fid, 1, '*uint8');
        switch fixation
            case 0
                data.fixation = 'center';
            case 1
                data.fixation = 'onh';                
            case 2
                data.fixation = 'macula';                
            case 3
                data.fixation = 'wide';
            case 16
                data.fixation = 'external';
            otherwise
                data.fixation = 'unknown';
        end
        
        data.unknown = fread(fid, 8, '*uint8');  % first is always three

        data.size_x  = fread(fid, 1, '*float64');
        data.size_y  = fread(fid, 1, '*float64');
        data.scale_z = fread(fid, 1, '*float64');
        
        zeros        = fread(fid, 16, '*uint8');
        unknown      = fread(fid, 1,  '*float64');
        zeros        = fread(fid, 1,  '*uint8');        
        unknown      = fread(fid, 1,  '*uint16');
        data.n_bscan = fread(fid, 1,  '*uint16');
        
    case '@PATIENT_INFO_02'
        id          = char(fread(fid, 32, '*uint8'))';
        name        = char(fread(fid, 32, '*uint8'))';
        surname     = char(fread(fid, 32, '*uint8'))';
        zeros1      = fread(fid, 8, '*uint8');
        
        date_flag   = fread(fid, 1, '*uint8'); % gender?
        
        birth_year  = fread(fid, 1, '*uint16');
        birth_month = fread(fid, 1, '*uint16');
        birth_day   = fread(fid, 1, '*uint16');
        
        zeros2  = fread(fid, 504, '*uint16');        
        
    case '@PATIENT_INFO_03'
        % Apparently encoded to hide patient information
        
        unknown  = fread(fid, 2, '*uint16');
        id       = fread(fid, 1, '*uint16');
                    
    otherwise
        warning(['Chunk:' chunk_name ' is unknown. Unable to read']);
end