function [header, seg, bscan, fundus] = read_fda(file)
%FUNCTIONNAME Summary of this function goes here
%
%   Usage example OUT = template(IN1)
%   Detail explanation goes here
%
%   Input arguments:
%  
%   'ARG1'           Description of the argument. Type and purpose.          
%  
%                    Accepted values
%
%                    Default: 
%            
%  
%  
%   Output arguments:
%  
%   'ARG1'           Description of the argument. Type and purpose.          
%  
%
%   
%   Notes
%   -----
%   Important usage informationAnother name for a gray-level co-occurrence matrix is a gray-level
%   spatial dependence matrix.
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
%   Example
%   ---------      
%   % Read fda file
%
%     [GLCMS,SI] = graycomatrix(I,'NumLevels',9,'G',[])
%     
%
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2022

% close all;clc;clearvars;
% file = '../data_private/149989.fda'; % star
% file = '../data_private/149990.fda'; % raster

fid = fopen(file);

% Preliminary values (not very useful for analysis)
type     = fread(fid, 4, '*char')';
fixation = fread(fid, 3, '*char')';
unknown1 = fread(fid, 1, '*uint32');
unknown2 = fread(fid, 1, '*uint32');

% Discover all chunks (name-size-position)
chunks = struct;
i = 1;

more_chunks = true;
while more_chunks
    chunk_name_size = fread(fid,1,'*uint8');
    
    if chunk_name_size == 0
        break; 
    else
       chunks(i).name = char(fread(fid, chunk_name_size, '*uint8')'); 
       disp(chunks(i).name);
       chunks(i).size = fread(fid, 1, '*uint32'); % Bytes!
       chunks(i).pos = ftell(fid);       
       fseek(fid, chunks(i).pos + chunks(i).size, 'bof');
       i = i + 1;
    end    
end
chunks = struct2table(chunks);

% Read header-related things first
header = read_header(fid, chunks);
if nargout == 1
    return
end

% Read segmentation
seg = read_segmentation(fid, chunks);
if nargout == 2
    return
end

% Read B-Scan
bscan = read_bscan(fid, chunks);
if nargout == 3
    return
end

% Read fundus
fundus = read_fundus(fid, chunks);
fclose(fid);

function header = read_header(fid, chunks)

data = read_chunk(fid, chunks,'@CAPTURE_INFO_02');
header.eye = data.eye;
header.scan_date = datetime(data.year, data.month, data.day, ...
                            data.hour,data.minute,data.second);

data = read_chunk(fid, chunks,'@FDA_FILE_INFO');
A = 2;
       

function seg = read_segmentation(fid, chunks)

data = read_chunk(fid, chunks,'@CONTOUR_INFO');
seg = data.seg;
        
function bscan = read_bscan(fid, chunks)

data = read_chunk(fid, chunks,'@IMG_JPEG');
bscan = data.bscan;

function data = read_chunk(fid, chunks, chunk_name)

% Locate chunk
idx = find(strcmp(chunks.name, chunk_name));
fseek(fid, chunks.pos(idx), 'bof');

data = struct();
switch chunk_name
    case '@ALIGN_INFO'

    case '@ANTERIOR_CALIB_INFO'

    case '@CAPTURE_INFO_02'
        eye = fread(fid, 1, '*uint8');
        if eye==0
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
        data.size    = fread(fid, 1, '*uint32'); % useful?
        
        n_voxel = data.n_ascan * data.n_bscan;
        if type == 0
            seg = fread(fid, n_voxel, '*uint16');
        else
            seg = fread(fid, n_voxel, '*float64');
        end
                
        data.seg = reshape(seg, n_ascan, n_bscan); % reshape it        
        data.seg_version = fread(fid, 32, '*char')';
        
    case '@CONTOUR_MASK_INFO'
        
    case '@EFFECTIVE_SCAN_RANGE'
        
    case '@FAST_Q2_INFO'
            
    case '@FDA_DISC_SEGMENTATION'
        
    case '@FDA_FILE_INFO'
        data.info1   = fread(fid, 1, '*uint32');
        data.info2   = fread(fid, 1, '*uint32');
        data.version = fread(fid, 32, '*char')';        

    case '@FOCUS_MOTOR_INFO'
        
    case '@GLA_LITTMANN_01'
        
    case '@HW_INFO_03'
        data.model_name = char(fread(fid, 16, '*uint8')');        
        data.serial_number = fread(fid, 16, '*uint8');        
        data.zeros = fread(fid, 32, '*uint8');        
        data.version = fread(fid, 16, '*uint8');        
        data.build_year = fread(fid, 1, '*uint16');        
        data.build_month = fread(fid, 1, '*uint16');        
        data.build_day = fread(fid, 1, '*uint16');        
        data.build_hour = fread(fid, 1, '*uint16');        
        data.build_minute = fread(fid, 1, '*uint16');        
        data.build_second = fread(fid, 1, '*uint16');        
        data.zeros2 = fread(fid, 8, '*uint8');        
        data.version_numbers = fread(fid, 5*16, '*uint8');        
        
    case '@IMG_EN_FACE'
        
    case '@IMG_FUNDUS'
        data.width          = fread(fid, 1, '*uint32');
        data.height         = fread(fid, 1, '*uint32');
        data.bits_per_pixel = fread(fid, 1, '*uint32');
        data.num_slices     = fread(fid, 1, '*uint32');
        data.val            = fread(fid, 1, '*uint32');        
        data.size           = fread(fid, 1, '*uint32');        
        data.fundus_bytes   = fread(fid, size, '*uint8');

        % Ugly trick to decode jp2
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
        data.type     = fread(fid, 1, '*uint8');
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

    case '@IMG_PROJECTION'

    case '@IMG_TRC_02'
        % Grayscale fundus image (lower quality)
        % It is store twice
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
        
    case '@MAIN_MODULE_INFO'

    case '@PARAM_EN_FACE'

    case '@PARAM_OBS_02'
        
    case '@PARAM_SCAN_04'
        data.unknown = fread(fid, 6, '*uint16');
        
        data.size_x = fread(fid, 1, '*float64');
        data.size_y = fread(fid, 1, '*float64');
        data.scale_z = fread(fid, 1, '*float64');
%         scale_x = fread(fid, 4, '*float64');
        
%         unknown2 = fread(fid, 2, '*float64');
%         unknown3 = fread(fid, 2, '*float64');
        
    case '@PARAM_TRC'
        
    case '@PARAM_TRC_02'        

    case '@PATIENTEXT_INFO'
        
    case '@PATIENT_INFO_02'
        id          = char(fread(fid, 32, '*uint8'))';
        name        = char(fread(fid, 32, '*uint8'))';
        surname     = char(fread(fid, 32, '*uint8'))';
        zeros1        = fread(fid, 8, '*uint8');
        
        date_flag   = fread(fid, 1, '*uint8'); % gender?
        
        birth_year  = fread(fid, 1, '*uint16');
        birth_month = fread(fid, 1, '*uint16');
        birth_day   = fread(fid, 1, '*uint16');
        
        zeros2  = fread(fid, 504, '*uint16');
        
        
    case '@PATIENT_INFO_03'
        % Apparently encoded to hide patient information
        
        unknown  = fread(fid, 2, '*uint16');
        id       = fread(fid, 1, '*uint16');
            
    case '@REF_IMG_SCAN'
        
    case '@REGIST_INFO'
               
    case '@REPORT_INFO'
    
    case '@RESULT_CORNEA_CURVE'

    case '@SCAN_POS_COMP_DATA'
        
    case '@THUMBNAIL'
        
    case '@TOPQEXT_INFO'
        
    otherwise
        warning(['Chunk:' chunk_name ' is unknown. Unable to read']);
end