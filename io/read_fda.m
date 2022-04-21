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

close all;clc;clearvars;
file = '../data_private/149989.fda'; % star
file = '../data_private/149990.fda'; % raster


fid = fopen(file);

type = char(fread(fid, 4, '*uchar')');
fixation = char(fread(fid, 3, '*uchar')');
version1 = fread(fid,1,'*uint32');
version2 = fread(fid,1,'*uint32');

more_chunks = true;

while more_chunks
    chunk_name_size = fread(fid,1,'*uint8');
    
    if chunk_name_size == 0
        break; 
    else
       chunk_name = char(fread(fid, chunk_name_size, '*uint8')'); 
       disp(chunk_name);
       chunk_size = fread(fid, 1, '*uint32'); % Bytes!
       chunk_pos = ftell(fid);
       data = read_chunk(fid, chunk_name, chunk_pos, chunk_size);
       
       fseek(fid, chunk_pos + chunk_size, 'bof');
    end    
end
end

function data = read_chunk(fid, chunk_name, chunk_pos, chunk_size)

data = struct();
fseek(fid, chunk_pos, 'bof');

switch chunk_name
    case '@ALIGN_INFO'

    case '@ANTERIOR_CALIB_INFO'

    case '@CAPTURE_INFO_02'
        eye = fread(fid, 1, '*uint8');
        
        unknown = fread(fid,1,'*uint8');
        zero = fread(fid, 52, '*uint16');
        
        year = fread(fid, 1, '*uint16');
        month = fread(fid, 1, '*uint16');
        day = fread(fid, 1, '*uint16');
        hour = fread(fid, 1, '*uint16');
        minute = fread(fid, 1, '*uint16');
        second = fread(fid, 1, '*uint16');
        
    case '@CONTOUR_INFO'
        % Segmentation
        id = char(fread(fid, 20, '*uint8')');
        
        type = fread(fid, 1, '*uint16');
        n_ascan = fread(fid, 1, '*uint32');  % width
        n_bscan = fread(fid, 1, '*uint32');  % height
        size = fread(fid, 1, '*uint32'); % useful?
        
        if type == 0
            seg = fread(fid, n_ascan*n_bscan, '*uint16');
        else
            seg = fread(fid, n_ascan*n_bscan, '*float64');
        end
        
        % Reshape it
        seg = reshape(seg, n_ascan, n_bscan);
        
        seg_version = char(fread(fid, 32, '*uint8')');
        
    case '@CONTOUR_MASK_INFO'
        
    case '@EFFECTIVE_SCAN_RANGE'
        
    case '@FAST_Q2_INFO'
            
    case '@FDA_DISC_SEGMENTATION'
        
    case '@FDA_FILE_INFO'
        data.info1 = fread(fid, 1, '*uint32');
        data.info2 = fread(fid, 1, '*uint32');
        data.version = char(fread(fid, 32, '*uint8'))';        

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
        width = fread(fid, 1, '*uint32');
        height = fread(fid, 1, '*uint32');
        bits_per_pixel = fread(fid, 1, '*uint32');
        num_slices = fread(fid, 1, '*uint32');
        val = fread(fid, 1, '*uint32');
        
        size = fread(fid, 1, '*uint32');
        
        fundus_bytes = fread(fid, size, '*uint8');

        % Ugly trick to decode jp2
        temp_file = 'temp.jpg';
        fid2 = fopen(temp_file, 'wb');
        fwrite(fid2, fundus_bytes, 'uint8');             
        fclose(fid2);        
        fundus = imread(temp_file);          
        delete(temp_file);
        
        % Reorient fundus properly
        fundus = flip(fundus, 1); % Upside-down
        fundus = flip(fundus, 3); % BGR to RGB
        
    case '@IMG_JPEG' 
        return
        data.type = fread(fid, 1, '*uint8');
        data.unknown1 = fread(fid, 1, '*uint32');
        data.unknown2 = fread(fid, 1, '*uint32');
        data.n_ascan = fread(fid, 1, '*uint32');  % width
        data.n_axial = fread(fid, 1, '*uint32');  % height
        data.n_bscan = fread(fid, 1, '*uint32');  % num_slices
        data.constant = fread(fid, 1, '*uint32');
        
        temp_file = 'temp.jpg';
        
        bscan = nan(data.n_axial, data.n_ascan, data.n_bscan);
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
            bscan(:,:,i) = imread(temp_file);              
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
        unknown = fread(fid, 6, '*uint16');
        
        size_x = fread(fid, 1, '*float64');
        size_y = fread(fid, 1, '*float64');
        scale_z = fread(fid, 1, '*float64');
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

end