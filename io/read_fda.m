function [header, oct] = read_fda()
close all;clc;clearvars;
file = '149989.fda';

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
       chunk_size = fread(fid, 1, '*uint32');
       chunk_pos = ftell(fid);
       data = read_chunk(fid, chunk_name, chunk_pos);
       
       fseek(fid, chunk_pos + chunk_size, 'bof');
    end    
end
end

function data = read_chunk(fid, chunk_name, chunk_pos)

data = struct();
fseek(fid, chunk_pos, 'bof');

switch chunk_name
    case '@FDA_FILE_INFO'
        data.info1 = fread(fid, 1, '*uint32');
        data.info2 = fread(fid, 1, '*uint32');
        data.info3 = fread(fid, 32, '*uint8');        
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
    case '@PATIENT_INFO_02'
        
    case '@CAPTURE_INFO_02'
        
    case '@IMG_JPEG' 
        data.type = fread(fid, 1, '*uint8');
        data.unknown1 = fread(fid, 1, '*uint32');
        data.unknown2 = fread(fid, 1, '*uint32');
        data.width = fread(fid, 1, '*uint32');
        data.height = fread(fid, 1, '*uint32');
        data.num_slices = fread(fid, 1, '*uint32');
        data.constant = fread(fid, 1, '*uint32');
        
        for i=1:data.num_slices
            slice_size = fread(fid, 1, '*uint32');
            bscan_bytes = fread(fid, slice_size, '*uint8');
            
%             jImg = javax.imageio.ImageIO.read(java.io.ByteArrayInputStream(bscan_bytes));

            fid = fopen('buffer.jpg', 'wb');
            fwrite(fid, bscan_bytes, 'uint8'); 
            fclose(fid);
            I = imread('buffer.jpg');              
%             imshow(I);
        end
        
    case '@PARAM_SCAN_04'
        
    case '@IMG_TRC_02'
        
    case '@PARAM_TRC'
        
    case '@IMG_FUNDUS'
        
    case '@PARAM_OBS_02'
        
    case '@IMG_MOT_COMP_03'
        
    case '@IMG_PROJECTION'
        
    case '@REGIST_INFO'
        
    case '@ALIGN_INFO'
        
    case '@CONTOUR_INFO'
        
    case '@FAST_Q2_INFO'
        
    case '@THUMBNAIL'
        
    case '@GLA_LITTMANN_01'
        
    case '@PATIENTEXT_INFO'
        
    case '@EFFECTIVE_SCAN_RANGE'
        
    case '@MAIN_MODULE'
        
    case '@RESULT_CORNEA_CURVE'
        
    case '@ANTERIOR_CALIB_INFO'
        
    case '@REPORT_INFO'
        
    otherwise
        warning(['Chunk:' chunk_name ' is unknown. Unable to read']);
end

end