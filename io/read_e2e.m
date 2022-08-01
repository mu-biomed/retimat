% function [header, segment, bscan, slo] = read_e2e(file, varargin)
%read_e2e Read .e2e file exported from Spectralis OCT (Heidelberg Engineering)
%
%   [header, segment, bscan, slo] = read_e2e(file, options)
%
%   This function reads the header, segmentation and image information 
%   contained in the .e2e/.E2E files. 
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
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2022

close all;clc;clearvars;

% file = '/home/david/GITHUB/retimat/data_private/oct_1.e2e';
file = 'C:/Users/dromero/Desktop/GITHUB/retimat/data_private/oct_1.e2e';

fid = fopen(file, 'rb', 'l');
 
% Read header
magic1   = string(fread(fid, 12, '*char')');
version1 = fread(fid, 1, '*uint32');
unknown1 = fread(fid, 9, '*uint16');
unknown2 = fread(fid, 1, '*uint16');

% Directory
magic2      = string(fread(fid, 12, '*char')');
version2    = fread(fid, 1, '*uint32');
unknown3    = fread(fid, 9, '*uint16');
unknown4    = fread(fid, 1, '*uint16');
num_entries = fread(fid, 1, '*uint32');
current     = fread(fid, 1, '*uint32');
zeros2      = fread(fid, 1, '*uint32');
unknown5    = fread(fid, 1, '*uint32');

chunks = struct;
i_main = 1;
i_chunk = 1;
while current ~=0% List of chunks
    fseek(fid, current, -1);

    magic3       = string(fread(fid, 12, '*char')');
    version3     = fread(fid, 1, '*uint32');
    unknown6     = fread(fid, 9, '*uint16');
    unknown7     = fread(fid, 1, '*uint16');
    num_entries2 = fread(fid, 1, '*uint32');
    unknown8     = fread(fid, 1, '*uint32');
    prev         = fread(fid, 1, '*uint32');
    unknown9     = fread(fid, 1, '*uint32');

    current = prev;

    for i=1:num_entries2
        pos        = fread(fid, 1, '*uint32');
        start      = fread(fid, 1, '*uint32');
        size       = fread(fid, 1, '*uint32');
        zero       = fread(fid, 1, '*uint32');
        patient_id = fread(fid, 1, '*uint32');
        study_id   = fread(fid, 1, '*uint32');
        series_id  = fread(fid, 1, '*uint32');
        slice_id   = fread(fid, 1, '*uint32');
        unknown    = fread(fid, 1, '*uint16');
        zero2      = fread(fid, 1, '*uint16');
        type       = fread(fid, 1, '*uint32');
        unknown11  = fread(fid, 1, '*uint32');

        if size > 0
            chunks(i_chunk).type  = type;
            chunks(i_chunk).pos   = pos;
            chunks(i_chunk).start = start;
            chunks(i_chunk).size  = size;

            i_chunk = i_chunk + 1;
        end
    end

    disp(['Read element ' num2str(i_main)]);
    
    i_main = i_main + 1;
end

chunks = struct2table(chunks);

chunk_header_size = 60;

match = search_by_value(fid, chunks, 49, {'*uint32'});
match = struct2table(match);

chunk_id = 512;
idx = find(chunks.type == chunk_id);

for i=1:length(idx)
    fseek(fid, chunks.start(idx(i)), -1);
    fread(fid, chunk_header_size, '*uint8');
    x = fread(fid, chunks.size(idx(i)), '*char')';
    disp(x);
    
    fseek(fid, chunks.start(idx(i)), -1);
    fread(fid, chunk_header_size, '*uint8');
    x8 = fread(fid, chunks.size(idx(i)), '*uint8')';
    x16 = fread(fid, chunks.size(idx(i))/2, '*uint16')';
    x32 = fread(fid, chunks.size(idx(i))/3, '*uint32')';
   %     disp(laterality');
    
%     parse_chunk(fid, chunk_id);

    
end


function match = search_by_value(fid, chunks, value, types)

match = struct;
i_match = 1;
for i=1:size(chunks,1)
    chunk_type = chunks.type(i);
    start   = chunks.start(i);
    n_bytes = chunks.size(i);

    disp(['Chunk ' num2str(i) '/' num2str(size(chunks,1))]);
    
    fseek(fid, start, -1);
    fread(fid, 60, '*uint8');
    
    x = {};
    types_all = {};
    for j=1:length(types)
        switch types{j}
            case '*uint8'
                x{end+1} = fread(fid, n_bytes, '*uint8');
                types_all{end+1} = '*uint8';
            case '*uint16'
                x{end+1} = fread(fid, n_bytes/2, '*uint16');
                
                fseek(fid, start, -1);
                fread(fid, 61, '*uint8');
                x{end+1} = fread(fid, n_bytes/2, '*uint16');
                
                types_all(end+1:end+2) = {'*uint16_0','*uint16_1'};
            case '*uint32'
                x{end+1} = fread(fid, n_bytes/4, '*int32');
                
                fseek(fid, start, -1);
                fread(fid, 61, '*uint8');                
                x{end+1} = fread(fid, n_bytes/4, '*int32');
                
                fseek(fid, start, -1);
                fread(fid, 62, '*uint8');                
                x{end+1} = fread(fid, n_bytes/4, '*int32');
                
                fseek(fid, start, -1);                
                fread(fid, 63, '*uint8');
                x{end+1} = fread(fid, n_bytes/4, '*int32');
                
                
                types_all(end+1:end+4) = {'*uint32_0','*uint32_1','*uint32_2','*uint32_3'};
        end                
    end
    
    for j=1:length(x)
        idx = find(x{j} == value);
        
        for ii=1:length(idx)
            match(i_match).idx = i;
            match(i_match).chunk_type = chunk_type;
            match(i_match).data_type = types_all{j};
            match(i_match).pos = idx(ii);
                        
            i_match = i_match + 1;
        end
    end
end
end

function data = parse_chunk(fid, type)

magic4     = string(fread(fid, 12, '*char')');
unknown    = fread(fid, 2, '*uint32');
pos        = fread(fid, 1, '*uint32');
c_size     = fread(fid, 1, '*uint32');
zero       = fread(fid, 1, '*uint32');
patient_id = fread(fid, 1, '*uint32');
study_id   = fread(fid, 1, '*uint32');
series_id  = fread(fid, 1, '*uint32');
slice_id   = fread(fid, 1, '*uint32');
ind        = fread(fid, 1, '*uint16');
unknown2   = fread(fid, 1, '*uint16');
c_type     = fread(fid, 1, '*uint32');
unknown3   = fread(fid, 1, '*uint32');

disp(series_id);

switch type
    case 3
        text = fread(fid, 12, '*char');
        
    case 7
        eye = fread(fid, 1, '*char');
        
    case 9 % patient info
        name       = deblank(string(fread(fid, 31, '*char')'));
        surname    = deblank(string(fread(fid, 66, '*char')'));
        birth_date = fread(fid, 1, '*uint32');
        sex        = fread(fid, 1, '*char');
        
        birth_date = (birth_date/64) - 14558805;  % to julian date
        birth_date = datetime(birth_date, 'ConvertFrom', 'juliandate');
                 
    case 11
        unknown    = fread(fid, 14, '*char');
        laterality = fread(fid, 1, '*char');
        unknown    = fread(fid, 14, '*uint8');
  
    case 13
        device = fread(fid, 260, '*char')';
        
    case 52
        code = fread(fid, 97, '*char');
    case 53
        code = fread(fid, 97, '*char');        
    case 10019  % segmentation
        unknown = fread(fid, 1, '*uint32');

    case 1073741824  % image data
        size       = fread(fid, 1, '*int32');
        image_type = fread(fid, 1, '*int32');
        n_pixel    = fread(fid, 1, '*int32');
        width      = fread(fid, 1, '*int32');
        height     = fread(fid, 1, '*int32');
        
        switch image_type
            case 33620481 % fundus
                bytes = fread(fid, n_pixel, '*uint8');
                I = reshape(bytes, [height width]);
                permute(I, [2 1]);
                imagesc(I);
                
            case 35652097 % b-scan                
                bytes = fread(fid, n_pixel, '*uint16');
                
                bin      = dec2bin(bytes);
                exponent = bin2dec(bin(:, 1:6));
                mantissa = bin2dec(bin(:, 7:end));
                a        = (1 + mantissa) / (2^10);
                b        = 2 .^ (exponent - 63);
                I        = a .* b;
                
                I = reshape(I, [height width]);
                I = permute(I, [2 1]);
                close all;
                imagesc(I.^0.25);
                disp('');
            otherwise
                warning('Unknown image type');
        end
    case  1073751825 % looks like a time series
%         unknown = fread(fid, 300,'*uint8');
        
    otherwise
        error("Unknown chunk type");
end
end