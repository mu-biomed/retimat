% function [header, segment, bscan, slo] = read_e2e(file, varargin)
%read_vol Read .e2e file exported from Spectralis OCT (Heidelberg Engineering)
%
%   [header, segment, bscan, slo] = read_e2e(file, options)
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

file = '/home/david/GITHUB/retimat/data_private/oct_1.e2e';
% file = '/home/david/GITHUB/retimat/data_private/octa_2.e2e';

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

% List of chunks
fseek(fid, current, -1);

magic3       = string(fread(fid, 12, '*char')');
version3     = fread(fid, 1, '*uint32');
unknown6     = fread(fid, 9, '*uint16');
unknown7     = fread(fid, 1, '*uint16');
num_entries2 = fread(fid, 1, '*uint32');
unknown8     = fread(fid, 1, '*uint32');
prev         = fread(fid, 1, '*uint32');
unknown9     = fread(fid, 1, '*uint32');

chunks = struct;

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
    
    chunks(i).type = type;
    chunks(i).pos = pos;
    chunks(i).start = start;
    chunks(i).size = size;
end

chunks = struct2table(chunks);

idx = find(chunks.type == 3);

for i=1:length(idx)
    fseek(fid, chunks.start(idx(i)), -1);
    laterality = fread(fid, chunks.size(idx(i)), '*uchar');
    disp(char(laterality)');
    
    
end
parse_chunk(fid, 11);


function data = parse_chunk(fid, type)

switch type
    case 3
        text = fread(fid, 12, '*char');
        
    case 10019  % segmentation
        unknown = fread(fid, 1, '*uint32');
        
    case 11
        unknown    = fread(fid, 14, '*uint8');
        laterality = fread(fid, 1, '*uint8');
        unknown    = fread(fid, 14, '*uint8');
        
    case 1073741824  % image data
        size    = fread(fid, 1, '*int32');
        type    = fread(fid, 1, '*int32');
        unknown = fread(fid, 1, '*int32');
        width   = fread(fid, 1, '*int32');
        height  = fread(fid, 1, '*int32');
        
        
    otherwise
        error("Unknown chunk type");
end
end