function [header, seg, bscan, fundus] = read_e2e(file, varargin)
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
%                    'verbose': Display header info during read.
%
%
%   Output arguments:
%  
%   'header'         Structure with .vol file header values.          
%  
%   'seg'            Segmenation data stored in the .vol file.
%
%   'bscan'          3D single image with B-Scans.
%
%   'fundus'         2D fundus image.
%   
%
%   Notes
%   -----
%   This function was developed based on the previous reverse engineering 
%   attempts [1-3] and is not an official file reader. Therefore, some of 
%   the information retrieved can be incorrect/incomplete.
%   
%   Scan focus, scale_x, scale_y and acquisition pattern not found yet.
%
%   Spectralis OCT data can be exported into both E2E and vol format. We
%   recommend using the latter as it is easier to parse and it only stores
%   a single eye and acquisition. 
%
%   References
%   ----------
%   [1] uocte documentation
%   https://bitbucket.org/uocte/uocte/wiki/Topcon%20File%20Format
%
%   [2] OCT-Converter, https://github.com/marksgraham/OCT-Converter
%   
%   [3] LibE2E, https://github.com/neurodial/LibE2E
%   
%   Examples
%   ---------      
%   % Read all the information in a .e2e/.E2E file
%
%     file = 'my_oct.e2e';
%     [header, segment, bscan, fundus] = read_e2e(file)
%     
%
%   % Read only the header (faster) of the .e2e/.E2E file
%     file = 'my_oct.e2e';
%     header = read_e2e(file)
%
%
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2022

verbose = any(strcmp('verbose', varargin));

global SEG_FLAG IMAGE_FLAG NA_FLAG BSCAN_META_FLAG PATIENT_FLAG EYE_FLAG 
PATIENT_FLAG    = 9;
EYE_FLAG        = 11;
BSCAN_META_FLAG = 10004;
SEG_FLAG        = 10019;
IMAGE_FLAG      = 1073741824;
NA_FLAG         = 4294967295; % = 2^32 - 1 (all ones, means not applicable)
CHUNK_HEADER_SIZE = 60;

fid = fopen(file, 'rb', 'l');
 
chunks = discover_chunks(fid, verbose);

% Get number of patients 
patients = unique(chunks.patient_id);
patients(patients == NA_FLAG) = [];
n_patient = length(patients);

if n_patient > 1
    error(['Number of patients in file is > 1 (' num2str(n_patient) ')']);        
end

% Get number of series per subject
series_id = unique(chunks.series_id);
series_id(series_id == NA_FLAG) = [];
n_series = length(series_id);

if verbose
    disp(['Number of series (images) in file: ' num2str(n_series)]);
end

% Read patient data first
patient_header = read_patient_data(fid, chunks);

% Parse each series separately
header = cell(1, n_series);
bscan  = cell(1, n_series);
seg    = cell(1, n_series);
fundus = cell(1, n_series);

for i_series = 1:n_series
    chunks_series = chunks(chunks.series_id == series_id(i_series),:);

    header{i_series} = read_header(fid, chunks_series, patient_header);
    
    bscan{i_series}  = read_bscan(fid, chunks_series, verbose);

    seg{i_series}    = read_segmentation(fid, chunks_series);
     
    fundus{i_series} = read_fundus(fid, chunks_series);      
end

function chunks = discover_chunks(fid, verbose)
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
        bscan_id   = fread(fid, 1, '*uint32');
        unknown    = fread(fid, 1, '*uint16');
        zero2      = fread(fid, 1, '*uint16');
        type       = fread(fid, 1, '*uint32');
        unknown11  = fread(fid, 1, '*uint32');

        if size > 0            
            chunks(i_chunk).patient_id = patient_id;
            chunks(i_chunk).series_id  = series_id;
            chunks(i_chunk).bscan_id   = bscan_id;
            chunks(i_chunk).type       = type;
            chunks(i_chunk).pos        = pos;
            chunks(i_chunk).start      = start;
            chunks(i_chunk).size       = size;

            i_chunk = i_chunk + 1;
        end
    end

    if verbose
        disp(['Read element ' num2str(i_main)]);
    end
    
    i_main = i_main + 1;
end

chunks = struct2table(chunks);

function header = read_patient_data(fid, chunks)
global PATIENT_FLAG

start = chunks.start(chunks.type == PATIENT_FLAG);
fseek(fid, start, -1);
header = parse_chunk(fid, PATIENT_FLAG);

function header = read_header(fid, chunks, header)
global EYE_FLAG BSCAN_META_FLAG

start = chunks.start(chunks.type == EYE_FLAG);
fseek(fid, start(1), -1);
data = parse_chunk(fid, EYE_FLAG);

header.eye = data.eye;

bscan_id = find(chunks.type == BSCAN_META_FLAG);
[~, idx] = sort(chunks.bscan_id(bscan_id));
bscan_id = bscan_id(idx);

for ii=1:length(bscan_id)
    start = chunks.start(bscan_id(ii));
    fseek(fid, start(1), -1);
    data = parse_chunk(fid, BSCAN_META_FLAG);

    if ii == 1
        header.n_bscan   = data.n_bscan;
        header.n_ascan   = data.n_ascan;
        header.n_axial   = data.n_axial;
        header.scale_z   = data.scale_z;
        header.fixation  = data.fixation;
        header.scan_date = data.scan_date; % we can have it for each bscan
    end
    
    header.quality(ii)   = data.quality;
    header.n_average(ii) = data.n_average;
end

function bscan = read_bscan(fid, chunks, verbose)
global IMAGE_FLAG NA_FLAG

is_image = chunks.type == IMAGE_FLAG;
is_bscan = chunks.bscan_id ~= NA_FLAG;
     
bscan_id = chunks.bscan_id(is_image & is_bscan);     
n_bscan = length(bscan_id);

if verbose
    disp(['Reading ' num2str(n_bscan) ' bscans']);
end
     
bscan_id = sort(bscan_id);
     
for i_bscan=1:n_bscan
    is_bscan = chunks.bscan_id == bscan_id(i_bscan);

    start = chunks.start(is_bscan & is_image);
    fseek(fid, start, -1);

    data = parse_chunk(fid, IMAGE_FLAG);

    % This can be eliminated if the header is built before
    if i_bscan == 1
       bscan = nan(data.n_axial, data.n_ascan, n_bscan);             
    end

    bscan(:, :, i_bscan) = data.bscan;
end  
    
function fundus = read_fundus(fid, chunks)
global IMAGE_FLAG NA_FLAG

is_image     = chunks.type == IMAGE_FLAG;
is_not_bscan = chunks.bscan_id == NA_FLAG;

start = chunks.start(is_image & is_not_bscan);
fseek(fid, start, -1);

data   = parse_chunk(fid, IMAGE_FLAG);
fundus = data.fundus;
     
function seg = read_segmentation(fid, chunks)
global SEG_FLAG NA_FLAG
layer_names = {'ILM',     'BM',      'RNFL_GCL', 'GCL_IPL', ...
               'IPL_INL', 'INL_OPL', 'OPL_ONL', 'unknown', 'ELM', ...
               'unknown', 'unknown', 'unknown','unknown','unknown', ...
               'MZ_EZ',   'OSP_IZ',  'IZ_RPE'};
           
is_seg = chunks.type == SEG_FLAG;

bscan_id = unique(chunks.bscan_id);
bscan_id(bscan_id == NA_FLAG) = [];
n_bscan  = length(bscan_id);

seg = struct;
for i_bscan=1:n_bscan
    is_bscan = chunks.bscan_id == bscan_id(i_bscan);
    
    start = chunks.start(is_seg & is_bscan);
    n_layer = length(start);
    for i_layer=1:n_layer
        fseek(fid, start(i_layer), -1);
        data = parse_chunk(fid, SEG_FLAG);        
        
        layer = layer_names{data.layer_id + 1};
        seg.(layer)(i_bscan, :) = data.seg;
    end
end

% Invalid values coded as 0 --> better have nan to avoid wrong computations
layer_names = fields(seg);
for i_layer=1:length(layer_names)
    layer = layer_names{i_layer};
    seg.(layer)(seg.(layer) == 0) = nan;
end

function data = parse_chunk(fid, type)
global SEG_FLAG IMAGE_FLAG PATIENT_FLAG BSCAN_META_FLAG EYE_FLAG

% We can probably skip chunk header as that info is not really useful
magic4     = string(fread(fid, 12, '*char')');
unknown    = fread(fid, 2, '*uint32');
pos        = fread(fid, 1, '*uint32');
c_size     = fread(fid, 1, '*uint32');
zero       = fread(fid, 1, '*uint32');
patient_id = fread(fid, 1, '*int32');
study_id   = fread(fid, 1, '*int32');
series_id  = fread(fid, 1, '*int32');
bscan_id   = fread(fid, 1, '*uint32');
ind        = fread(fid, 1, '*uint16');
unknown2   = fread(fid, 1, '*uint16');
c_type     = fread(fid, 1, '*uint32');
unknown3   = fread(fid, 1, '*uint32');

switch type
    case 3
        text = fread(fid, 12, '*char');
        
    case 7
        eye = fread(fid, 1, '*char');
        
    case PATIENT_FLAG % patient info
        name       = deblank(fread(fid, 31, '*char')');
        surname    = deblank(fread(fid, 66, '*char')');
        birth_date = fread(fid, 1, '*uint32');
        sex        = fread(fid, 1, '*char');
        
        birth_date = (birth_date/64) - 14558805;  % to julian date
        birth_date = datetime(birth_date, 'ConvertFrom', 'juliandate');
        
        data.name       = name;
        data.surname    = surname;
        data.birth_date = datestr(birth_date);
        
        if sex == 'M'
            data.sex = 'male';
        elseif sex == 'F'
            data.sex = 'female';
        else
            data.sex = 'unknown';
        end        
        
    case EYE_FLAG
        unknown = fread(fid, 14, '*char');
        eye     = fread(fid, 1, '*char');  % Laterality
        unknown = fread(fid, 14, '*uint8');
  
        if eye == 'L'
            data.eye = 'OS';
        elseif eye == 'R'
            data.eye = 'OD';
        else
            warning("Unknown laterality value");
            data.eye = 'unknown';
        end
    case 13
        device = fread(fid, 260, '*char')';
        
    case 52
        code = fread(fid, 97, '*char');
    case 53
        code = fread(fid, 97, '*char');        
        
    case BSCAN_META_FLAG
        unknown1 = fread(fid, 1, '*uint32'); %% always 20(20 kHz scan speed maybe)
        
        n_axial   = fread(fid, 1, '*uint32');
        n_ascan1  = fread(fid, 1, '*uint32');        
        pos_x1    = fread(fid, 1, '*float');
        pos_y1    = fread(fid, 1, '*float');
        pos_x2    = fread(fid, 1, '*float');
        pos_y2    = fread(fid, 1, '*float');

        zero1     = fread(fid, 1, '*uint32');
        unknown2  = fread(fid, 1, '*float');  %
        
        scale_z   = fread(fid, 1, '*float');
        
        unknown3  = fread(fid, 1, '*float');
        zero2     = fread(fid, 1, '*uint32');        
        unknown4  = fread(fid, 2, '*float');
        zero3     = fread(fid, 1, '*uint32');
        
        n_ascan2  = fread(fid, 1, '*uint32');
        n_bscan   = fread(fid, 1, '*uint32');
        bscan_id  = fread(fid, 1, '*uint32');
        fixation  = fread(fid, 1, '*uint32');

        center_x  = fread(fid, 1, '*float');
        center_y  = fread(fid, 1, '*float');

        unknown5  = fread(fid, 1, '*uint32');
        scan_date = fread(fid, 1, '*int64'); %% exam_time in read_vol()
        unknown6  = fread(fid, 6, '*uint32');
        n_average = fread(fid, 1, '*uint32');
        unknown7  = fread(fid, 8, '*uint32');

        quality   = fread(fid, 1, '*single');
        
        data.n_ascan   = n_ascan1;
        data.n_bscan   = n_bscan;
        data.n_axial   = n_axial;        
        data.scale_z   = scale_z;
        data.quality   = quality;
        data.n_average = n_average;
        
        % Transform date to dates and add offset since 1600/12/31 (unsure)
        data.scan_date = datestr(double(scan_date)/(1e7*60*60*24) + 584755 + 2/24); 
        
        if fixation == 1
            data.fixation = 'macula';
        elseif fixation == 2
            data.fixation = 'disk';
            if data.n_bscan == 0  % don't know why peripapillar are saved as 0
                data.n_bscan = 1;
            end
        else
            data.fixation = 'unknown';
        end
                
    case SEG_FLAG  % segmentation
        unknown  = fread(fid, 1, '*uint32');
        layer_id = fread(fid, 1, '*uint32');
        unknown  = fread(fid, 1, '*uint32');
        n_ascan  = fread(fid, 1, '*uint32');
        
        seg      = fread(fid, n_ascan, '*float32');
        
        data.layer_id = layer_id;
        data.n_ascan  = n_ascan;
        data.seg      = seg;

    case IMAGE_FLAG  % image data
        data_size  = fread(fid, 1, '*int32');
        image_type = fread(fid, 1, '*int32');
        n_pixel    = fread(fid, 1, '*int32');
        
        switch image_type
            case 33620481 % fundus
                width  = fread(fid, 1, '*int32');
                height = fread(fid, 1, '*int32');
                bytes  = fread(fid, n_pixel, '*uint8');
                
                fundus = reshape(bytes, [width height]);
                fundus = permute(fundus, [2 1]);
                
                data.fundus     = fundus;
                data.n_x_fundus = width;
                data.n_y_fundus = height;
            case 35652097 % b-scan 
                n_axial = fread(fid, 1, '*int32');
                n_ascan = fread(fid, 1, '*int32');
                bytes   = double(fread(fid, n_pixel, '*uint16'));
                
                % Old implementation (intuitive but very slow)
                %  bin      = dec2bin(bytes);
                %  exponent = bin2dec(bin(:, 1:6));
                %  mantissa = bin2dec(bin(:, 7:end));
                
                exponent = floor(bytes / 2^10);
                mantissa = mod(bytes, 2^10);

                a     = 1 + mantissa/2^10;
                b     = 2 .^ (exponent - 63);
                bscan = a .* b;
                
                bscan = reshape(bscan, [n_ascan n_axial]);
                bscan = permute(bscan, [2 1]);

                data.n_ascan = n_ascan;
                data.n_axial = n_axial;
                data.bscan   = bscan;
            otherwise
                warning('Unknown image type');
        end
    case 1073751825 % looks like a time series
%         unknown = fread(fid, 300,'*uint8');
        
    case 1073751826
        unknown = fread(fid, 4, '*uint32');
        nx_fundus = fread(fid, 1, '*uint32');
        
    otherwise
        error("Unknown chunk type");
end
