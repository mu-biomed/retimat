function [header, seg, bscan, fundus] = read_tabs(folder, coordinates)
%READ_TABS Read TABS output
%
%   [header, seg, bscan, slo] = read_tabs(folder)
%   Load and organize the output of Topcon Advanced Boundary Segmentation
%   (TABS). All TABS output files must be in the in_folder folder.
%
%   Input arguments:
%  
%   'in_folder'      String with the directory containing TABS outputs.
%            
%  
%   Output arguments:
%  
%   'header'         Structure with metadata.
%
%   'seg'            Struct with the segmentation of retinal layers 
%                    performed by TABS.
%
%   'bscan'          3D matrix containing bscans.
%
%   'fundus'         Fundus image.
%  
%  
%   Notes
%   -----
%   This function has only been tested with a macular scans and may fail
%   with other acquisition patterns.
%
%
%   Example
%   ---------      
%   % Read TABS output
%
%     [header, seg, bscan, slo] = read_tabs('../tabs_out_folder')
%
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2022.


% Check arguments
if nargin == 1
    coordinates = false;
end

% Read meta-data
file = [folder '/FDSInfo.txt'];
header = read_metadata(file, coordinates);

% Read extra meta-data (foveal center)
file = [folder '/SegIndicators.txt'];
header = read_metadata2(file, header, coordinates);
if nargout == 1
    return;
end

% Read segmented boundary names
file = [folder '/TabilSegStat000.txt'];
header = read_seg_layers(file, header);

% Read segmentation
file = [folder '/TabilSegStat.dat'];
seg = read_seg(file, header);
if nargout == 2
    return; 
end

% Read bscans
file = [folder '/AlignedOCTImages.dat'];
bscan = read_bscan(file, header);
if nargout == 3 
    return
end

% Read fundus
file = [folder '/ColorFundus.dat'];
fundus = read_fundus(file, header);

function header = read_metadata(file, coordinates)
% Parameters are structured in categories and fields within each category.
% This function merges all the fields into a single header struct.

if ~exist(file,'file')
    error(['Unable to find metadata file: ' file]);
end

fid = fopen(file);

% Read all the text and convert it to lines
meta_data = char(fread(fid)');
lines = splitlines(meta_data);

% Remove empty lines (last one basically)
is_empty = cellfun(@isempty, lines);
lines(is_empty) = [];
n_line = length(lines);

% Loop through fields
header = struct;
for i_line=1:n_line
    
    line = lines{i_line};
    
    % If category
    if line(1) ~= ' '
        category = line;
        continue;
    end
    
    % If field
    idx = strfind(line,':');
    if isempty(idx)
        idx = strfind(line,'=');            
    end

    field_name = strtrim(line(1:idx(1)-1));        
    field_val = strtrim(line(idx(1)+1:end));
                
    switch field_name
        case 'Subject ID'
            header.subject_id = field_val;
        case 'Last name'
            header.last_name = field_val;                        
        case 'First name'
            header.first_name = field_val;
        case 'Gender'
            header.sex = field_val;                        
        case 'Date of Birth'
            header.date_of_birth = field_val;                        
        case 'Eye ID'
            header.eye = field_val;
        case '(x,y,z)'
            dims = strrep(field_val,'(','');            
            dims = strrep(dims,')','');            
            switch category
                case 'OCT image'
                    dims = strsplit(dims,',');
                    
                    header.n_ascan = str2num(dims{1});
                    header.n_bscan = str2num(dims{2});
                    header.n_axial = str2num(dims{3});                        
                case 'Scan parameters'
                    dims = strrep(dims,'mm','');
                    dims = strsplit(dims,',');

                    header.size_x = str2num(dims{1});
                    header.size_y = str2num(dims{2});
                    header.size_z = str2num(dims{3});      
                case 'Color fundus photo'
                    dims = strsplit(dims,',');

                    dims = cellfun(@str2num, dims);
                    header.fundus_dims = dims([2 1 3]);
                otherwise
                    error('Unknown category for (x,y,z)');
            end            
        case '(x,y)'
            dims = strrep(field_val,'(','');            
            dims = strrep(dims,')','');       
            dims = strsplit(dims,',');
            dims = cellfun(@str2num, dims);
            header.oct_fundus_ul = dims(1:2);
            header.oct_fundus_br = dims(3:4);
        case 'Image Quality'        
            header.quality = str2double(field_val);
        case 'Fixation Type'                        
            header.fixation_type = field_val;                                                                
        case 'Scan Protocol'        
            header.scan_protocol = field_val;                    
        case 'Scan mode'
            header.scan_mode = field_val;
        case 'Scan date'
            header.scan_date = field_val;
        case 'Scan time (24 HR)'
            header.scan_time = field_val;    
        case 'Version'
            header.seg_version = str2num(field_val);
        case 'Layers'
            header.n_layer = str2num(field_val);
        case '(Ascan,Frame)'
            vals = strrep(field_val,'(','');            
            vals = strrep(vals,')','');      
            vals = strsplit(vals, ',');
            vals = cellfun(@str2num, vals);
            
            switch category
                case 'Fovea Center'
                    header.fovea_center = vals;
                case 'Disc Center'
                    header.disc_center = vals;
                otherwise
                    error('Unknown category for (Ascan,Frame)');
            end
        otherwise
            field_name = strrep(field_name,'(',' ');
            field_name = strrep(field_name,')',' ');
            field_name = strrep(field_name,',','_');        
            field_name = strtrim(field_name);
            field_name = strrep(field_name,' ','_');        
            field_name = lower(field_name);

            header.(field_name) = field_val;
    end         
end
   
header.scale_x = header.size_x/(header.n_ascan-1);
header.scale_y = header.size_y/(header.n_bscan-1);
header.scale_z = header.size_z/(header.n_axial-1);

if coordinates
   
    x = linspace(-header.size_x/2, header.size_x/2, header.n_ascan);
    y = linspace(header.size_y/2, -header.size_y/2, header.n_bscan);

    [header.X_oct, header.Y_oct] = meshgrid(x,y);
end
    
function header = read_metadata2(file, header, coordinates)
if ~exist(file,'file')
    error(['Unable to find segmentation metadata file: ' file]);
end

fid = fopen(file);
meta_data = char(fread(fid)');

lines = splitlines(meta_data);
is_empty = cellfun(@isempty, lines);
lines(is_empty) = [];

n_line = length(lines);

for i_line=1:n_line
    line = lines{i_line};
    
    idx = strfind(line,':');
    if isempty(idx)
        idx = strfind(line,'=');
    end
    
    field_name = strtrim(line(1:idx(1)-1));
    field_val = strtrim(line(idx(1)+1:end));
    
    switch field_name
        case 'Data'
            header.file_id = field_val;
        case '(Macula_Center-frame, Macula_Center-aline)'
            field_val = strrep(field_val,')','');
            field_val = strrep(field_val,'(','');
            field_val = strsplit(field_val,',');
            
            header.fovea_center_px = cellfun(@str2num, field_val);
        case '(Disc_Center-frame, Disc_Center-aline)'
            field_val = strrep(field_val,')','');
            field_val = strrep(field_val,'(','');
            field_val = strsplit(field_val,',');
            
            header.onh_center_px = cellfun(@str2num, field_val);
        case 'ONH size'
            header.onh_size = str2num(field_val);
        case 'ILM indicator'
            header.ilm_indicator = str2num(field_val);
        case 'Valid degree'
            header.valid_degree = str2num(field_val);
        case 'Min motion correlation'
            header.min_motion_correlation = str2num(field_val);
        case 'Max motion delta'
            header.max_motion_delta = str2num(field_val);
        case 'Max motion factor'
            header.max_motion_factor = str2num(field_val);
        otherwise
    end
end

% Foveal center in mm
if coordinates
    header.x_fovea = header.X_oct(1, header.fovea_center_px(2));
    header.y_fovea = header.Y_oct(header.fovea_center_px(1),1);
end

fclose(fid);

function header = read_seg_layers(file, header)
if ~exist(file,'file')
    warning(['Unable to find segmentation .txt file: ',file,...
        ' layer names might be inaccurate']);
    header.layers = {'ILM','RNFL_GCL','GCL_IPL','IPL_INL','INL_OPL','ELM','ISOS','RPE','BM','CSI'};
    return;
end

fid = fopen(file);

data = char(fread(fid)');
lines = splitlines(data);
layers = strsplit(lines{1}, '\t');

% Remove first column
header.layers = layers(2:end);

fclose(fid);

function seg = read_seg(file, header)

if ~exist(file,'file')
    error(['Unable to find segmentation .dat file: ' file]);
end

fid = fopen(file);
seg_data = fread(fid,'*uint16');

dims = [header.n_ascan header.n_bscan header.n_layer];

seg_data = reshape(seg_data, dims);
seg_data = permute(seg_data, [2 1 3]);

seg = struct;
for i_layer=1:header.n_layer
    seg.(header.layers{i_layer}) = seg_data(:,:,i_layer);
end

function bscan = read_bscan(file, header)
if ~exist(file,'file')
    error(['Unable to find bscan data file: ' file]);
end

fid = fopen(file);

I = fread(fid, '*uint16');

% Remove useless data (not necessary if we read it as uint16)
% I = I(2:2:end); % other is just 0s 

dims = [header.n_ascan header.n_axial header.n_bscan];
% dims = [512 650 128]; % hard coded

bscan = reshape(I, dims);
bscan = permute(bscan, [3 2 1]);

fclose(fid);

function fundus = read_fundus(file, header)
if ~exist(file,'file')
    error(['Unable to find colour fundus file: ' file]);
end

fid = fopen(file);
fundus = fread(fid, '*uint8');

% fundus_dims = [3 2048 1536];
% Define colour fundus image dimensions
dims = flip(header.fundus_dims);

fundus = reshape(fundus, dims);
fundus = permute(fundus, [3 2 1]);
fundus = flip(fundus,3);

fclose(fid);
