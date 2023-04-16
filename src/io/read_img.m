function [header, bscan] = read_img(file, scan_size, varargin)
% Read OCT images from Cirrus .img files
%
%
% Input arguments
% --------------- 
% * **file**:        String with the path to the img file name to be read.
%           
% * **scan_size**:   Array of 3 numeric values defining the size of the 
%                    scanned region. Sizes must be given following [x y z]
%                    convention. By default a cube of [6 6 2] mm is 
%                    considered.
%
% * **varargin**:    optional string flags from the list:
%                       
%   - 'get_coordinates': retrieve fundus and A-Scan X, Y coordinates. If the scan pattern is unknown coordinates cannot be computed.
%  
%
% Output arguments
% ---------------- 
% * **header**: struct with metadata retrieved from filename and volume dimensions.
%
% * **bscan**: 3D matrix of size [n_axial x n_ascan x n_bscan] with image data.
%
%
% Notes
% -----
% img files do not include any spatial or metadata information. Therefore,
% to correctly reconstruct the image volume this function relies on two
% strategies:
%
% 1. Inspection of the filename.
% 2. Matching of the voxel count with this known acquisition patterns:
%
%   - Macular Cube: [200 x 200 x 1024] : 40960000
%   - Opitic Disc Cube: [512 x 128 x 1024] : 67108864
%   - 5 Line Raster: [1024 x 5 x 1024]
%   - Hidef: [512 x 2 x 1024] (not supported yet)
%
% Scans differing from the previous protocols might not be read properly. 
%
% Info from previous function:
% hidef scan consists of 2 orthogonal bscans intersecting in the center of
% the volume. The first is in the y direction (across B-scans), the second
% is in the x direction (aligned with the center B-scan). Unfortunately,
% they are not in good alignment with the volumetric data.
%
%
% Example 
% -------      
% .. code-block:: matlab
%
%   [header, bscan] = read_img('myfile.img')    

if nargin < 2 || isempty(scan_size)
    scan_size = [6 6 2]; % assume 6mm x 6mm x 2mm cube size
end

get_coordinates = any(strcmp('get_coordinates', varargin));

% Parse file name for useful information
[~, fname, ext] = fileparts(file);
if ~strcmp(ext,'.img')
    error(['Input file must be of type .img and not ' ext '.'])
end

header = data_from_filename(fname);               

% Read image data 
fid = fopen(file);
if fid == -1
    error('Could not open the file. Check the path is correct.');
end
A = fread(fid,'*uchar');
fclose(fid);

% Guess dimensions from voxel count
n_voxel = numel(A);
[dims_count_ok, dims] = guess_dimensions(n_voxel);

% Check if dimensions from filename are ok and make sense
if isempty(header.n_ascan)
    dims_header_ok = false;
else
    header.n_axial = n_voxel/(header.n_ascan * header.n_bscan);
    dims_header_ok = rem(header.n_axial,1) ==0;        
end

if ~dims_header_ok & ~dims_count_ok
    warning("Unable to compute image dimensions. Imposible to reconstruct the volume");
    bscan = A;
    return;
elseif ~dims_header_ok & dims_count_ok
    header.protocol = dims.protocol;
    header.n_ascan  = dims.n_ascan;
    header.n_bscan  = dims.n_bscan;
    header.n_axial  = dims.n_axial;
elseif dims_header_ok & dims_count_ok % Double check
    vals_fname = [header.n_ascan header.n_bscan header.n_axial];
    vals_vox = [dims.n_ascan dims.n_bscan dims.n_axial];
    if ~isequal(vals_fname, vals_vox)
        warning("Dimensions from filename and vol size do not match. Keeping dimensions from filename");
    end
end
% When header dims are ok and count dims not we do nothing.

% Form the volume
bscan = reshape(A, header.n_ascan, header.n_axial, header.n_bscan);
bscan = permute(bscan, [2 1 3]); % correct order: [z, x, y]

% Reshape - bscans are rotated
bscan = flip(bscan, 1); % upside down (z axis (depth))
bscan = flip(bscan, 2); % flip left and right (x axis)

% Resolution in each direction - assume cube begins and ends at voxel
% center
header.scale_x = scan_size(1)/(header.n_ascan - 1);
header.scale_y = scan_size(2)/(header.n_bscan - 1);
header.scale_z = scan_size(3)/(header.n_axial - 1);

if get_coordinates
    [header.X_oct, header.Y_oct] = get_ascan_coordinates(header);
end

function header = data_from_filename(fname)

header.pid = [];
header.scan_date = [];
header.n_ascan = [];
header.n_bscan = [];
header.eye = [];

% Divide filename into chunks based on underscores
% if ~isempty(strfind(name,'clinic'))
if contains(fname, 'clinic')
    C = textscan(fname, '%s %s %s %s %s %s %s %s %s %s', 'Delimiter', '_');
    C{1}{1} = [C{1}{1} '_' C{2}{1}];
    C(2) = [];
else
    C = textscan(fname, '%s %s %s %s %s %s %s %s %s', 'Delimiter', '_');
end

% 1. Patient ID
header.pid = C{1}{1}; % patient ID

% 2. Scan type and volume dimensions
if isempty(C{2}{1})
    warning("Could not obtain info from file name");
    return;
end
scan_data = C{2}{1};
C2 = textscan(scan_data, '%s %s %s %s');

if contains(scan_data, 'Macular Cube')
    header.protocol      = 'Macular Cube';
    header.fixation      = 'macula';
    header.bscan_pattern = 'raster';
    
    ind = strfind(C2{3}{1}, 'x');
    if length(ind)==1 & ind >1 & ind <length(C2{3}{1})
        header.n_ascan = str2double(C2{3}{1}(1:ind-1));
        header.n_bscan = str2double(C2{3}{1}(ind+1:end));
    end

elseif contains(scan_data, 'Optic Disc Cube')
    header.protocol      = 'Optic Disc Cube';
    header.fixation      = 'onh';
    header.bscan_pattern = 'raster';

    ind = strfind(C2{4}{1}, 'x');
    if length(ind)==1 & ind >1 & ind <length(C2{4}{1})
        header.n_ascan = str2double(C2{4}{1}(1:ind-1));
        header.n_bscan = str2double(C2{4}{1}(ind+1:end));
    end
    
elseif contains(scan_data, 'HD 5 Line Raster')
    header.protocol      = 'HD 5 Line Raster';   
    header.fixation      = 'unknown';  % macula probably
    header.bscan_pattern = 'raster';
    header.n_ascan = 1024;
    header.n_bscan = 5;

elseif contains(scan_data, 'HD 21 Line')
    header.protocol      = 'HD 21 Line';    
    header.fixation      = 'unkwnon';  % macula probably
    header.bscan_pattern = 'raster';
    header.n_ascan = 1024;
    header.n_bscan = 21;    

elseif contains(scan_data, 'Anterior Segment Cube')
    header.protocol      = 'Anterior Segment Cube';
    header.fixation      = 'anterior';
    header.bscan_pattern = 'raster';
    
    ind = strfind(C2{4}{1}, 'x');
    if length(ind)==1 & ind >1 & ind <length(C2{4}{1})
        header.n_ascan = str2double(C2{4}{1}(1:ind-1));
        header.n_bscan = str2double(C2{4}{1}(ind+1:end));
    end
else
    header.protocol      = 'unknown';
    header.fixation      = 'unknown';
    header.bscan_pattern = 'unknown';
    
    warning("Could not retrieve scan type from file name");
end

% 3./4. Scan-date
if isempty(C{3}) | isempty(C{4})
    warning("Could not retrieve full data from filename");
    return;
end
header.scan_date = [C{3}{1} ' ' C{4}{1}];

% 5. Eye (OD/OS)
if isempty(C{5})
    warning("Could not retrieve full data from filename");
    return;
end
if ~any(strcmp(C{5}{1}, {'OD','OS'}))
    warning("Retrieved eye from filename is neither OD nor OS. Caution");
end
header.eye = C{5}{1};

% 6. Scan number
header.scan_number = C{6}{1}; % 'snXXXX'

% 7./8. Extra info about file (useful only to discern cube_raw, cube_z, hidef)
n_field = sum(~cellfun(@isempty, C));
if n_field == 8
    header.last = [C{7}{1} '_' C{8}{1}]; % 'cube'
elseif n_field == 7
    header.last = C{7}{1}; % 'z' or 'raw'    
end

function  [dims_ok, dims] = guess_dimensions(n_voxel)
% Warning: anterior segment cube has the same number of voxels as macular_cube
% and is therefore indistinguishable based solely on voxel count.
dims_ok = true;
switch n_voxel
    case 67108864
        dims.protocol      = 'Macular Cube';
        dims.fixation      = 'macula';
        dims.bscan_pattern = 'raster';
        dims.n_ascan = 512;
        dims.n_bscan = 128;
        dims.n_axial = 1024;
    case 40960000
        dims.protocol     = 'Optic Disc Cube';
        dims.fixation      = 'onh';
        dims.bscan_pattern = 'raster';
        dims.n_ascan= 200;
        dims.n_bscan = 200;
        dims.n_axial = 1024;
    case 22020096
        dims.protocol      = 'HD 21 Line';
        dims.fixation      = 'unknown';  % macula probably
        dims.bscan_pattern = 'raster';
        dims.n_ascan = 1024;
        dims.n_bscan = 21;
        dims.n_axial = 1024;        
    case 20971520
        dims.protocol      = 'HD 5 Line Raster (wide)';
        dims.fixation      = 'unknown';  % macula probably
        dims.bscan_pattern = 'raster';
        dims.n_ascan = 4096;
        dims.n_bscan = 5;
        dims.n_axial = 1024;
    case 5242880
        dims.protocol      = 'HD 5 Line Raster';
        dims.fixation      = 'unknown';  % macula probably
        dims.bscan_pattern = 'raster';
        
        dims.n_ascan = 1024;
        dims.n_bscan = 5;
        dims.n_axial = 1024;
    otherwise
        dims_ok = false;
        dims = nan;
end
