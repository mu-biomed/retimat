function [header, bscan] = read_img(file, scan_size, get_coord)
%READ_IMG Read data from .img files from Zeiss Cirrus scans
%
%   [header, bscan] = read_img(file, scan_size, get_coord)
%   Reads the image volume stored in 'file' and tries to retrieve spatial 
%   information.
%
%   Input arguments:
%  
%   'file'           String with the path to the .img file name to be read.
%           
%   'scan_size'      Array of 3 numeric values defining the size of the 
%                    scanned region. Sizes must be given following [x y z]
%                    convention. By default a cube of [6 6 2] mm is 
%                    considered.
%
%   'get_coord'      If true it tries to compute a-scan coordinates.
%                    Default is false. If the scan patter is unknown 
%                    coordinates cannot be computed.
%  
%
%   Output arguments:
%  
%   'header'         Basic header information retrieved from filename and
%                    volume dimensions.
%
%   'bscan'          3D matrix of size [n_axial x n_ascan x n_bscan] with 
%                    image data.
%
%
%
%   Notes
%   -----
%   .img files do not include any spatial or metadata information. Therefore,
%   to correctly reconstruct the image volume this function relies on two
%   strategies:
%   1. Inspection of the filename, which has a 
%   2. Matching of the voxel count with this known acquisition patterns:
%     - Macular Cube: [200 x 200 x 1024] : 40960000
%     - Opitic Disc Cube: [512 x 128 x 1024] : 67108864
%     - 5 Line Raster: [1024 x 5 x 1024]
%     - Hidef: [512 x 2 x 1024] (not supported yet)
%
%   Scans differing from the previous might not be read properly. 
%
%  Info from previous function:
%  hidef scan consists of 2 orthogonal bscans intersecting in the center of
%  the volume. The first is in the y direction (across B-scans), the second
%  is in the x direction (aligned with the center B-scan). Unfortunately,
%  they are not in good alignment with the volumetric data.
%
%   References
%   ----------
%   [1] 
%
%   Example 
%   ---------      
%   % Read img file
%
%     [header, bscan] = read_img('myfile.img')    
%
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2021

if nargin < 2 || isempty(scan_size)
    scan_size = [6 6 2]; % assume 6mm x 6mm x 2mm cube size
end

if nargin < 3 || isempty(get_coord)
    get_coord = false;
end

% Parse file name for useful information
[~, fname, ext] = fileparts(file);
if ~strcmp(ext,'.img')
    error(['Input file must be of type .img and not ' ext '.'])
end

header = data_from_filename(fname);               

% Get volume dimensions
% if strcmp(cube,'hidef')
%     n_bscan = 2;
%     if strcmp(header.scan_type,'Macular')
%         n_ascan = 1024;
%     elseif strcmp(header.scan_type,'Optic Disc')
%         n_ascan = 1000;
%     end
% end


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
    header.scan_type = dims.scan_type;
    header.n_ascan = dims.n_ascan;
    header.n_bscan = dims.n_bscan;
    header.n_axial = dims.n_axial;
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
bscan = flip(bscan, 3); % flip front to back (y axis)

% Resolution in each direction - assume cube begins and ends at voxel
% center
header.scale_x = scan_size(1)/(header.n_ascan - 1);
header.scale_y = scan_size(2)/(header.n_bscan - 1);
header.scale_z = scan_size(3)/(header.n_axial - 1);

if get_coord
    x_max = scan_size(1)./2;
    y_max = scan_size(2)./2;
        
    % Pending to check carefully up-down or down-up
    y_range = linspace(y_max, -y_max, header.n_bscan);
    
    % Apparently
    % OD: already pointing nasal
    % OS: need to flip
    if strcmp(header.eye,'OD')
        x_range = linspace(-x_max, x_max, header.n_ascan);
    else
        x_range = linspace(x_max, -x_max, header.n_ascan);
    end
    [header.X, header.Y] = meshgrid(x_range, y_range);    
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
    warning("Could not obtain info from filename");
    return;
end
scan_data = C{2}{1};
C2 = textscan(scan_data, '%s %s %s %s');

if contains(scan_data, 'Macular Cube')
    header.scan_type = 'macular_cube';

    ind = strfind(C2{3}{1}, 'x');
    if length(ind)==1 & ind >1 & ind <length(C2{3}{1})
        header.n_ascan = str2double(C2{3}{1}(1:ind-1));
        header.n_bscan = str2double(C2{3}{1}(ind+1:end));
    end

elseif contains(scan_data, 'Optic Disc Cube')
    header.scan_type = 'optic_disc_cube';
    
    ind = strfind(C2{4}{1}, 'x');
    if length(ind)==1 & ind >1 & ind <length(C2{4}{1})
        header.n_ascan = str2double(C2{4}{1}(1:ind-1));
        header.n_bscan = str2double(C2{4}{1}(ind+1:end));
    end
    
elseif contains(scan_data, 'HD 5 Line Raster')
    header.scan_type = '5line_raster';    
    header.n_ascan = 1024;
    header.n_bscan = 5;

elseif contains(scan_data, 'HD 21 Line')
    header.scan_type = '21line_raster_wide';    
    header.n_ascan = 1024;
    header.n_bscan = 21;    

elseif contains(scan_data, 'Anterior Segment Cube')
    header.scan_type = 'anterior_segment_cube';
   
    ind = strfind(C2{4}{1}, 'x');
    if length(ind)==1 & ind >1 & ind <length(C2{4}{1})
        header.n_ascan = str2double(C2{4}{1}(1:ind-1));
        header.n_bscan = str2double(C2{4}{1}(ind+1:end));
    end
else
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
        dims.scan_type = 'macular_cube';
        dims.n_ascan = 512;
        dims.n_bscan = 128;
        dims.n_axial = 1024;
    case 40960000
        dims.scan_type = 'optic_disc_cube';
        dims.n_ascan= 200;
        dims.n_bscan = 200;
        dims.n_axial = 1024;
    case 22020096
        dims.scan_type = '21line_raster_wide';
        dims.n_ascan = 1024;
        dims.n_bscan = 21;
        dims.n_axial = 1024;        
    case 20971520
        dims.scan_type = '5line_raster_wide';
        dims.n_ascan = 4096;
        dims.n_bscan = 5;
        dims.n_axial = 1024;
    case 5242880
        dims.scan_type = '5line_raster';
        dims.n_ascan = 1024;
        dims.n_bscan = 5;
        dims.n_axial = 1024;
    otherwise
        dims_ok = false;
        dims = nan;
end
