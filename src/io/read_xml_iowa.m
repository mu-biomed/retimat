function [header, seg] = read_xml_iowa(file, varargin)
% Read segmentation computed by OCTExplorer (IOWA reference algorithm)
%
% Input arguments
% ---------------
% * **file**:    path to xml output file ending with Surfaces_Iowa.xml
%            
% * **varargin**: different string flags:
%
%   - 'get_coordinates': to return A-Scan coordinates as part of the header.
%   - 'keep_names': to use layer naming convention from OCTExplorer instead of RETIMAT.
%   - 'keep_nan': to keep regions with noisy segmentation as numerical values instead of nan.
%   - 'verbose': to plot additional info on read layers.
%  
% Output arguments
% ---------------- 
% * **header**: struct with metadata related with the segmentation: version, resolution, dimensions.
%
% * **seg**: struct with point segmentation values (in pixels).          
%  
% Notes
% -----
% This function has been only tested for raster acquisitions and
% OCTExplorer version version 3.8.0.
%
% The axes convention used in OCTExplorer differs from the one used here:
%
% * x: horizontal (temporal - nasal)
% * y: vertical (inferior - superior)
% * z: axial (depth)
%
% Layer naming convention returned by this function differs from the one
% used by OCTExplorer.
%
% OCTExplorer uses different ways to determine an invalid region:
%
% * Segmentation = 0
% * Segmentation NA (not sure)
% * Undefined region chunk (full square regions I believe)
%   
%
% References
% ----------
% [1] The Iowa Reference Algorithms (Retinal Image Analysis Lab, Iowa 
% Institute for Biomedical Imaging, Iowa City, IA), 
% https://www.iibi.uiowa.edu/oct-reference
%
% [2] Abramoff MD et al., Retinal Imaging and Image Analysis. IEEE 
% Reviews in Biomedical Engineering, 2010. doi:10.1109/RBME.2010.2084567
%   
% [3] Kang L et al., Optimal Surface Segmentation in Volumetric Images – 
% A Graph-Theoretic Approach. IEEE Transactions on Pattern Analysis and 
% Machine Intelligence, 2006. doi:10.1109/TPAMI.2006.19
%   
% [4] Garvin MK et al., Automated 3-D Intraretinal Layer Segmentation of
% Macular Spectral-Domain Optical Coherence Tomography Images. IEEE Trans 
% Med. Imaging, 2009. doi:10.1109/TMI.2009.2016958
%
%
% Example
% -------      
% .. code-block:: matlab
%
%   [header, seg] = read_xml_iowa('my_file_Surfaces_Iowa.xml')

get_coordinates = any(strcmp(varargin, 'get_coordinates'));
keep_names      = any(strcmp(varargin, 'keep_names'));
keep_nan        = any(strcmp(varargin, 'keep_nan'));
verbose         = any(strcmp(varargin, 'verbose'));

% Read whole xml as text
str = fileread(file);

% Get OCTExplorer and executable version, manufacturer and laterality (eye)
[idx1, idx2] = get_index(str, 'version');
header.version      = str(idx1(1)+9:idx2(1)-1);
header.exec_version = str(idx1(2)+9:idx2(2)-1);

header.manufacturer  = get_field(str, 'manufacturer', 'str');
header.eye           = get_field(str, 'laterality',   'str');
header.fixation      = lower(get_field(str, 'center_type',  'str'));
header.bscan_pattern = 'raster';

% Dimensions
[idx1, idx2] = get_index(str, 'size');
header.n_ascan = get_field(str(idx1:idx2), 'x', 'num');
header.n_axial = get_field(str(idx1:idx2), 'y', 'num');
header.n_bscan = get_field(str(idx1:idx2), 'z', 'num');

% Resolution
[idx1, idx2] = get_index(str, 'voxel_size');
header.scale_x = get_field(str(idx1:idx2), 'x', 'num');
header.scale_z = get_field(str(idx1:idx2), 'y', 'num');
header.scale_y = get_field(str(idx1:idx2), 'z', 'num');

% Surfaces (segmentation)
[idx1, idx2] = get_index(str, 'surface');
n_seg = length(idx1);

seg = struct;
for i_seg=1:n_seg
    % Extract surface data
    chunk_seg = str(idx1(i_seg):idx2(i_seg));
    
    % Surface name
    seg_name = get_field(chunk_seg, 'name', 'str');    
    if ~keep_names
        seg_name = change_name(seg_name);
    else
        seg_name = strrep(seg_name, '(', '_'); % to avoid struct errors
        seg_name = strrep(seg_name, ')', '_'); 
        seg_name = strrep(seg_name, ' ', '_'); 
        seg_name = strrep(seg_name, '-', '_');
        seg_name = strrep(seg_name, '/', '_');        
    end
    if verbose
        disp(['Reading layer: ' seg_name]);
    end
    
    % Get B-scans
    [idx_b1, idx_b2] = get_index(chunk_seg, 'bscan');
    n_bscan = length(idx_b1);
    
    seg.(seg_name) = nan(header.n_bscan, header.n_ascan);
    for i_bscan=1:n_bscan
        chunk_bscan = chunk_seg(idx_b1(i_bscan):idx_b2(i_bscan));
        
        % Find A-scan
        [idx_a1, idx_a2] = get_index(chunk_bscan, 'y');
        n_ascan = length(idx_a1);
        
        % Fast method
        n_len = idx_a2 - idx_a1 - 3;
        n_len_unique = unique(n_len);
        
        if length(n_len_unique) == 1
            vals = str2num(chunk_bscan(idx_a1'+repmat(3:2+n_len_unique,n_ascan,1)));            
        else        
            vals = cellstr(chunk_bscan(idx_a1'+repmat(3:2+max(n_len),n_ascan,1)));
            vals = cellfun(@(x,y) str2double(x(1:n_len(y))), vals, num2cell(1:n_ascan)');
        end
        
        seg.(seg_name)(i_bscan,:) = vals;
        
        % Slow method
%         for i_ascan=1:n_ascan
%             seg.(seg_name)(i_bscan,i_ascan) = str2double(chunk_bscan(idx_a1(i_ascan)+3:idx_a2(i_ascan)-1));
%         end
    end
end

% Set to nan rois labeled as 'undefined region' due to poor quality
% This is slow and most likely improvable
if ~keep_nan    
    [idx1, idx2] = get_index(str, 'undefined_region');
    chunk = str(idx1:idx2);
    [idx_a1, idx_a2] = get_index(chunk, 'x');
    [idx_b1, idx_b2] = get_index(chunk, 'z');

    n_nan = length(idx_a1);
    is_nan = false(n_bscan, n_ascan);
    for i_nan=1:n_nan
        ascan_nan = str2double(chunk((idx_a1(i_nan)+3) : (idx_a2(i_nan)-1))) + 1;
        bscan_nan = str2double(chunk((idx_b1(i_nan)+3) : (idx_b2(i_nan)-1))) + 1;
        is_nan(bscan_nan, ascan_nan) = true;
    end
    
    layers = fields(seg);
    for i_layer=1:n_seg
        seg.(layers{i_layer})(is_nan) = nan;
    end
end

if get_coordinates
    [header.X_oct, header.Y_oct] = get_ascan_coordinates(header);
end

function [idx1, idx2] = get_index(text, label)

idx1 = strfind(text, ['<' label '>']);
idx2 = strfind(text, ['</' label '>']);

function seg_name = change_name(seg_name)

switch seg_name
    case 'ILM (ILM)'
        seg_name = 'ILM';
    case 'RNFL-GCL (RNFL-GCL)'
        seg_name = 'RNFL_GCL';
    case 'GCL-IPL (GCL-IPL)'
        seg_name = 'GCL_IPL';
    case 'IPL-INL (IPL-INL)'
        seg_name = 'IPL_INL';
    case 'INL-OPL (INL-OPL)'
        seg_name = 'INL_OPL';
    case 'OPL-Henles fiber layer (OPL-HFL)'
        seg_name = 'OPL_ONL';
    case 'Boundary of myoid and ellipsoid of inner segments (BMEIS)'
        seg_name = 'MZ_EZ';
    case 'IS/OS junction (IS/OSJ)'
        seg_name = 'EZ_OSP';
    case 'Inner boundary of OPR (IB_OPR)'
        seg_name = 'OSP_IZ';
    case 'Inner boundary of RPE (IB_RPE)'
        seg_name = 'IZ_RPE';
    case 'Outer boundary of RPE (OB_RPE)'
        seg_name = 'BM';
    otherwise
        warning(['Unrecognized layer: ' seg_name]);
end
    
function val = get_field(text, label, val_type)

n_char = length(label);

[idx1, idx2] = get_index(text, label);

if length(idx1) >1 || length(idx2) >1
    error(["Field " label " is present more than once in text"]);
end

if isempty(idx1) || isempty(idx2)
    error(["Field " label " not found entirely in text"]);
end

switch val_type
    case 'str'
        val = text(idx1+n_char+2:idx2-1);
    case 'num'
        val = str2double(text(idx1+n_char+2:idx2-1));
    otherwise
        error("Erroneous data type");
end
