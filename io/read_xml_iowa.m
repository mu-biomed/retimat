function [header, seg] = read_xml_iowa(file)

% Read whole xml as text
str = fileread(file);

% Get OCTExplorer and executable version
[idx1, idx2] = get_index(str, 'version');
header.version = str(idx1(1)+9:idx2(1)-1);
header.exec_version = str(idx1(2)+9:idx2(2)-1);

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

% Surfaces
[idx1, idx2] = get_index(str, 'surface');
n_seg = length(idx1);

seg = struct;
for i_seg=1:n_seg
    % Extract surface data
    chunk_seg = str(idx1(i_seg):idx2(i_seg));
    
    % Surface name
    seg_name = get_field(chunk_seg, 'name', 'str');    
    seg_name = change_name(seg_name);
    
    disp(['Reading segmentation:' seg_name]);
    % Get B-scans
    [idx_b1, idx_b2] = get_index(chunk_seg, 'bscan');
    n_bscan = length(idx_b1);
    
    seg.(seg_name) = nan(header.n_bscan, header.n_ascan);
    for i_bscan=1:n_bscan
        chunk_bscan = chunk_seg(idx_b1(i_bscan):idx_b2(i_bscan));
        
        % Find A-scan
        [idx_a1, idx_a2] = get_index(chunk_bscan, 'y');
        n_ascan = length(idx_a1);
        
        for i_ascan=1:n_ascan
            seg.(seg_name)(i_bscan,i_ascan) = str2double(chunk_bscan(idx_a1(i_ascan)+3:idx_a2(i_ascan)-1));
        end
    end
end

disp('finished');

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
        seg_name = 'OPL_HFL';
    case 'Boundary of myoid and ellipsoid of inner segments (BMEIS)'
        seg_name = 'BMEIS';
    case 'IS/OS junction (IS/OSJ)'
        seg_name = 'ISOSJ';
    case 'Inner boundary of OPR (IB_OPR)'
        seg_name = 'IB_OPR';
    case 'Inner boundary of RPE (IB_RPE)'
        seg_name = 'IB_RPE';
    case 'Outer boundary of RPE (OB_RPE)'
        seg_name = 'OB_RPE';
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
