function [header, segment, bscan, slo] = read_vol(file, options)

% Read Heidelberg Engineering (HE) OCT raw files (VOL ending)
% [header, segment, bscan, slo] = read_vol(file, options)
% This function performs volume OCT data (xxxx.vol) reading. 
% HEADER: Header information as described by HE. Struct with each entry
%   named by the HE conventions. 
% BSCANHEADER: B-scan header information. Struct with each entry named by 
%   the HE conventions. The entries consist of vectors/matrices, with one 
%   field per B-Scan. 
% SLO: Slo image as unsigned integers.
% BScans: BScans. 3D matrix with floating point data of the B-Scans.
% file: Filename of the VOL-file to read, with ending.
%
% OPTIONS: Various output possibilites, 
%   written in the options variable as string text, i.e. 'visu writemeta'
%   Possibilities: 
%       visu: The read-in is visualized.
%       visuseg: The HE segmentation is also visualized.
%       header: Only the Header and BScansHeader Information is read, 
%            not the image data  
%       nodisp: nothing is diplayed during read in
%
% Originally writen by Radim Kolar, Brno University, Czech Republic
% Modified by Markus Mayer, Pattern Recognition Lab, University of
% Erlangen-Nuremberg
% Modified by Kris Sheets, Retinal Cell Biology Lab, 
% Neuroscience Center of Excellence, LSU Health Sciences Center, 
% New Orleans
% 
% Modified by Andrew Lang, Image Analysis and Communications Lab, Johns
% Hopkins University - Modifications to increase efficiency - 5/17/2012 
% 
% Modified by David Romero, Biomedical Engineering Department, Mondragon
% Unibertsitatea, 11/08/2021
%
% You may use this code as you want. I would be grateful if you would go to
% our homepage look for articles that you find worth citing in your next
% publication:
% http://www.vutbr.cz/lide/radim-kolar-2796/publikace
% http://www5.informatik.uni-erlangen.de/en/our-team/mayer-markus
% Thanks, Radim & Markus


% Configure visualization options
if nargin==1
    options = '';
end

visu = 0;
visuseg = 0;

if numel(strfind(options, 'visu')) ~= 0
    visu = 1;
end
if numel(strfind(options, 'visuseg')) ~= 0
    visuseg = 1;
end

% Open the file
fid = fopen(file);
 
% Read header
version = fread(fid, 12, '*int8');
n_ascan = fread(fid, 1, '*int32');
n_bscan = fread(fid, 1, '*int32');
n_axial = fread(fid, 1, '*int32');
scale_x = fread(fid, 1, '*double');
scale_y = fread(fid, 1, '*double');
scale_z = fread(fid, 1, '*double');
size_x_slo = fread(fid, 1, '*int32');
size_y_slo = fread(fid, 1, '*int32');
scale_x_slo = fread(fid, 1, '*double');
scale_y_slo = fread(fid, 1, '*double');
fov_slo = fread(fid, 1, '*int32');
scan_focus = fread(fid, 1, '*double');
eye = char(fread(fid, 4, '*uchar')');
exam_time = fread(fid, 1, '*int64');
scan_pattern = fread(fid, 1, '*int32');
bscan_hdr_size = fread(fid, 1, '*int32');
id = char(fread(fid, 16, '*uchar')');
reference_id = char(fread(fid, 16, '*uchar')');
pid = fread(fid, 1, '*int32');
patient_id = char(fread(fid, 21, '*uchar')');
padding = fread(fid, 3, '*int8');
dob = fread(fid, 1, '*double');
vid = fread(fid, 1, '*int32');
visit_id = char(fread(fid, 24, '*uchar')');
visit_date = fread(fid, 1, '*double');
grid_type = fread(fid, 1, '*int32');
grid_offset = fread(fid, 1, '*int32');
spare = fread(fid, 1832, '*int8');
 
if numel(strfind(options, 'nodisp')) == 0
    disp(['---------------------------------------------']);
    disp(['           Version: ' char(version')]);
    disp(['         NumAScans: ' num2str(n_ascan)]);
    disp(['         NumBScans: ' num2str(n_bscan)]);
    disp(['             SizeZ: ' num2str(n_axial)]);
    disp(['            ScaleX: ' num2str(scale_x) ' mm']);
    disp(['            ScaleY: ' num2str(scale_y) ' mm']);
    disp(['            ScaleZ: ' num2str(scale_z) ' mm']);
    disp(['          SizeXSlo: ' num2str(size_x_slo)]);
    disp(['          SizeYSlo: ' num2str(size_y_slo)]);
    disp(['         ScaleXSlo: ' num2str(scale_x_slo) ' mm']);
    disp(['         ScaleYSlo: ' num2str(scale_y_slo) ' mm']);
    disp(['FieldSizeSlo (FOV): ' num2str(fov_slo) 'deg']);
    disp(['         ScanFocus: ' num2str(scan_focus) ' dpt']);
    disp(['               Eye: ' char(eye)]);
    disp(['          ExamTime: ' datestr(double(exam_time(1)/(1e7*60*60*24)+584755+(2/24)))]);
    disp(['       ScanPattern: ' num2str(scan_pattern)]);
    disp(['      BScanHdrSize: ' num2str(bscan_hdr_size) ' bytes']);
    disp(['                ID: ' char(id)]);
    disp(['       ReferenceID: ' char(reference_id)]);
    disp(['               PID: ' num2str(pid)]);
    disp(['         PatientID: ' char(patient_id)]);
    disp(['               DOB: ' datestr(dob+693960)]);
    disp(['               VID: ' num2str(vid)]);
    disp(['           VisitID: ' char(visit_id)]);
    disp(['         VisitDate: ' datestr(double(visit_date+693960))]);
    disp(['          GridType: ' num2str(grid_type)]);
    disp(['        GridOffset: ' num2str(grid_offset)]);
    disp(['---------------------------------------------']);
end

if nargout == 1
    fclose(fid);
    return
end
status = fseek(fid, 2048, -1 );


% Read fundus image (slo)
if(~any(strcmp(options, 'header')))
    slo = fread(fid, size_x_slo*size_y_slo, '*uint8');
    slo = reshape(slo, size_x_slo, size_y_slo);
    slo = slo';
else
    slo = [];
end
 
% Display fundus image (slo)
if visu==1
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 0 scrsz(3) scrsz(4)-70]);
    subplot(1,2,1);
    imshow(slo);
end


% Read BScans, A-Scan coordinates and boundary segmentation
status = fseek(fid, 2048+(size_x_slo*size_y_slo), -1);
 
if(~any(strcmp(options, 'header')))
    bscan=zeros(n_axial, n_ascan, n_bscan, 'single');
else
    bscan= [];
end

start_x = zeros(1, n_bscan, 'double');
start_y = zeros(1, n_bscan, 'double');
end_x = zeros(1, n_bscan, 'double');
end_y = zeros(1, n_bscan, 'double');
n_seg = zeros(1, n_bscan, 'int32');
quality = zeros(1, n_bscan, 'single');
shift = zeros(1, n_bscan, 'int32');
off_seg = zeros(1, n_ascan, 'single');

ILM = zeros(n_bscan, n_ascan, 'single');
BM = zeros(n_bscan, n_ascan, 'single');
NFL_GCL = zeros(n_bscan, n_ascan, 'single');
GCL_IPL = zeros(n_bscan, n_ascan, 'single');
IPL_INL = zeros(n_bscan, n_ascan, 'single');
INL_OPL = zeros(n_bscan, n_ascan, 'single');
OPL_ONL = zeros(n_bscan, n_ascan, 'single');
ELM = zeros(n_bscan, n_ascan, 'single');
MZ_EZ = zeros(n_bscan, n_ascan, 'single');
PHROS_IDZ = zeros(n_bscan, n_ascan, 'single');
IDZ_RPE = zeros(n_bscan, n_ascan, 'single');

for i_bscan = 1:n_bscan
    ii = i_bscan - 1;
   
    % Read Bscan/Ascan positions in SLO image    
    status = fseek(fid, 16 + 2048 + (size_x_slo*size_y_slo)+(ii*(bscan_hdr_size+n_ascan*n_axial*4)), -1);
    start_x(i_bscan) = fread(fid, 1, '*double'); 
    start_y(i_bscan) = fread(fid, 1, '*double');  
    end_x(i_bscan) = fread(fid, 1, '*double');
    end_y(i_bscan) = fread(fid, 1, '*double');  
    
          
    if visu==1
        start_x_px = round(start_x(i_bscan)/scale_x_slo); %StartX in pixels
        start_y_px = round(start_y(i_bscan)/scale_y_slo); %StartY in pixels
        end_x_px = round(end_x(i_bscan)/scale_x__slo); %EndX in pixels
        end_y_px = round(end_y(i_bscan)/scale_y_slo); %EndY in pixels
        subplot(1,2,1); hold on, line([start_x_px end_x_px],[start_y_px end_y_px],'color','b');
    end
        
    % Read Bscan info 
    n_seg(i_bscan) = fread(fid, 1, '*int32'); % number of segmentated boundaries
    off_seg(i_bscan) = fread(fid, 1, '*int32'); % just in case it is useful
    quality(i_bscan) = fread(fid, 1, '*float32');
    shift(i_bscan) = fread(fid, 1, '*int32');
        
    % Read Bscan (images)    
    if(~any(strcmp(options, 'header')))
        status = fseek( fid, bscan_hdr_size + 2048 + (size_x_slo*size_y_slo) + (ii*(bscan_hdr_size+n_ascan*n_axial*4)), -1);
        oct = fread(fid, n_ascan*n_axial, '*float32');
        oct = reshape(oct, n_ascan, n_axial);
        oct = oct.^0.25;  % necessary to enhance the contrast but not sure why is that value concretely.
        bscan(:,:,i_bscan)=oct';
    end
         
    if visu==1
        subplot(1,2,2);
        imshow(bscan(:,:,i_bscan),[0 1]);
        drawnow
    end
 
    % Read segmentation     
    status = fseek(fid, 256 + 2048 + (size_x_slo*size_y_slo) + (ii*(bscan_hdr_size+n_ascan*n_axial*4)), -1 );
    seg = (fread(fid, n_seg*n_ascan, '*float' ))';
    
    ILM(i_bscan,:) = seg(1:n_ascan);
    BM(i_bscan,:) = seg((1:n_ascan) + n_ascan);
    NFL_GCL(i_bscan,:) = seg((1:n_ascan) + 2*n_ascan);
    
    if n_seg > 3
        GCL_IPL(i_bscan,:) = seg((1:n_ascan) + 3*n_ascan);
        IPL_INL(i_bscan,:) = seg((1:n_ascan) + 4*n_ascan);
        INL_OPL(i_bscan,:) = seg((1:n_ascan) + 5*n_ascan);
        OPL_ONL(i_bscan,:) = seg((1:n_ascan) + 6*n_ascan);   

        ELM(i_bscan,:) = seg((1:n_ascan) + 8*n_ascan);

        MZ_EZ(i_bscan,:) = seg((1:n_ascan) + 14*n_ascan);
        PHROS_IDZ(i_bscan,:) = seg((1:n_ascan) + 15*n_ascan);
        IDZ_RPE(i_bscan,:) = seg((1:n_ascan) + 16*n_ascan);
    end

end

% Display the segmentation
if visuseg == 1
    figure;
    mesh(double(ILM));
    hold on
    mesh(double(BM));
    if n_seg == 3
        hold on
        mesh(double(NFL_GCL));
    end
end

fclose(fid);

% Build the header
header.version = version;

header.n_ascan = n_ascan;
header.n_bscan = n_bscan;
header.n_axial = n_axial;
header.scale_x = scale_x;
header.scale_y = scale_y;
header.scale_z = scale_z;

header.size_x_slo = size_x_slo;
header.size_y_slo = size_y_slo;
header.scale_x_slo = scale_x_slo;
header.scale_y_slo = scale_y_slo;
header.fov_slo = fov_slo;
header.scan_focus = scan_focus;
header.eye = eye;
header.exam_time = exam_time;
header.scan_pattern = scan_pattern;
header.bscan_hdr_size = bscan_hdr_size;
header.id = id;
header.reference_id = reference_id;
header.pid = pid;
header.patient_id = patient_id;
header.padding = padding;
header.dob = dob;
header.vid = vid;
header.visit_id = visit_id;
header.visit_date = visit_date;
header.grid_type = grid_type;
header.grid_offset = grid_offset;
header.spare = spare;

header.quality = quality;
header.off_seg = off_seg;

header.start_x = start_x;
header.start_y = start_y;
header.end_x = end_x;
header.end_y = end_y;
header.n_seg = n_seg;
header.shit = shift;
header.off_seg = off_seg;

segment.ILM = ILM;
segment.BM = BM;
segment.NFL_GCL = NFL_GCL;
segment.GCL_IPL = GCL_IPL;
segment.IPL_INL = IPL_INL;
segment.INL_OPL = INL_OPL;
segment.OPL_ONL = OPL_ONL;
segment.ELM = ELM;
segment.MZ_EZ = MZ_EZ;
segment.PHROS_IDZ = PHROS_IDZ;
segment.IDZ_RPE = IDZ_RPE;