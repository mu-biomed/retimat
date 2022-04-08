close all;clc;clear all;

addpath(genpath('C:\Users\dromero\Desktop\GITHUB\retimat'));

%% Read fundus
fname = '1000084_21011_0_0/ColorFundus.dat';

fid = fopen(fname);
fundus = fread(fid, '*uint8');

n_size = [3 2048 1536];
fundus = reshape(fundus, n_size);
fundus = permute(fundus, [3 2 1]);
fundus = flip(fundus,3);
%% Read B-Scans
fname = '1000084_21011_0_0/AlignedOCTImages.dat';

fid = fopen(fname);
I = fread(fid, '*uint16');

% Remove useless data
% I = I(2:2:end); % other is just 0s (maybe unit16)

n_size = [512 650 128];
bscan = reshape(I, n_size);
bscan = permute(bscan, [3 2 1]);

en_face = squeeze(mean(bscan,2));
%% Read segmentation
n_layer = 10;
n_ascan = 512;
n_bscan = 128;

fname = '1000084_21011_0_0/TabilSegStat.dat';
fid = fopen(fname);
seg = fread(fid,'*uint16');

n_size = [n_ascan n_bscan n_layer];

seg = reshape(seg, n_size);

seg = permute(seg, [2 1 3]);

TRT = seg(:,:,9) - seg(:,:,1);
GCIPL = seg(:,:,4) - seg(:,:,2);
%% Read meta-data
fname = '1000084_21011_0_0/FDSInfo.txt';
fid = fopen(fname,'r');
% seg = fscanf(fid,'%f',2);
meta_data = char(fread(fid)');
lines = splitlines(meta_data);

header = struct;
is_empty = cellfun(@isempty, lines);
lines(is_empty) = [];
n_line = length(lines);

idx_cat = find(cellfun(@(x) x(1) ~= ' ', lines));
n_cat = length(idx_cat);

% Loop through fields
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

                    header.fundus_dims = cellfun(@str2num, dims);
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
   
fclose(fid);

%% Read meta-data 2
fname = '1000084_21011_0_0/SegIndicators.txt';
fid = fopen(fname, 'r');
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
            
            header.macula_center = cellfun(@str2num, field_val);
        case '(Disc_Center-frame, Disc_Center-aline)'
            field_val = strrep(field_val,')','');
            field_val = strrep(field_val,'(','');
            field_val = strsplit(field_val,',');
            
            header.onh_center = cellfun(@str2num, field_val);
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

x = linspace(-header.size_x/2, header.size_x/2, n_ascan);
y = linspace(header.size_y/2, -header.size_y/2, n_bscan);

[X, Y] = meshgrid(x,y);
header.scale_z = header.size_z/(header.n_axial-1);

% Foveal center
x_cen = X(1,header.macula_center(2));
y_cen = Y(header.macula_center(1),1);

% Align the scans
X_cen = X - x_cen;
Y_cen = Y - y_cen;

% Pixel to um
TRT = double(TRT) * header.scale_z * 1e3;
GCIPL = double(GCIPL) * header.scale_z * 1e3;
%% Plotting
tiledlayout(3,3,'TileSpacing','compact');

% Fundus
nexttile;
imagesc(fundus);
hold on;
x = [header.oct_fundus_ul(1) header.oct_fundus_br(1)];
y = [header.oct_fundus_ul(2) header.oct_fundus_br(2)];

scatter(x,y,'g');
plot([x(1) x(2) x(2) x(1) x(1)],[y(1) y(1) y(2) y(2) y(1)],'--g');
axis off;

% En face
nexttile;
imagesc(en_face);
colormap(gca,gray);
hold on;
scatter(header.macula_center(2), header.macula_center(1),15,'r','filled');
axis off;
title('En face');

% B-scan
% i_bscan = round(linspace(1,n_bscan,8));
i_bscan = round(n_bscan/2);
for i=1:length(i_bscan)
    nexttile;

    imagesc(squeeze(bscan(i_bscan(i),:,:)));hold on;
    for i_layer=1:n_layer
        plot(squeeze(seg(i_bscan(i),:,i_layer)));
    end
    colormap(gca,gray);
    axis off;
    title(['Bscan:' num2str(i_bscan(i))]);
end 

% TRT
nexttile;
surf(X,Y,TRT,'EdgeColor','none');view(0,90);
hold on;
scatter3(x_cen, y_cen, max(TRT(:)), 15, 'r', 'filled');
axis([-3 3 -3 3]);
title('TRT');
colorbar;

% TRT ETDRS
nexttile;
[TRT_etdrs, sect] = sectorize_map(X_cen,Y_cen,TRT,'mean','etdrs');
plot_sectors(TRT_etdrs, sect);
title('TRT ETDRS');
axis([-3 3 -3 3]);
axis on;

% TRT rings
nexttile;
[TRT_ring, sect] = sectorize_map(X_cen,Y_cen,TRT,'mean','ring',linspace(0,3,10));
plot_sectors(TRT_ring, sect);
title('TRT rings');
axis([-3 3 -3 3]);
axis on;

% GCIPL
nexttile;
surf(X,Y,GCIPL,'EdgeColor','none');view(0,90);
hold on;
scatter3(x_cen, y_cen, max(TRT(:)), 15, 'r', 'filled');
axis([-3 3 -3 3]);
title('GCIPL');
colorbar;

% GCIPL ETDRS
nexttile;
[GCIPL_etdrs, sect] = sectorize_map(X_cen,Y_cen,GCIPL,'mean','etdrs');
plot_sectors(GCIPL_etdrs, sect);
title('GCIPL ETDRS');
axis([-3 3 -3 3]);
axis on;

% GCIPL rings
nexttile;
[GCIPL_ring, sect] = sectorize_map(X_cen,Y_cen,GCIPL,'mean','ring',linspace(0,3,10));
plot_sectors(GCIPL_ring, sect);
title('GCIPL rings');
axis([-3 3 -3 3]);
axis on;

sgtitle('1000084_21011_0_0','Interpreter','none');