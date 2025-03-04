function fig = generate_report(bscan, seg, fundus, layers, varargin)
% Create a summary figure with fundus, en-face reflectance map, thickness
% maps and bscans. 
%
%
% Input arguments (mandatory)
% ---------------------------
% * **bscan**:       Raw b-scans with shape [n_axial, n_ascan, n_bscan]
%
% * **seg**:         Struct with segmentation data. Each field corresponds to a different boundary. Expected shape of each field [n_bscan, n_ascan]
%
% * **fundus**:      Fundus image.
% 
% * **layers**:      List of layers for which to plot thickness maps.
%
%
% Input arguments (optional)
% --------------------------
% * **n_plot_bscan**: How many bscans to plot. Default: maximum between 5 and n_bscan.
%  
% * **n_col_max**:    Maximum number of columns in the figure. Default: 5
%
% * **file_name**:    File name to save the plot. Default: [] (plot is not saved).
%
% * **visible**:      'on' to show the figure window and 'off' to hide it. Default: 'on'
%
% * **title**:        text to us as title. Default: [] (no title).
%
%
% Output arguments
% ---------------- 
% * **fig**:           Created figure handle.
%   
%
% Notes
% -----
% This reports are useful to identify segmentation errors or poor signal
% quality.
%
%
% Example
% -------      
% .. code-block:: matlab
% 
%     [~, segment, bscan, fundus] = read_vol(file)
%     generate_report(bscan, segment, fundus, {'TRT','RNFL'})

args = parse_inputs(varargin);

% Fundus
if isempty(fundus)
    n_fundus = 0; 
else
    n_fundus = 1; 
end

% Segmentation
if isempty(seg)
    n_layer = 0;
else
    n_layer = length(layers);    
end

% Calculate the number of bscans
if isempty(bscan)
    args.n_plot_bscan = 0;
    n_en_face = 0;
else
    n_en_face = 1;
    n_bscan = size(bscan, 3);
    if args.n_plot_bscan > n_bscan
        args.n_plot_bscan = n_bscan;
        msg = ['n_plot_bscan is larger than number of b-scans in the image.', ...
               ' Only the maximum number will be plotted'];
        warning(msg);
    end
    idx_bscan = round(linspace(1, n_bscan, args.n_plot_bscan));
end

% Number of plots
n_plot = args.n_plot_bscan + n_layer + n_fundus + n_en_face;

if n_plot == 0
    warning('Not enough data to plot anything');
    return;
end

% Define figure
n_col = min([args.n_col_max n_plot]);
n_row = ceil(n_plot/n_col);

height = n_row * 300;
width  = 1500;
% width = n_col * 300;

% Plotting
fig = figure('Position',[0 0 width height],'Visible', args.visible);
tiledlayout(n_row, n_col, 'TileSpacing', 'tight');

if n_fundus == 1
    nexttile;
    if ismatrix(fundus)
        imagesc(fundus);
        colormap(gca, 'gray');
    else
        imshow(fundus);
    end
    title('Fundus');
    axis off;
end

if n_en_face == 1
    en_face = squeeze(mean(bscan, 1, 'omitnan')).';
    nexttile;
    imagesc(en_face);
    title('Reflectance');
    colormap(gca, 'gray');
    axis off;
end

if n_layer > 0
    Thick = compute_thickness(seg, layers);
    for i=1:n_layer
        nexttile;
        imagesc(Thick.(layers{i}));
        title(layers{i}, 'Interpreter', 'none');
        colormap(gca, 'jet');
        axis off;
    end
end

if args.n_plot_bscan > 0
    boundaries = fields(seg);
    n_bound = length(boundaries);
    for i=idx_bscan
        nexttile;
        imagesc(bscan(:,:,i));
        hold on;
        colormap(gca,'gray');
        axis off;
        for i_bound=1:n_bound
            plot(seg.(boundaries{i_bound})(i,:));
        end
        title(['bscan: ' num2str(i)]);
    end
end

if ~isempty(args.title)
    sgtitle(args.title, 'Interpreter','None');
end

if ~isempty(args.file_name)
    saveas(fig, args.file_name);
end

function args = parse_inputs(extra_args)

n_args = length(extra_args)/2;

if mod(n_args, 1) ~= 0
    error("Optional argument number must be even"); 
end

% Check arguments
args.n_plot_bscan = 5;
args.n_col_max    = 5;
args.file_name    = [];  
args.visible      = 'on';
args.title        = [];

for i=1:n_args
    arg = extra_args{2*i - 1};
    
    if ~ischar(arg) & ~isstring(arg)
        warning('Argument names must be char or string');
        continue;
    end
    
    if any(strcmp(arg, {'n_plot_bscan','n_col_max','file_name','visible','title'}))
        args.(arg) = extra_args{2*i};
    else        
        msg = ['Unknown argument. Valid options are: n_col_max ', ...
               'n_plot_bscan','file_name','visible','title'];
        warning(msg);
    end
end