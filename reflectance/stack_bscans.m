function stacked = stack_bscans(bscan, seg, layers)
%stacked = STACK_BSCANS(bscan, layers)
%
%   stacked = stack_bscans(bscan, seg, layers)
%   Stack the bscan portion showing the same layer into the a single image as 
%   proposed in [1].
%
%   Input arguments:
%  
%   'bscan'          Matrix with B-scans of dimensions n_axial x n_ascan x 
%                    n_bscan.
%            
%   'seg'            Struct with segmentation data.
%
%   'layers'         String or a cell array of strings with the name of the
%                    layer to be stacked. Note that, depending on the layer of
%                    interest, the segmentation must follow a certain naming
%                    convention (see the beggining of the function).
%  
%
%   Output arguments:
%  
%   'stacked'        Struct with a stacked b-scan for each layer in layers.          
%  
%
%   
%   Notes
%   -----
%   The stacking process is sensitive to NaN values in the segmentation. This
%   function performs a rough extrapolation of NaN values that may not work if
%   the number of NaN values is high.
%
%
%   References
%   ----------
%   [1] Tazarjani, Retinal OCT Texture Analysis for Differentiating Healthy
%   Controls from Multiple Sclerosis (MS) with/without Optic Neuritis,
%   BioMed Research International, 2021. https://doi.org/10.1155/2021/5579018
%
%
%   Example
%   ---------      
%   % Generate stacked b-scans of RNFL and GCIP layers
%
%     [~, seg, bscan] = read_vol(file);
%     layers = {'RNFL', 'GCIP'};
%     stacked = stack_bscans(bscan, seg, layers)
%     
%
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2022


layer_top_bottom = {'RNFL',     'ILM',        'RNFL_GCL';
                    'GCL',      'RNFL_GCL',   'GCL_IPL';
                    'IPL',      'GCL_IPL',    'IPL_INL';
                    'INL',      'IPL_INL',    'INL_OPL';
                    'OPL',      'INL_OPL',    'OPL_ONL';
                    'ONL',      'OPL_ONL',    'ELM';
                    'MZ',       'ELM',        'MZ_EZ';
                    'EZ',       'MZ_EZ',      'EZ_OSP';
                    'OSP',      'EZ_OSP',     'OSP_IZ';
                    'IZ',       'OSP_IZ',     'IZ_RPE';                
                    'RPE',      'IZ_RPE',     'BM';
                    
                    % Composite layers
                    'TRT',      'ILM',        'BM';
                    'GCIP',     'RNFL_GCL',   'IPL_INL';
                    'GCC',      'ILM',        'IPL_INL';   
                    'IRL',      'ILM',        'IPL_INL';
                    'ONPL',     'INL_OPL',    'ELM';
                    'EZOSP',    'MZ_EZ',      'OSP_IZ';
                    'ELM_BM',   'ELM',        'BM'};
                

if ischar(layers)
    if strcmp(layers, 'all')
        layers = layer_top_bottom(:,1);
    else
        layers = {layers};
    end
end

n_layer = length(layers);

if ndims(bscan) ~=3
    error('Input bscan must be a 3D volume');
end

[~, n_ascan, n_bscan] = size(bscan);

% Grid in case we need to extrapolate segmentation
[X,Y] = meshgrid(linspace(0,1,n_ascan),linspace(0,1,n_bscan));

stacked = struct;
for i_layer=1:n_layer
    layer = layers{i_layer};
    
    ind = find(strcmp(layer_top_bottom(:,1), layer));
    
    if length(ind)~=1 
        error(['Unknown layer:' layer '. Accepted values: ' strjoin(layer_top_bottom(:,1))]);
    end
    
    % Get top and bottom layer names    
    top = layer_top_bottom{ind, 2};
    bottom = layer_top_bottom{ind, 3};
    
    if ~isfield(seg, top)
        warning(['Boundary ' top ' not found in segmentation. Unable to compute thickness for ' layer ' layer']);
        continue
    elseif ~isfield(seg, bottom)
        warning(['Boundary ' bottom ' not found in segmentation. Unable to compute thickness for ' layer ' layer']);
        continue
    end
    
    % Get segmentation
    z_top = double(seg.(top));
    z_bottom = double(seg.(bottom));
    
    % Remove NaNs
    mask_nan = isnan(z_top);
    if any(mask_nan(:)) 
        warning(['NaN values in the segmentation of ' top ' layer.']);           
        interpol = scatteredInterpolant(X(~mask_nan),Y(~mask_nan),z_top(~mask_nan));
        z_top = reshape(interpol(X(:),Y(:)),size(X));
    end
    
    mask_nan = isnan(z_bottom);
    if any(mask_nan(:)) 
        warning(['NaN values in the segmentation of ' bottom ' layer.']);         
        interpol = scatteredInterpolant(X(~mask_nan),Y(~mask_nan),z_bottom(~mask_nan));
        z_bottom = reshape(interpol(X(:),Y(:)),size(X));
    end
    
    % Round to indexes
    z_top = round(z_top);
    z_bottom = round(z_bottom);
    
    % If last layer we can get one pixel more otherwise we need -1 to avoid
    % including the same pixel in a more than 1 layer when a boundary is the
    % top and bottom boundary of two different layers. Not a big deal probably.
    if strcmp(bottom, 'BM')
        idx = 0;
    else
        idx = 1;
    end
    
    % Stack b-scans (loop through all b-scans for each a-scan)    
    stack = cell(1,n_ascan);
    for a=1:n_ascan
        stack{a} = [];
        for b=1:n_bscan
            stack{a} = [stack{a} ; bscan(z_top(b,a):z_bottom(b,a)-idx,a,b)];
        end
    end
    
    % Get maximum height (number of rows of the end stacked image)
    max_height = max(cellfun(@length, stack));
    
    % Padd columns with NaN and obtain a matrix
    stack = cellfun(@(x) [nan(max_height-length(x),1) ; x], stack, ...
                    'UniformOutput', false);
    stacked.(layer) = cell2mat(stack);
end                
