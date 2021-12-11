function Thickness = compute_thickness(seg, layers, scale_z)
%COMPUTE_THICKNESS Compute thickness for several layers
%
%   Thickness = compute_thickness(seg, layers, scale_z)
%   Compute the thickness of the layers specified by 'layers' using the
%   segmentation data stored in 'seg' struct.
%
%   Input arguments:
%  
%   'Seg'            Struct with the segmentation of the boundaries. Each field
%                    must correspond to a specific boundary.
%
%   'layers'         String or a cell array of strings indicating the layers
%                    for which to compute thickness.
%                    
%   'scale_z'        Optional. Axial (depth) resolution of the images in mm. If 
%                    provided, it is used to transform thicknes values from
%                    pixel to um units.
%  
%  
%   Output arguments:
%  
%   'Thickness'      Struct with thickness values for each layer.
%
%   
%   Notes
%   -----
%   The naming convention of the boundaries and layers might differ from one 
%   scan to the other. This function assumes the convention specified below.
%
%
%   Example
%   ---------      
%   % Compute thickness for the NFL and GCL layers
%
%   [header, seg, ~, ~] = read_vol(file,'verbose', 'coordinates');
%   Thickness = compute_thickness(seg, {'NFL','GCL'}, header.scale_z);
%
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2021

% Definition of layer and boundary names.
% Each row: [layer, top boundary, bottom boundary]
layer_top_bottom = {'TRT','ILM','BM';
    'NFL','ILM','NFL_GCL';
    'GCL','NFL_GCL','GCL_IPL';
    'IPL','GCL_IPL','IPL_INL';
    'GCIPL','NFL_GCL','IPL_INL';
    'INL','IPL_INL','INL_OPL';
    'OPL','INL_OPL','OPL_ONL';
    'ONL','OPL_ONL','ELM';
    'OPL_ONL','INL_OPL','ELM';
    'PHR1','ELM','MZ_EZ';
    'PHR2','MZ_EZ','PHROS_IDZ';
    'PHR3','PHROS_IDZ','IDZ_RPE';                
    'RPE','IDZ_RPE','BM';
    'ELM_BM','ELM','BM'};
    
if nargin == 1
    error("At least 2 input arguments must be provided");
elseif nargin == 2
    scale_z = 1;  % if not provided --> in pixels
end

if ischar(layers)
    layers = {layers}; 
end

for i=1:length(layers)
    ind = find(strcmp(layer_top_bottom(:,1), layers{i}));
    
    if length(ind)~=1 
        error(['Unknown layer. Accepted values: ',...
            'TRT,NFL,GCL,IPL,GCIPL,INL,OPL,ONL,OPL_ONL,PHR1,PHR2,PHR3,RPE,ELM_BM']);
    end
    
    % Get top and bottom layer names    
    top = layer_top_bottom{ind, 2};
    bottom = layer_top_bottom{ind, 3};
    
    if ~isfield(seg, top)
        warning(['Boundary ' top ' not found in Seg. Unable to compute thickness for ' layers{i} ' layer']);
    end
    if ~isfield(seg, bottom)
        warning(['Boundary ' bottom ' not found in Seg. Unable to compute thickness for ' layers{i} ' layer']);
    end
    
    % Compute thickness
    Thickness.(layers{i}) = 1e3 * scale_z * abs(double(seg.(top)) - double(seg.(bottom)));
end