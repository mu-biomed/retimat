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
%   The naming convention of the boundaries and layers might differ both in
%   the literature and in different segmentation algorithms. Here we mostly
%   follow the APOSTEL 2.0 recommendations ([1]) as specified below:
%    - BM:    Bruch's membrane
%    - ELM:   external limiting membrane
%    - EZ:    ellipsoid zone
%    - GCL:   ganglion cell layer
%    - GCIP:  ganglion cell and inner plexiform layer (composite)
%    - ILM:   inner limiting membrane
%    - INL:   inner nuclear layer 
%    - IPL:   inner plexiform layer
%    - IRL:   inner retinal layers (composite)
%    - IZ:    interdigitation zone
%    - MZ:    myoid zone
%    - ONL:   outer nuclear layer
%    - ONPL:  outer nuclear - plexiform layer (composite)
%    - OPL:   outer plexiform layer
%    - OSP:   outer segment of the photoreceptors
%    - RNFL:  retinal nerve fiber layer
%    - RPE:   retinal pigment epithelium
%   For boundaries and layer composites not specified in [1] we have used
%   underscores to denote top/bottom boundaries. To consider the following
%   naming case:
%   - EZ_OSP: boundary between EZ and OSP
%   - EZOSP: composite layer with EZ + OSP
%
%
%   References
%   ----------
%   [1] Aytulun et al., "APOSTEL 2.0 Recommendations for Reporting 
%   Quantitative Optical Coherence Tomography Studies", Neurology, 2021
%   https://doi.org/10.1212/WNL.0000000000012125 
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
% Each row:         Layer      Top boundary  Bottom boundary
%
%                   Single layers
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
                    'IRL',      'ILM',        'IPL_INL';
                    'ONPL',     'INL_OPL',    'ELM';
                    'EZOSP',   'MZ_EZ',      'OSP_IZ';
                    'ELM_BM',   'ELM',        'BM'};

if nargin == 1
    error("At least 2 input arguments must be provided");
elseif nargin == 2
    scale_z = 1;  % if not provided --> in pixels
end

if ischar(layers)
    if strcmp(layers, 'all')
        layers = layer_top_bottom(:,1);
    else
        layers = {layers};
    end
end

Thickness = struct;

for i=1:length(layers)
    ind = find(strcmp(layer_top_bottom(:,1), layers{i}));
    
    if length(ind)~=1 
        error(['Unknown layer:' layers{i} '. Accepted values: ' strjoin(layer_top_bottom(:,1))]);
    end
    
    % Get top and bottom layer names    
    top = layer_top_bottom{ind, 2};
    bottom = layer_top_bottom{ind, 3};
    
    if ~isfield(seg, top)
        warning(['Boundary ' top ' not found in segmentation. Unable to compute thickness for ' layers{i} ' layer']);
        continue
    elseif ~isfield(seg, bottom)
        warning(['Boundary ' bottom ' not found in segmentation. Unable to compute thickness for ' layers{i} ' layer']);
        continue
    end
    
    % Compute thickness
    Thickness.(layers{i}) = 1e3 * scale_z * abs(double(seg.(top)) - double(seg.(bottom)));
end