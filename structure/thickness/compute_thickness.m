function Thickness = compute_thickness(header,layers)
% get_thickness - Compute thickness from from .vol data
%
% Thickness = compute_thickness(header,layers)
%
% Input arguments:
%   header: struct with .vol header
%   layers: cell array specifying the layers to compute
%
% Output arguments:
%   Thickness: struct with a field for each layer
%
% David Romero-Bascones 
% dromero@mondragon.edu
% 2021, Mondragon Unibertsitatea, Biomedical Engineering Department
% ------------------------------------------------------------------------

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
    

% Get necessary header parameters
sizeZ  = double(header.SizeZ);
scaleZ = double(header.ScaleZ);

for i=1:length(layers)
    ind = find(strcmp(layer_top_bottom(:,1),layers{i}));
    
    if length(ind)~=1 
        error('wrong number of layers detected');
    end
    
    % Get top and bottom layer names    
    top = layer_top_bottom{ind,2};
    bottom = layer_top_bottom{ind,3};
    
    % Get segmented boundaries in um
    top_thick = (sizeZ - double(header.(top)))*scaleZ;
    bottom_thick = (sizeZ - double(header.(bottom)))*scaleZ;
    
    % Compute thickness
    Thickness.(layers{i}) = 1e3*(top_thick - bottom_thick);
end