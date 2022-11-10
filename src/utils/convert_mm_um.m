function [Z, unit_in] = convert_mm_um(Z, unit)
%CONVERT_MM_UM Convert mm to um and viceversa in a safe way
%
%   [Z, unit_in] = convert_mm_um(Z, unit)
%   Transform the values in Z to the units define by 'unit'
%
%   Input arguments:
%  
%   'Z'              Input values in either mm or um.          
%  
%   'unit'           Intended output unit.
%                    Options: ['mm', 'um']
%  
%  
%   Output arguments:
%  
%   'Z'              Output values in the unit defined by 'unit'
%
%   'unit_in'        Unit of the provided inpyt Z. 'mm' or 'um'
%
%
%   Example
%   ---------      
%   % Convert thickness from mm to um
%
%   [header, seg, ~, ~] = read_vol(file);
%   Thickness_mm = compute_thickness(seg, 'TRT', header.scale_z);
%   Thickness_um = convert_mm_um(Thickness_mm, 'um');
%
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2021

% Get initial units
if mean(abs(Z(:)), 'omitnan') > 10 % um
    unit_in = 'um';
else
    unit_in = 'mm';
end

% Convert units
if isequal(unit_in, 'mm') & isequal(unit, 'um')
    Z = 1e3*Z;
elseif isequal(unit_in, 'um') & isequal(unit, 'mm')
    Z = 1e-3*Z;
end
end