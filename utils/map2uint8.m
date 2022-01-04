function I = map2uint8(Z, gray_limits)
%MAP2UINT8 Convert a distribution or matrix into a uint8 format
%
%   Z = map2uint8(Z)
%   Transforms the original range into [0,255] 
%
%   Input arguments:
%  
%   'Z'              Original 2D map (usually double).
%            
%   'gray_limits'    Array with minimum and maximum levels of Z to be mapped to
%                    0 and 255 correspondly. If not provided, then the minimum
%                    and maximum values are used.
%  
%  
%   Output arguments:
%  
%   'I'              2D map transformed into uint8
%
%   
%   Example
%   ---------      
%   % Conversion of a random 2D map
%
%     Z = randn(1000,1000);
%     I = map2uint8(Z);
%     
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2021

if nargin == 1
    gray_limits = [min(Z(:)) max(Z(:))];
elseif nargin == 2
    if length(gray_limits) ~= 2
        error('gray_limits is expected to be an array with 2 elements'); 
    end
else 
    error('Number of input arguments must be 1 or 2'); 
end

I = (Z - gray_limits(1))/(gray_limits(2) - gray_limits(1));  %  Normalize [0,1]

I = uint8(round(255*I));  %  Normalize [0, 255]



