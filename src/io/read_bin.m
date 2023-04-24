function fundus = read_bin(file, im_size)
% Read grayscale slo (fundus) images from Cirrus .bin files
%
%  
% Input arguments
% --------------- 
% * **file**: path of the .bin file to read.
%
% * **im_size**: size of the slo image. If not provided it is assumed to be [664 512];
%
%                    
% Output arguments
% ---------------- 
% * **fundus**: fundus grayscale image.
%
%   
% Notes
% -----
% This code was developed by reverse engineering a few .bin files and may
% not work for every .bin image.
%
%
% Example
% -------      
% .. code-block:: matlab
%
%   I = read_bin('myfile.bin')

if nargin==1
   im_size = [664 512];
end

% Attempt to read file
fid = fopen(file);
if fid == -1
    error('Could not open the file. Check the path is correct.');
end
fundus = fread(fid, '*uint8');
fclose(fid);

fundus = reshape(fundus, im_size);

fundus = flip(fundus, 2); % swap left - right


% Reverse engineering process
% 1- get factors of the number of elements
% 2- try combinations of those
%
% n = numel(A);
% f = factors(n)
% n_x = 8*83;
% n_y = floor(n/n_x);
% n_est = n_x*n_y;
% fundus = reshape(I(1:n_est),[n_x n_y]);
%
% imshow(fundus);caxis([0 80]);