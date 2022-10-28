function fundus = read_bin(file, im_size)
%READ_BIN Read slo images from Cirrus .bin files
%
%   fundus = read_bin(file, im_size)
%
%   Reads the grayscale slo image stored in a .bin Cirrus file.
%
%
%   Input arguments:
%  
%   'file'           Path of the .bin file to read.
%
%   'im_size'        Size of the slo image. If not provided it is assumed 
%                    to be [664 512];
%
%                    
%   Output arguments:
%  
%   'fundus'              Fundus grayscale image.
%
%   
%   Notes
%   -----
%   This code was developed by reverse engineering a few .bin files and may
%   not work for every .bin image.
%
%
%   Example
%   ---------      
%   % Read fda file
%
%     I = read_bin(file.bin)
%     
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2022


% slo size reverse engineering

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