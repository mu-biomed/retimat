function I = read_bin(file, size)
% Function that reads slo images from Cirrus .bin files

% slo size reverse engineering

if nargin==1
   size = [664 512];
end

% Attempt to read file
fid = fopen(file);
if fid == -1
    error('Could not open the file. Check the path is correct.');
end
I = fread(fid, '*uint8');
fclose(fid);

I = reshape(I, size);


% Reverse engineering process
% 1- get factors of the number of elements
% 2- try combinations of those
%
% n = numel(A);
% f = factors(n)
% n_x = 8*83;
% n_y = floor(n/n_x);
% n_est = n_x*n_y;
% I = reshape(I(1:n_est),[n_x n_y]);
%
% imshow(I);caxis([0 80]);