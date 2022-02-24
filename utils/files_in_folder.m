function files = files_in_folder(folder, extension)
%FILES_IN_FOLDER Return the name of the files in folder. 
%
%   files = files_in_folder(folder, extension)
%   Return files in folder with the defined extension.
%
%   Input arguments:
%  
%   'folder'         String defining the folder to examine.
%            
%   'extension'      String or cell array of strings defining the extension
%                    of the files to return. If not provided, all files in
%                    folder are returned. 
%  
%  
%   Output arguments:
%  
%   'files'          Cell array containing the name of the files in folder.
%
%   
%   Example 1
%   ---------      
%   % Get all files in folder
%
%     vol_files = files_in_folder(folder);
%     
%   Example 2
%   ---------      
%   % Get .vol files in folder
%
%     vol_files = files_in_folder(folder, '.vol');
%     
%   Example 3
%   ---------      
%   % Get both png and jpg files from folder
%
%     files = files_in_folder(folder, {'.png','.jpg'});
%     
%
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2022

folder_data = dir(folder);

if isempty(folder_data)
    error(['Folder ' folder ' not found. Check that it exists']); 
end

is_dir = vertcat(folder_data.isdir);
files = vertcat({folder_data(~is_dir).name});

if isempty(files)
    warning('No files found in folder');  
    return 
end

% If no extension given return all
if nargin == 1
   return  
end

% Get file extension
file_extension = cellfun(@get_extension, files, 'UniformOutput', false);

% Retrieve only files with the defined extension
has_extension = cellfun(@(x) any(strcmp(x, extension)), file_extension);

files = files(has_extension);

function ext = get_extension(x)
% Find dots
idx_dot = strfind(x, '.');

if isempty(idx_dot)
    ext = '';  % Files without extension  
else
    ext = x(idx_dot(end):end); % (end) in case of multiple dots
end
