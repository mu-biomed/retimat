function parsave(file,varargin)
% parload - save variables in file within a parfor loop
%
% parsave(filename,varargin)
%
% Input arguments:
%   filename: .mat file to save variables
%   varargin: name of the variables to be saved
%
%  2021, Mondragon Unibertsitatea, Biomedical Engineering Department


var_names = cell(1,nargin-1);

% Get input variable name and create them again with their value
for n=1:(nargin-1)
    var_names{n} = inputname(n+1);
    eval([var_names{n},'=varargin{n};']); % Asign variable value to a variable with the same name
end

% Save data (switch cause I didn't know how to concatenate variables)
switch (nargin-1)
    case 1
        save(file,var_names{1});
    case 2
        save(file,var_names{1},var_names{2});
    case 3
        save(file,var_names{1},var_names{2},var_names{3});
    case 4
        save(file,var_names{1},var_names{2},var_names{3},var_names{4});
    case 5
        save(file,var_names{1},var_names{2},var_names{3},var_names{4},...
            var_names{5});
    case 6
        save(file,var_names{1},var_names{2},var_names{3},var_names{4},...
            var_names{5},var_names{6});
    case 7
        save(file,var_names{1},var_names{2},var_names{3},var_names{4},...
            var_names{5},var_names{6},var_names{7});
    case 8
        save(file,var_names{1},var_names{2},var_names{3},var_names{4},...
            var_names{5},var_names{6},var_names{7},var_names{8});
    case 9
        save(file,var_names{1},var_names{2},var_names{3},var_names{4},...
            var_names{5},var_names{6},var_names{7},var_names{8},var_names{9});
    otherwise
        error('Wrong number of input arguments');
end

end