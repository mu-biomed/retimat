function parsave(filename,varargin)
% parload - save variables in file within a parfor loop
%
% parsave(filename,varargin)
%
% Input arguments:
%   filename: .mat file to save variables
%   varargin: name of the variables to be saved
%
%  2021, Mondragon Unibertsitatea, Biomedical Engineering Department


varNames = cell(1,nargin-1);

% Get input variable name and create them again with their value
for n=1:(nargin-1)
    varNames{n} = inputname(n+1);
    eval([varNames{n},' = varargin{n};']); % Asign variable value to a variable with the same name
end

% Save data (switch cause I didn't know how to concatenate variables)
switch (nargin-1)
    case 1
        save(filename,varNames{1});
    case 2
        save(filename,varNames{1},varNames{2});
    case 3
        save(filename,varNames{1},varNames{2},varNames{3});
    case 4
        save(filename,varNames{1},varNames{2},varNames{3},varNames{4});
    case 5
        save(filename,varNames{1},varNames{2},varNames{3},varNames{4},...
            varNames{5});
    case 6
        save(filename,varNames{1},varNames{2},varNames{3},varNames{4},...
            varNames{5},varNames{6});
    case 7
        save(filename,varNames{1},varNames{2},varNames{3},varNames{4},...
            varNames{5},varNames{6},varNames{7});
    case 8
        save(filename,varNames{1},varNames{2},varNames{3},varNames{4},...
            varNames{5},varNames{6},varNames{7},varNames{8});
    case 9
        save(filename,varNames{1},varNames{2},varNames{3},varNames{4},...
            varNames{5},varNames{6},varNames{7},varNames{8},varNames{9});
    otherwise
        error('Wrong number of input arguments');
end

end