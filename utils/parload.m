function varargout = parload(file,varargin)
% parload - load variables from file within a parfor loop
%
% varargout = parload(filename,varargin)
%
% Input arguments:
%   file: .mat file to be open
%   varargin: name of the variables to be loaded
%
% Output arguments:
%   varargout: loaded variables
%
%  2021, Mondragon Unibertsitatea, Biomedical Engineering Department


switch (nargin-1)
    case 1
        S = load(file,varargin{1});
    case 2
        S = load(file,varargin{1},varargin{2});
    case 3
        S = load(file,varargin{1},varargin{2},varargin{3});
    case 4
        S = load(file,varargin{1},varargin{2},varargin{3},varargin{4});
    case 5
        S = load(file,varargin{1},varargin{2},varargin{3},varargin{4},...
            varargin{5});
    case 6
        S = load(file,varargin{1},varargin{2},varargin{3},varargin{4},...
            varargin{5},varargin{6});
    case 7
        S = load(file,varargin{1},varargin{2},varargin{3},varargin{4},...
            varargin{5},varargin{6},varargin{7});
    case 8
        S = load(file,varargin{1},varargin{2},varargin{3},varargin{4},...
            varargin{5},varargin{6},varargin{7},varargin{8});
    case 9
        S = load(file,varargin{1},varargin{2},varargin{3},varargin{4},...
            varargin{5},varargin{6},varargin{7},varargin{8},varargin{9});
    otherwise
        error('Wrong number of input arguments');
end

for n=1:(nargin-1)
    varargout{n} = S.(varargin{n});
end

end