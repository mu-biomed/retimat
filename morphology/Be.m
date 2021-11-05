function B = Be(t,i)
% Be - Evaluate Bernstein polynomial
%
% B = Be(t,i)
%
% Input arguments:
%   t: evaluation point (in range [0,1])
%   i: polynomial order
%
% Output arguments:
%   B: Evaluated polynomial value
%
%  2021, Mondragon Unibertsitatea, Biomedical Engineering Department

switch i
    case 0
        B = (1-t).^3;
    case 1
        B = 3.*t.*(1-t).^2;
    case 2
        B = 3.*t.^2.*(1-t);
    case 3
        B = t.^3;
end