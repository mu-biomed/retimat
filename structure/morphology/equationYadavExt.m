function yf = equationYadavExt(x,alfa,P2x,P2y,P0,P3)
% equationYadavExt - evaluate the model presented by Yadav et al. (2017)
% to model de external part of the B-Scan (beyond rim)
%
% yf = equationYadavExt(x,alfa,P2x,P2y,P0,P3)
%
% Input arguments:
%   alfa,P2x,P2Y,P0,P3: model coefficients (defining Bezier curves)
%   x: evaluating point
%
% Output arguments:
%   yf: evaluated value
%
%  2021, Mondragon Unibertsitatea, Biomedical Engineering Department

t = linspace(0,1,length(x));

T = [1;0];

P1 = P0 + alfa*T;
P2 = [P2x;P2y];

Q = P0*Be(t,0) + P1*Be(t,1) + P2*Be(t,2) + P3*Be(t,3);
yf = interp1(Q(1,:),Q(2,:),x);
