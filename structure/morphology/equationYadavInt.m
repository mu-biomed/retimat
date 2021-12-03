function yf = equationYadavInt(x,alfa,beta,P0,P3)
% equationYadavInt - evaluate the model presented by Yadav et al. (2017)
% to model de inner part of the B-Scan (foveal center to rim)
%
% yf = equationYadavInt(x,alfa,beta,P0,P3)
%
% Input arguments:
%   alfa,beta,P0,P3: model coefficients (defining Bezier curves)
%   x: evaluating point
%
% Output arguments:
%   yf: evaluated value
%
%  2021, Mondragon Unibertsitatea, Biomedical Engineering Department

t = linspace(0,1,length(x));

T = [1;0];

P1 = P0 + alfa*T;
P2 = P3 - beta*T;

Q = P0*Be(t,0) + P1*Be(t,1) + P2*Be(t,2) + P3*Be(t,3);

yf = interp1(Q(1,:),Q(2,:),x);