function [TRT_fit,fitCoeffs] = modelDing(theta,rho,TRT)
% modelBreher - fit the foveal pit model presented in Ding et al. ,2014 
%
% [TRT_fit,fitCoeff] = modelDing(theta,rho,TRT)
%
% Input arguments:
%   theta,rho: matrixes with polar coordinates of TRT points
%   TRT: matrix with TRT points (each row is an angular direction)
%
% Output arguments:
%   TRT_fit: fitted TRT values
%   fitCoeffs: estimated model coefficients
%
%  2021, Mondragon Unibertsitatea, Biomedical Engineering Department


[nAngles, nPoints] = size(TRT);

TRT_fit = nan(nAngles,nPoints);

% Get cartesian coords
[X,Y] = pol2cart(theta,rho);

% Remove pit value (repeated in all radial scans)
Z = TRT(:);
Xf = X;
Yf = Y;

for n=2:nAngles
    Z(n,1) = nan; 
    Xf(n,1) = nan; 
    Yf(n,1) = nan; 
end

Z = Z(:);Z(isnan(Z))=[];
Xf = Xf(:);Xf(isnan(Xf))=[];
Yf = Yf(:);Yf(isnan(Yf))=[];

% Define fitting configuration
x0 = [360 8 -10 1.5 -8 -130 0.5*1e3 0.5*1e3]./1000;

fitfun = fittype( @(A0,A11,A12,A21,A22,K,s1,s2,x,y) ...
    A0+A11*x+A21*y+A12*x.^2+A22*y.^2 +K*exp(-x.^2/(2*s1^2)-y.^2/(2*s2^2)),...
    'independent', {'x', 'y'},'dependent', 'z' );

options = fitoptions('StartPoint',x0,'Method','NonlinearLeastSquares',...
    'TolFun',1e-6,'TolX',1e-6,'MaxIter',1000,'Display','off');

% Fit the model
fitted = fit([Xf Yf],Z,fitfun,options);
    
% Reconstruct radially
for n=1:nAngles
    TRT_fit(n,:) = fitted([X(n,:)' Y(n,:)'])';
end

fitCoeffs = [fitted.A0 fitted.A11 fitted.A12 fitted.A21 fitted.A22 fitted.K fitted.s1 fitted.s2];
