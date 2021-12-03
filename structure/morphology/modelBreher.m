function  [TRT_fit,fitCoeff] = modelBreher(rho,TRT)
% modelBreher - fit the foveal pit model presented in Breher et al. ,2019 
%
% [TRT_fit,fitCoeff] = modelBreher(rho,TRT)
%
% Input arguments:
%   rho: matrix with rho coordinates
%   TRT: matrix with TRT points (each row is an angular direction)
%
% Output arguments:
%   TRT_fit: fitted TRT values
%   fitCoeff: estimated model coefficients
%
%  2021, Mondragon Unibertsitatea, Biomedical Engineering Department

[nAngles, nPoints] = size(TRT);

TRT_fit = nan(nAngles,nPoints);

nBScans = nAngles/2;

fitCoeff = nan(nBScans,9);

% Fitting configuration
x0 = [0.3 -0.05 -0.02 0.3 0.04 -0.47 4.8 0.4 0.6];

lower = [0 -Inf -Inf -0.5 -0.5 -1 -Inf -Inf -Inf];
upper = [1 Inf Inf    0.5  0.5  3 Inf Inf Inf];
fitfun = fittype(@(a1,a2,a3,b1,b2,b3,c1,c2,c3,x) ...
    a1.*exp(-((x-b1)./c1).^2) + a2.*exp(-((x-b2)./c2).^2) + a3.*exp(-((x-b3)./c3).^2));
options = fitoptions('StartPoint',x0,'Method','NonlinearLeastSquares',...
    'Lower',lower,'Upper',upper,'TolFun',1e-6,'TolX',1e-6,'MaxIter',1000,'Display','off');

    
for n=1:nBScans
    
    x = [-fliplr(rho(n+nBScans,:)) rho(n, 2:end)];
    y = [fliplr(TRT(n+nBScans,:)) TRT(n, 2:end)]; % Dont get center two times
    
    % Fit the model
    fitted = fit(x',y',fitfun,options);
    
    % Reconstruct fitted curve
    yFit = fitted(x)';
    
    % From fitted B-Scan to radial again    
    TRT_fit(n,:) = yFit((end-nPoints+1):end);% Right part of the B-Scan
    TRT_fit(n+nBScans,:) = yFit(nPoints:-1:1);
    
    % Store Coefficients
    fitCoeff(n,:) = [fitted.a1 fitted.a2 fitted.a3 fitted.b1 fitted.b2 fitted.b3 fitted.c1 fitted.c2 fitted.c3];
end
