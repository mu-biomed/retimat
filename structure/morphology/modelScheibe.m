function [TRT_fit,fitCoeff] = modelScheibe(rho,TRT)
% modelScheibe - fit the foveal pit model presented in Scheibe et al. ,2014
%
% [TRT_fit,fitCoeff] = modelScheibe(rho,TRT)
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


% Get layer parameters
[nAngles, nPoints] = size(TRT);

TRT_fit = nan(nAngles,nPoints);

fitCoeff = nan(nAngles,4);

% Define fitting configuration
x0 = [1 0.5 1.5 0.1]; % Inital estimation
lower = [0 0.1 0 0];
upper = [40 10 40 10];
fitfun = fittype(@(mu,sigma,gamma,alfa,x) ...
    mu*sigma^2.*x.^gamma.*exp(-mu.*x.^gamma)+alfa.*(1-exp(-mu.*x.^gamma)));    
options = fitoptions('StartPoint',x0,'Method','NonlinearLeastSquares',...
    'Lower',lower,'Upper',upper,'TolFun',1e-6,'TolX',1e-6,'MaxIter',1000,'Display','off');

% Get value at the centre
pitVal = TRT(1,1);

% Fit the model in a radial fashion
for iAngle = 1:nAngles
    
        x = rho(iAngle,:);        
        y = TRT(iAngle,:);
        
        % Normalize all points to that centre value setting the pit to 0
        y = y - pitVal;
        
        % Fit 
        [fitted,~] = fit(x',y',fitfun,options);
        fitCoeff(iAngle,:) = [fitted.mu fitted.sigma fitted.gamma fitted.alfa];

        % Reconstruct curve
        yFit = fitted(x)';

        % Add pit thickness 
        TRT_fit(iAngle,:) = yFit + pitVal;    
end

end
