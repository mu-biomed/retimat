function [TRT_fit,fitCoeffs] = modelYadav(rho,TRT)
% modelYadav - fit the foveal pit model presented in Yadav et al. ,2017
%
% [TRT_fit,fitCoeff] = modelYadav(rho,TRT)
%
% Input arguments:
%   rho: matrix with rho coordinates
%   TRT: matrix with TRT points (each row is an angular direction)
%
% Output arguments:
%   TRT_fit: fitted TRT values
%   fitCoeffs: estimated model coefficients
%
%  2021, Mondragon Unibertsitatea, Biomedical Engineering Department


[nAngles, nPoints] = size(TRT);

TRT_fit = nan(nAngles,nPoints);

fitCoeffs = nan(nAngles,7);

% Interior segment equation
fitfun.int = fittype('equationYadavInt(x,alfa,beta,P0,P3)',...
    'coefficients',{'alfa','beta'},'problem',{'P0','P3'});

% Exterior segment fitting configuration
fitfun.ext = fittype('equationYadavExt(x,alfa,P2x,P2y,P0,P3)',...
    'coefficients',{'alfa','P2x','P2y'},'problem',{'P0','P3'});

for iAngle=1:nAngles

    x = rho(iAngle,:);
    y = TRT(iAngle,:);
    
    % Get Max
    [~,indmax] = max(y);
    indmax = indmax(1);
    
    % Fit interior
    xInt = x(1:indmax);
    yInt = y(1:indmax);

    P0 = [xInt(1) ;yInt(1)];
    P3 = [xInt(end) ;yInt(end)];

    % initialize alfa,beta to half width
    alfa0 = (P3(1)-P0(1))/2;
    beta0 = (P3(1)-P0(1))/2;
    x0 = [alfa0 beta0]; 
    lower = [0 0];
    upper = [2 2];

    options.int = fitoptions('StartPoint',x0,'Method','NonlinearLeastSquares',...
        'Lower',lower,'Upper',upper,'TolFun',1e-6,'TolX',1e-6,'MaxIter',1000,'Display','off');
    
    fitted = fit(xInt',yInt',fitfun.int,options.int,'problem',{P0,P3});
    yIntFit = fitted(xInt)';
    
    % Fit exterior  
    xExt = x(indmax:end);
    yExt = y(indmax:end);

    P0 = [xExt(1) ;yExt(1)];
    P3 = [xExt(end) ;yExt(end)];
    
    alfa0 = (P3(1)-P0(1))/2; % half width of the exterior segment
    P2x0 = P3(1)-(P3(1)-P0(1))/2; % P3 - half width
    P2y0 = (P0(2)+P3(2))/2; % in the middle between P0y and P3y
    x0 = [alfa0  P2x0 P2y0]; 
    lower = [0 0 0.200];
    upper = [2 3 0.400];
 
    options.ext = fitoptions('StartPoint',x0,'Method','NonlinearLeastSquares',...
        'Lower',lower,'Upper',upper,'Robust','LAR','TolFun',1e-6,'TolX',1e-6,'MaxIter',1000,'Display','off');

    fitted = fit(xExt',yExt',fitfun.ext,options.ext,'problem',{P0,P3});
    yExtFit = fitted(xExt)';
    
    % Force first point to have the highest thickness
    yExtFit(yExtFit>yExtFit(1)) = yExtFit(1);
    
    % Recombine segments
    yFit = [yIntFit(1:end-1) yExtFit];

    % Store 
    TRT_fit(iAngle,:) = yFit;
    
end





