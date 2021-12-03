function [layerFit,fitCoeff] = fit_pit_model(theta,rho,TRT,pitModel)
% fitPitModel - Adjust a mathematical foveal pit model to TRT data
%
% [layerFit,fitCoeff] = fit_pit_model(theta,rho,TRT,pitModel)
%
% Input arguments:
%   theta, rho: matrixes with radial coordinates of TRT curve points
%   TRT: matrix with TRT curve points
%   pitModel: string defining the model to be used
%
% Output arguments:
%   layerFit: matrix with fitted TRT curves
%   fitCoeff: struct with fitted coefficients of the model
%
%  2021, Mondragon Unibertsitatea, Biomedical Engineering Department


% Check the presence of nan values
if sum(isnan(TRT(:))) > 0
    warn('NaN values in layer');
end

% Moving average
if contains(pitModel,'movAvg')
    ind = strfind(pitModel,'_');
    nMoving = str2double(pitModel(ind+1:end));    
    [layerFit,fitCoeff] = movingAverage(TRT,nMoving);
    return;
end

% P-splines
if contains(pitModel,'pSpline')
    ind = strfind(pitModel,'_');
    lambda = str2double(pitModel(ind+1:end));    
    [layerFit,fitCoeff] = pSpline(rho,TRT,lambda);
    return;
end

% Loess 
if contains(pitModel,'loess')
    ind = strfind(pitModel,'_');
    span = str2double(pitModel(ind+1:end));    
    [layerFit,fitCoeff] = loess(rho,TRT,span);
    return;
end

switch pitModel
   
    case 'Scheibe'
        TRT = 1e-3*TRT;
        [layerFit,fitCoeff] = modelScheibe(rho,TRT); 
        layerFit = 1e3*layerFit;
    case 'Yadav'
        [layerFit,fitCoeff] = modelYadav(rho,TRT);
    case 'Ding'
        [layerFit,fitCoeff] = modelDing(theta,rho,TRT);
    case 'Breher'
        [layerFit,fitCoeff] = modelBreher(rho,TRT);
    case 'Dubis'
        [layerFit,fitCoeff] = modelDubis(rho,TRT);
    case 'Liu'
        [layerFit,fitCoeff] = modelLiu(rho,TRT);
    case 'none'
        layerFit = TRT;
        fitCoeff = nan;
    case 'otherwise'
        error('Wrong model name');
end
