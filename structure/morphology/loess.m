function [TRT_fit,fitparams] = loess(rho,TRT,span)
% loess - smooth the TRT curves based on LOESS smoothing
%
% [TRT_fit,fitparams] = loess(rho,TRT,span)
%
% Input arguments:
%   rho: matrix with rho coordinates
%   TRT: matrix with TRT points (each row is an angular direction)
%   span: % of points used for smoothing in LOESS
%
% Output arguments:
%   TRT_fit: smoothed TRT values
%   fitparams: span value
%
%  2021, Mondragon Unibertsitatea, Biomedical Engineering Department

[nAngles, nPoints] = size(TRT);

nBScans = nAngles/2;

TRT_fit = nan(nAngles,nPoints);

for n=1:nBScans

    % Reconstruct the B-Scan
    x = [-fliplr(rho(n+nBScans,:)) rho(n,2:end)];
    y = [fliplr(TRT(n+nBScans,:)) TRT(n,2:end)]; % Dont get center two times

    % Smooth the signal    
    ySmooth = smooth(x,y,span/100,'loess')';

    % From fitted B-Scan to radial again
    TRT_fit(n,:) = ySmooth((end-nPoints+1):end);% Right part of the B-Scan
    TRT_fit(n+nBScans,:) = ySmooth(nPoints:-1:1);
end

fitparams = span;
