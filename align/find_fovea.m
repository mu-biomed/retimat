function [x_fovea,y_fovea] = find_fovea(X,Y,TRT,interpMethod,method)
% find_fovea - Align the TRT thickness map by locating the foveal center
% 
% [x_fovea,y_fovea] = find_fovea(X,Y,TRT,interpMethod,alignMethod)
%
% Input arguments:
%   X: matrix of A-Scan x coordinates
%   Y: matrix of A-Scan y coordinates
%   TRT: matrix of TRT values corresponding to each coordinate
%   interpMethod: method used for interpolating in resample + min and 
% smooth + min
%   alignMethod: method used to locate the foveal center
%
% Output arguments:
%   x_fovea: x coordinate of the located foveal center
%   y_fovea: y coordinate of the located foveal center
%
%  IMPORTANT: the smooth+min method is based on the foveafinder.m function
%  of AURA Tools. If you use it, please provide appropriate credit to the
%  original work (https://www.nitrc.org/projects/aura_tools/).
%
%  2021, Mondragon Unibertsitatea, Biomedical Engineering Department


switch method
    case 'none'     % default (without any alignemnt)
        x_fovea = 0;
        y_fovea = 0;
    
    case 'min'
        % Search only in the central region
        maxR = 0.85; 
        R = sqrt(X.^2+Y.^2);
        maskCentral = R <= maxR;
        TRT_aux = TRT;
        TRT_aux(~maskCentral) = nan;
        
        % get min
        [~,indmin] = min(TRT_aux(:));
        [indx,indy] = ind2sub(size(X),indmin);        
        x_fovea = X(indx,indy);
        y_fovea = Y(indx,indy);        
     case 'resample_min'
        % Interpolate to small regular grid 
        maxD = 0.85; 
        step = 0.02; 
        [Xg,Yg] =meshgrid(-maxD:step:maxD,-maxD:step:maxD);        
                
        TRT_grid = reshape(griddata(X(:),Y(:),TRT(:),Xg(:),Yg(:),...
            interpMethod),size(Xg));
  
        [~,indmin] = min(TRT_grid(:));
        [indx, indy] = ind2sub(size(TRT_grid),indmin);
        x_fovea = Xg(indx,indy);
        y_fovea = Yg(indx,indy);   
    case 'smooth_min'
        % Based on AURA foveaFinder.m function
        
        % Interpolate to small regular grid 
        maxD = 0.85; % 
        step = 0.02; 
        [Xg,Yg] =meshgrid(-maxD:step:maxD,-maxD:step:maxD);        
                
        TRT_grid = reshape(griddata(X(:),Y(:),TRT(:),Xg(:),Yg(:),...
            interpMethod),size(Xg));
        
        % Minimum thickness average in 0.05 mm radius circle
        maxR = 0.05;
        R = Xg.^2+Yg.^2;
        m = R < maxR;
        TRT_filt = imfilter(TRT_grid,double(m)/sum(m(:)),'replicate');

        [~,indmin] = min(TRT_filt(:));
        [indx, indy] = ind2sub(size(TRT_filt),indmin);
        x_fovea = Xg(indx,indy);
        y_fovea = Yg(indx,indy);

    otherwise
        error('Unsupported alignment method');
end