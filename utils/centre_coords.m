function [X,Y,offset,ILM,BM] = centre_coords(X,Y,scanType,ILM_0,BM_0)
% CentreCoordinates - Center A-Scan coordinates so that the origin is at the scanning
%
% [X,Y,offset,ILM,BM] = centreCoordinates(X,Y,scanType,ILM_0,BM_0)
%
% Input arguments:
%   X: matrix of A-Scan x coordinates (from .vol)
%   Y: matrix of A-Scan y coordinates (from .vol)
%   scanType: raster, star or onh
%   ILM_0: matrix of ILM segmentation
%   BM_0: matrix of BM segmentation
%
% Output arguments:
%   X: matrix of A-Scan x coordinates after centering
%   Y: matrix of A-Scan y coordinates after centering
%   offset: 1x2 array of the x and y distance between the original and
%   centered
%   ILM: matrix of ILM segmentation centered
%   BM: matrix of BM segmentation centered
% 
% This function is customly designed to deal with Spectralis acquisition so
% it might not be of interest when using other equipment.
%
%  2021, Mondragon Unibertsitatea, Biomedical Engineering Department


    
switch scanType
    case 'raster'
        nBScans = size(X,1);
        nAScans = size(X,2);


        % Get central Yc
        midBScan = (nBScans+1)/2; % y nBScans=25 -> midBScan=13

        % Check nBScans is odd
        if mod(nBScans,2) ~=0
            % B-Scan and A-Scan index in the middle of the raster pattern
            Yc = mean(Y(midBScan,:));
        else
            error('Even number of B-Scans in raster');
        end

        % Get central Xc            
        midAScan = (nAScans+1)/2;
        if mod(nAScans,2) ~=0 % If odd get the central directly
            Xc = X(midBScan,midAScan);
        else % if odd there is no center point (find it between previous and posterior)
            xc_prev = X(midBScan,floor(midAScan));
            xc_post = X(midBScan,ceil(midAScan));            
            Xc = (xc_prev + xc_post)/2;
        end

        % Centre B-Scan
        X = X - Xc;
        Y = Y - Yc;

        offset = [Xc Yc];
        
        ILM = ILM_0;
        BM = BM_0;
        
    case 'star'
        nBScans = size(X,1);
        nAScans = size(X,2);
        nAngles = 2*nBScans;

        iAScan = ceil(nAScans/2); % number of A-scans in every direction

        % Remap B-Scans to a radial pattern
        Xradial = nan(nAngles,iAScan);
        Yradial = nan(nAngles,iAScan);     
        
        % Right part of each B-Scan
        Xradial(1:nBScans,:) = X(:,(end-iAScan+1):end);
        Yradial(1:nBScans,:) = Y(:,(end-iAScan+1):end);
        
        % Left part of each B-Scan (needs to be inverted)
        Xradial((nBScans+1):end,:) = X(:,iAScan:-1:1);
        Yradial((nBScans+1):end,:) = Y(:,iAScan:-1:1);  

        % Remap also layer data to a radial pattern
        ILM(1:nBScans,:) = ILM_0(:,(end-iAScan+1):end);
        ILM((nBScans+1):nAngles,:) = ILM_0(:,iAScan:-1:1);    
        
        BM(1:nBScans,:) = BM_0(:,(end-iAScan+1):end);
        BM((nBScans+1):nAngles,:) = BM_0(:,iAScan:-1:1);            

        % Get centre Xc based on the first (vertical) B-Scan            
        Xc = mean(X(1,:)); % first vertical B-Scan

        % Get centre Yc based on the horizontal B-Scan
        horBScan = nBScans/2 + 1;
        Yc = mean(Y(horBScan,:));

        X = Xradial - Xc;
        Y = Yradial - Yc;

        % Set to zero horizontal and vertical
        X(abs(X)<1e-6) = 0;
        Y(abs(Y)<1e-6) = 0;

        offset = [Xc Yc];

end
