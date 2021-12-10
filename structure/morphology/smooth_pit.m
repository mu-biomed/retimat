function Z_smooth = smooth_pit(rho, Z, method, varargin)

if nargin < 3
    error("At least 4 input arguments are expected");
end

% Check the presence of nan values
if sum(isnan(Z(:))) > 0
    warning('NaN values in layer');
end

[n_angle, n_point] = size(Z);
n_bscan = n_angle/2;

Z_smooth = nan(n_angle, n_point);

switch method
    case 'mov_average'
        n_moving = varargin{1};
        Z_smooth = moving_average(Z, n_moving);
    case 'p_spline'
        lambda = varargin{1};
        Z_smooth = pSpline(rho, Z, lambda);        
    case 'loess'
        span = varargin{1}; 

        for n=1:n_bscan
            
            % Reconstruct the B-Scan
            x = [-fliplr(rho(n+n_bscan,:)) rho(n,2:end)];
            y = [fliplr(Z(n+n_bscan,:)) Z(n,2:end)]; % Dont get center two times
            
            % Smooth the signal
            y_smooth = smooth(x, y, span/100, 'loess')';
            
            % From fitted B-Scan to radial again
            Z_smooth(n,:) = y_smooth((end-n_point+1):end);% Right part of the B-Scan
            Z_smooth(n+n_bscan,:) = y_smooth(n_point:-1:1);
        end

    case 'otherwise'
        error('Wrong model name');
end
end
