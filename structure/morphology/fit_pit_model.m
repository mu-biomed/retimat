function [Z_fit, fit_coeff] = fit_pit_model(theta, rho, Z, pit_model, varargin)
% fitPitModel - Adjust a mathematical foveal pit model to foveal pit surface
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

if nargin < 4
    error("At least 4 input arguments are expected");
end

% Check the presence of nan values
if sum(isnan(Z(:))) > 0
    warning('NaN values in layer');
end

[n_angle, n_point] = size(Z);
n_bscan = n_angle/2;

Z_fit = nan(n_angle, n_point);

switch pit_model
   
    case 'Ding'
        [X,Y] = pol2cart(theta,rho);
        
        % Remove pit value (repeated in all radial scans)
        Z = Z(:);
        Xf = X;
        Yf = Y;
        
        for n=2:n_angle
            Z(n,1) = nan;
            Xf(n,1) = nan;
            Yf(n,1) = nan;
        end
        
        Z = Z(~isnan(Z)); 
        Xf = Xf(~isnan(Xf));
        Yf = Yf(~isnan(Yf)); 
        
        % Define fitting configuration
        x0 = [360 8 -10 1.5 -8 -130 0.5*1e3 0.5*1e3]./1000;
        
        fun = fittype( @(A0, A11, A12, A21, A22, K, s1, s2, x, y) ...
            A0 + A11*x + A21*y +A12*x.^2 +A22*y.^2 + ...
            K * exp(-x.^2/(2*s1^2) - y.^2/(2*s2^2)), ...
            'independent', {'x', 'y'}, 'dependent', 'z' );        
        
        opt_int = fitoptions('Method','NonlinearLeastSquares',...
            'StartPoint',x0,...
            'TolFun',1e-6,...
            'TolX',1e-6,...
            'MaxIter',1000,...
            'Display','off');
        
        % Fit the model
        fitted = fit([Xf Yf], Z, fun, opt_int);
        
        % Reconstruct radially
        for n=1:n_angle
            Z_fit(n,:) = fitted([X(n,:)' Y(n,:)'])';
        end
        
        fit_coeff = [fitted.A0 fitted.A11 fitted.A12 fitted.A21 fitted.A22 ...
            fitted.K fitted.s1 fitted.s2];

    case 'Breher'                
        fit_coeff = nan(n_bscan, 9);
        
        % Fitting configuration
        x0 = [0.3 -0.05 -0.02 0.3 0.04 -0.47 4.8 0.4 0.6];        
        lower = [0 -Inf -Inf -0.5 -0.5 -1 -Inf -Inf -Inf];
        upper = [1 Inf Inf    0.5  0.5  3 Inf Inf Inf];
        
        fun_int = fittype(@(a1, a2, a3, b1, b2, b3, c1, c2, c3,x) ...
            a1.*exp(-((x-b1)./c1).^2) + ...
            a2.*exp(-((x-b2)./c2).^2) + ...
            a3.*exp(-((x-b3)./c3).^2));
        
        opt_int = fitoptions('Method', 'NonlinearLeastSquares', ...
            'StartPoint', x0, ...
            'Lower', lower, ...
            'Upper', upper, ...
            'TolFun', 1e-6, ...
            'TolX',1e-6, ...
            'MaxIter',1000, ...
            'Display','off');
                
        for n=1:n_bscan            
            x = [-fliplr(rho(n+n_bscan,:)) rho(n, 2:end)];
            y = [fliplr(Z(n+n_bscan,:)) Z(n, 2:end)]; % Dont get center two times
            
            % Fit the model
            fitted = fit(x',y',fun_int,opt_int);
            
            % Reconstruct fitted curve
            y_fit = fitted(x)';
            
            % From fitted B-Scan to radial again
            Z_fit(n,:) = y_fit((end-n_point+1):end);% Right part of the B-Scan
            Z_fit(n+nBScans,:) = y_fit(n_point:-1:1);
            
            % Store Coefficients
            fit_coeff(n, :) = [fitted.a1 fitted.a2 fitted.a3 fitted.b1 ...
                fitted.b2 fitted.b3 fitted.c1 fitted.c2 fitted.c3];
        end
    case 'Dubis'
        error('Not supported yet');
    case 'Liu'
        error('Not supported yet');         
    case 'Scheibe'
        Z = 1e-3*Z;
        x0 = [1 0.5 1.5 0.1];
        lower = [0 0.1 0 0];
        upper = [40 10 40 10];

        fit_coeff = nan(n_angle, 4);

        fun = fittype(@(mu, sigma, gamma, alfa,x) ...
            mu * sigma^2 .* x.^gamma .* exp(-mu.*x.^gamma) + ...
            alfa .* (1 - exp(-mu.*x.^gamma)));
        
        opt_int = fitoptions('Method','NonlinearLeastSquares',...
            'StartPoint',x0,...
            'Lower',lower,...
            'Upper',upper,...
            'TolFun',1e-6,...
            'TolX',1e-6,...
            'MaxIter',1000,...
            'Display','off');

        % Get value at the centre
        pit_val = Z(1,1);

        % Fit the model in a radial fashion
        for i_angle = 1:n_angle
            
            x = rho(i_angle,:);
            y = Z(i_angle,:);
            
            % Normalize all points to that centre value setting the pit to 0
            y = y - pit_val;
            
            % Fit
            [fitted,~] = fit(x',y',fun, opt_int);
            fit_coeff(i_angle,:) = [fitted.mu fitted.sigma fitted.gamma fitted.alfa];
            
            % Reconstruct curve
            y_fit = fitted(x)';
            
            % Add pit thickness
            Z_fit(i_angle,:) = y_fit + pit_val;
        end
        Z_fit = 1e3*Z_fit;
    case 'Yadav'
        fit_coeff = nan(n_angle, 7);
        
        % Interior segment equation
%         fun_int = fittype('equation_yadav_int(x, alfa, beta, P0, P3)',...
%             'coefficients', {'alfa', 'beta'},...
%             'problem', {'P0', 'P3'});
        
        fun = @(alfa, beta, P0, P3,x) equation_yadav_int(x, alfa, beta, P0, P3);
        fun_int = fittype(fun,...
            'coefficients', {'alfa', 'beta'},...
            'problem', {'P0', 'P3'});
        
        % Exterior segment fitting configuration
%         fun_ext = fittype('equation_yadav_ext(x, alfa, P2x, P2y, P0, P3)',...
%             'coefficients', {'alfa', 'P2x', 'P2y'},...
%             'problem', {'P0', 'P3'});
        fun = @(alfa, beta, P2x, P2y, P0, P3, x) equation_yadav_ext(x, alfa, beta, P2x, P2y, P0, P3);
        fun_ext = fittype(fun, ...
            'coefficients', {'alfa', 'P2x', 'P2y'},...
            'problem', {'P0', 'P3'});

        for i_angle=1:n_angle
            
            x = rho(i_angle,:);
            y = Z(i_angle,:);
            
            % Get Max
            [~,ind_max] = max(y);
            ind_max = ind_max(1);
            
            % Fit interior
            x_int = x(1:ind_max);
            y_int = y(1:ind_max);
            
            P0 = [x_int(1) ;y_int(1)];
            P3 = [x_int(end) ;y_int(end)];
            
            % initialize alfa,beta to half width
            alfa0 = (P3(1) - P0(1))/2;
            beta0 = (P3(1) - P0(1))/2;
            x0 = [alfa0 beta0];
            lower = [0 0];
            upper = [2 2];
            
            opt_int = fitoptions('Method','NonlinearLeastSquares',...
                'StartPoint',x0,...
                'Lower',lower,...
                'Upper',upper,...
                'TolFun',1e-6,...
                'TolX',1e-6,...
                'MaxIter',1000,...
                'Display','off');
            
            fitted = fit(x_int', y_int', fun_int, opt_int, 'problem', {P0, P3});
            y_int_fit = fitted(x_int)';
            
            % Fit exterior
            x_ext = x(ind_max:end);
            y_ext = y(ind_max:end);
            
            P0 = [x_ext(1) ;y_ext(1)];
            P3 = [x_ext(end) ;y_ext(end)];
            
            alfa0 = (P3(1)-P0(1))/2; % half width of the exterior segment
            P2x0 = P3(1)-(P3(1)-P0(1))/2; % P3 - half width
            P2y0 = (P0(2)+P3(2))/2; % in the middle between P0y and P3y
            x0 = [alfa0  P2x0 P2y0];
            lower = [0 0 0.200];
            upper = [2 3 0.400];
            
            opt_ext = fitoptions('Method','NonlinearLeastSquares',...
                'StartPoint',x0,...
                'Lower',lower,...
                'Upper',upper,...
                'TolFun',1e-6,...
                'TolX',1e-6,...
                'MaxIter',1000,...
                'Display','off');
            
            fitted = fit(x_ext',y_ext', fun_ext, opt_ext, 'problem', {P0, P3});
            y_ext_fit = fitted(x_ext)';
            
            % Force first point to have the highest thickness
            y_ext_fit(y_ext_fit>y_ext_fit(1)) = y_ext_fit(1);
            
            % Recombine segments            
            Z_fit(i_angle, :) = [y_int_fit(1:end-1) y_ext_fit];            
        end
    case 'none'
        Z_fit = Z;
        fit_coeff = nan;
    case 'otherwise'
        error('Wrong model name');
end
end

function y = equation_yadav_int(x, alfa, beta, P0, P3)
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

t = linspace(0,1,length(x));

T = [1;0];

P1 = P0 + alfa*T;
P2 = P3 - beta*T;

Q = P0*Be(t,0) + P1*Be(t,1) + P2*Be(t,2) + P3*Be(t,3);

y = interp1(Q(1,:),Q(2,:),x);
end

function y = equation_yadav_ext(x, alfa, P2x, P2y, P0, P3)
% equationYadavExt - evaluate the model presented by Yadav et al. (2017)
% to model de external part of the B-Scan (beyond rim)
%
%
% Input arguments:
%   alfa,P2x,P2Y,P0,P3: model coefficients (defining Bezier curves)
%   x: evaluating point

t = linspace(0,1,length(x));

T = [1;0];

P1 = P0 + alfa*T;
P2 = [P2x;P2y];

Q = P0*Be(t,0) + P1*Be(t,1) + P2*Be(t,2) + P3*Be(t,3);
y = interp1(Q(1,:),Q(2,:),x);
end

function B = Be(t,i)
% Be - Evaluate Bernstein polynomial
%
% B = Be(t,i)
%
% Input arguments:
%   t: evaluation point (in range [0,1])
%   i: polynomial order
%
% Output arguments:
%   B: Evaluated polynomial value

switch i
    case 0
        B = (1-t).^3;
    case 1
        B = 3.*t.*(1-t).^2;
    case 2
        B = 3.*t.^2.*(1-t);
    case 3
        B = t.^3;
end
end