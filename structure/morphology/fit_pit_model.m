function [Z_fit, Fit_coeff] = fit_pit_model(theta, rho, Z, pit_model, varargin)
%   [Z_fit, Fit_coeff] = fit_pit_model(theta, rho, Z, pit_model)
%   Detail explanation goes here
%
%   Input arguments (required):
%  
%   'theta'          Matrix with theta coordinates (polar). It expects data in
%                    format n_directions x n_points/direction.
%
%   'rho'            Matrix with rho value of polar coordinates.It expects data in
%                    format n_directions x n_points/direction.
%
%   'Z'              Thickness map.
%
%   'pit_model'      String defining the mathematical model to be used.
%                    Options:['Breher','Ding','Dubis','Liu','Scheibe','Yadav']
%  
%   Input arguments (name-value pairs):
%
%   'max_iter'       Maximum number of iterations for each fitting.
%                    Default: 1000
%
%   'tol_x'          Tolerance of the errors during fitting.
%                    Default: 1e-6
%
%   'tol_fun'        Tolerance of the function during fitting.
%                    Default: 1e-6
%
%   'x0'             Initial iteration coefficient values. By default 
%                    specific values related to the model are used.
%   
%   'x_low'          Inferior limits of model coefficients
%   
%   'x_sup'          Superior limits of model coefficients
%
%
%   Output arguments:
%  
%   'Z_fit'          Matrix with fitted values
%
%   'Fit_coeff'      Struct with estimated model coefficients
%  
%
%   
%   Notes
%   -----
%   Yadav model does not work yet.
%
%
%   References
%   ----------
%   [1] Romero-Bascones et al., Foveal Pit Morphology Characterization: A 
%   Quantitative Analysis of the Key Methodological Steps, Entropy, 2021
%   https://doi.org/10.3390/e23060699
%
%   [1] Breher K. et al., Direct Modeling of Foveal Pit Morphology from 
%   Distortion-Corrected OCT Images, Biomedical Optics Express, 2019.
%
%   [2] Ding Y. et al., Application of an OCT Data-Based Mathematical Model of  
%   the Foveal Pit in Parkinson Disease, Journal fo Neural Transmission, 2014.
%
%   [3] Dubis A.M. et al., "Reconstructing Foveal Pit Morphology from Optical 
%   Coherence Tomography Imaging", British Journal of Ophthalmology, 2009.
%
%   [4] Liu L. et al., Sloped Piecemeal Gaussian Model for Characterising 
%   Foveal Pit Shape, Ophthalmic Physiological Optics, 2016.
%
%   [5] Scheibe P. et al., Parametric Model for the 3D Reconstruction of
%   Individual Fovea Shape from OCT Data, Experimental Eye Research, 2014.
%
%   [6] Yadav S. et al., CuBe: Parametric Modeling of 3D Foveal Shape Using 
%   Cubic Bézier, Biomedical Optics Express, 2017.
%
%
%   Example 1
%   ---------      
%   % Example description
%
%     I = [1 1 5 6 8 8;2 3 5 7 0 2; 0 2 3 5 6 7];
%     [GLCMS,SI] = graycomatrix(I,'NumLevels',9,'G',[])
%     
%
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2021


% Parse inputs
p = inputParser;
addRequired(p, 'theta', @(x)validateattributes(x,{'double'}, {'nonnan','nonempty'}));
addRequired(p, 'rho', @(x)validateattributes(x,{'double'}, {'nonnan', 'nonempty','nonnegative'}));
addRequired(p, 'Z', @(x)validateattributes(x,{'double'}, {'nonempty','positive'}));
addRequired(p, 'pit_model', @(x)validateattributes(x,{'char'}, {'scalartext'}));
addParameter(p, 'max_iter', 1000, @(x)validateattributes(x,{'double'}, {'integer','scalar'}));
addParameter(p, 'tol_x', 1e-6, @(x)validateattributes(x,{'double'}, {'positive','real','scalar'}));
addParameter(p, 'tol_fun', 1e-6, @(x)validateattributes(x,{'double'}, {'positive','real','scalar'}));
parse(p,theta, rho, Z, pit_model, varargin{:})
max_iter = p.Results.max_iter;
tol_x = p.Results.tol_x;
tol_fun = p.Results.tol_fun;

% Check the presence of nan values
if sum(isnan(Z(:))) > 0
    warning('NaN values in layer. Fitting might result in errors.');
end

[n_angle, n_point] = size(Z);
n_bscan = n_angle/2;

Z_fit = nan(n_angle, n_point);

switch pit_model   
    case 'Breher'                
        params = {'a1', 'a2', 'a3', 'b1', 'b2', 'b3', 'c1', 'c2', 'c3'};
        fit_type = 'bscan';
        fit_unit = 'mm';                
        fun = fittype(@(a1, a2, a3, b1, b2, b3, c1, c2, c3,x) ...
            a1.*exp(-((x-b1)./c1).^2) + ...
            a2.*exp(-((x-b2)./c2).^2) + ...
            a3.*exp(-((x-b3)./c3).^2));
                        
        x_low = [0   -Inf  -Inf  -0.5 -0.5   -1   -Inf -Inf -Inf];
        x0 =    [0.3 -0.05 -0.02  0.3  0.04 -0.47  4.8  0.4  0.6];        
        x_up =  [1    Inf   Inf   0.5  0.5    3    Inf  Inf  Inf];
        
   case 'Ding'   
        params = {'A0', 'A11', 'A12', 'A21', 'A22', 'K', 's1', 's2'};
        fit_type = '3d';
        fit_unit = 'mm';                        
        fun = fittype( @(A0, A11, A12, A21, A22, K, s1, s2, x, y) ...
            A0 + A11*x + A21*y +A12*x.^2 +A22*y.^2 + ...
            K * exp(-x.^2/(2*s1^2) - y.^2/(2*s2^2)), ...
            'independent', {'x', 'y'}, 'dependent', 'z' );         
       
        x_low = -inf(1, length(params));
        x0 = [360 8 -10 1.5 -8 -130 500 500]/1000;
        x_up = inf(1, length(params));
        
    case 'Dubis'
        params = {'a1', 'a2', 'b1', 'b2', 'c1', 'c2'};
        fit_type = 'bscan';
        fit_unit = 'mm';                        
        fun = fittype(@(a1,a2,b1,b2,c1,c2,z,x) ...
            a1.*exp((x-b1).^2./(-2*c1)) - a2.*exp((x-b2).^2./(-2*c2)) + z);       
        x_low = [ 0   0   -0.3 -0.3 0   0    0];
        x0 =    [0.2 0.14   0    0  5  0.09 0.8];
        x_up =  [ 1  0.3   0.3  0.1 30  10   2];
              
    case 'Liu'
        params = {'a', 'f', 'g', 'lambda', 'mu', 'sigma'};
        fit_type = 'bscan';
        fit_unit = 'um';                        
        myfun = @(a, f, g, lambda, mu, sigma, x)  equation_liu(x,a,f,g,lambda,mu,sigma);
        fun = fittype(myfun);
        x_low = [10 -10 200 0 -1 0.1];
        x0 =    [80   0 320 0  0 0.5];
        x_up =  [200 10 450 5  1 1.0];                

    case 'Scheibe'
        params = {'mu', 'sigma', 'gamma', 'alpha'};
        fit_type = 'radial';
        fit_unit = 'mm';
        fun = fittype(@(mu, sigma, gamma, alpha,x) ...
            mu * sigma^2 .* x.^gamma .* exp(-mu.*x.^gamma) + ...
            alpha .* (1 - exp(-mu.*x.^gamma)));
        
        x_low = [0  0.1  0    0  ];
        x0 =    [1  0.5  1.5  0.1];
        x_up =  [40  10  40   10 ];        
        
    case 'Yadav'
        [Z, unit_in] = convert_mm_um(Z, 'mm');
        
        % Interior/exterior segment equations        
        fun = @(alpha, beta, P0, P3,x) equation_yadav_int(x, alpha, beta, P0, P3);
        fun_int = fittype(fun, 'coefficients', {'alpha', 'beta'},...
                               'problem', {'P0', 'P3'});
        
        fun = @(alpha, P2x, P2y, P0, P3, x) equation_yadav_ext(x, alpha, P2x, P2y, P0, P3);
        fun_ext = fittype(fun, 'coefficients', {'alpha', 'P2x', 'P2y'},...
                               'problem', {'P0', 'P3'});

        for i_angle=1:n_angle
            
            x = rho(i_angle,:);
            y = Z(i_angle,:);
            
            % Get Max
            [~, ind_max] = max(y);
            ind_max = ind_max(1);
            
            % Fit interior
            x_int = x(1:ind_max);
            y_int = y(1:ind_max);
            
            P0 = [x_int(1) ;y_int(1)];
            P3 = [x_int(end) ;y_int(end)];
            
            % initialize alpha,beta to half width
            alpha0 = (P3(1) - P0(1))/2;
            beta0 = (P3(1) - P0(1))/2;
            x0 = [alpha0 beta0];
            x_low = [0 0];
            x_up = [2 2];
            
            opt = fitoptions('Method','NonlinearLeastSquares',...
                             'StartPoint',x0,...
                             'Lower',x_low,...
                             'Upper',x_up,...
                             'TolFun',tol_fun,...
                             'TolX',tol_x,...
                             'MaxIter',max_iter,...
                             'Display','off');
             
            fitted = fit(x_int', y_int', fun_int, opt, 'problem', {P0, P3});
            y_int_fit = fitted(x_int)';
            
            Fit_coeff.alpha_int(i_angle) = fitted.alpha;
            Fit_coeff.beta(i_angle) = fitted.beta;
            
            % Fit exterior
            x_ext = x(ind_max:end);
            y_ext = y(ind_max:end);
            
            P0 = [x_ext(1) ;y_ext(1)];
            P3 = [x_ext(end) ;y_ext(end)];
            
            alpha0 = (P3(1)-P0(1))/2; % half width of the exterior segment
            P2x0 = P3(1)-(P3(1)-P0(1))/2; % P3 - half width
            P2y0 = (P0(2)+P3(2))/2; % in the middle between P0y and P3y
            x0 = [alpha0  P2x0 P2y0];
            x_low = [0 0 0.200];
            x_up = [2 3 0.400];
            
            opt_ext = fitoptions('Method','NonlinearLeastSquares',...
                                 'StartPoint',x0,...
                                 'Lower',x_low,...
                                 'Upper',x_up,...
                                 'TolFun',tol_fun,...
                                 'TolX',tol_x,...
                                 'MaxIter',max_iter,...
                                 'Display','off');

            fitted = fit(x_ext',y_ext', fun_ext, opt_ext, 'problem', {P0, P3});
            y_ext_fit = fitted(x_ext)';
            
            Fit_coeff.alpha_ext(i_angle) = fitted.alpha;
            Fit_coeff.P2x(i_angle) = fitted.P2x;
            Fit_coeff.P2y(i_angle) = fitted.P2y;

            % Force first point to have the highest thickness
            y_ext_fit(y_ext_fit > y_ext_fit(1)) = y_ext_fit(1);
            
            % Recombine segments            
            Z_fit(i_angle, :) = [y_int_fit(1:end-1) y_ext_fit];            
        end
        Z_fit = convert_mm_um(Z_fit, unit_in);

        return
    case 'none'
        Z_fit = Z;
        Fit_coeff = nan;
        return
    otherwise
        error('Unknown model name');
end


opt = fitoptions('Method','NonlinearLeastSquares',...
                 'StartPoint',x0,...
                 'Lower',x_low,...
                 'Upper',x_up,...
                 'TolFun',tol_fun,...
                 'TolX',tol_x,...
                 'MaxIter',max_iter,...
                 'Display','off');

[Z, unit_in] = convert_mm_um(Z, fit_unit);

% Fit the model
switch fit_type
    case 'radial'
       
        % Set pit value to 0
        pit_val = Z(1,1);
        Z = Z - pit_val;
        
        % Fit the model in a radial fashion
        for i_angle = 1:n_angle
            
            x = rho(i_angle,:);
            y = Z(i_angle,:);
            [fitted,~] = fit(x, y, fun, opt);
            
            y_fit = fitted(x)';  %  Reconstruct curve
            
            Z_fit(i_angle,:) = y_fit + pit_val;  %  Add pit thickness
            
            for i=1:length(params)
                Fit_coeff.(params{i})(i_angle) = fitted.(params{i});
            end
        end
        
    case 'bscan'                
        for n=1:n_bscan            
            x = [-fliplr(rho(n+n_bscan,:)) rho(n, 2:end)];
            y = [fliplr(Z(n+n_bscan,:)) Z(n, 2:end)]; % Dont get center two times
            
            % Fit the model
            fitted = fit(x',y',fun,opt);
            
            % Reconstruct fitted curve
            y_fit = fitted(x)';
            
            % From fitted B-Scan to radial again
            Z_fit(n,:) = y_fit((end-n_point+1):end);% Right part of the B-Scan
            Z_fit(n+n_bscan,:) = y_fit(n_point:-1:1);
            
            % Store Coefficients
            for i=1:length(params)
                Fit_coeff.(params{i})(n,:) = fitted.(params{i});
            end

        end     
        
    case '3d'       
        [X, Y] = pol2cart(theta, rho);
        
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
   
        % Fit the model
        fitted = fit([Xf Yf], Z, fun, opt);
        
        % Reconstruct radially
        for n=1:n_angle
            Z_fit(n,:) = fitted([X(n,:)' Y(n,:)'])';
        end
        
        for i=1:length(params)
            Fit_coeff.(params{i}) = fitted.(params{i});
        end
    otherwise
        error("Unknown fit type");
end

Z_fit = convert_mm_um(Z_fit, unit_in);

end

function y = equation_yadav_int(x, alpha, beta, P0, P3)
% equationYadavInt - evaluate the model presented by Yadav et al. (2017)
% to model de inner part of the B-Scan (foveal center to rim)
%
% yf = equationYadavInt(x,alpha,beta,P0,P3)
%
% Input arguments:
%   alpha,beta,P0,P3: model coefficients (defining Bezier curves)
%   x: evaluating point
%
% Output arguments:
%   yf: evaluated value

t = linspace(0,1,length(x));

T = [1;0];

P1 = P0 + alpha*T;
P2 = P3 - beta*T;

Q = P0*Be(t,0) + P1*Be(t,1) + P2*Be(t,2) + P3*Be(t,3);

y = interp1(Q(1,:),Q(2,:),x);
end

function y = equation_yadav_ext(x, alpha, P2x, P2y, P0, P3)
% equationYadavExt - evaluate the model presented by Yadav et al. (2017)
% to model de external part of the B-Scan (beyond rim)
%
%
% Input arguments:
%   alpha,P2x,P2Y,P0,P3: model coefficients (defining Bezier curves)
%   x: evaluating point

t = linspace(0,1,length(x));

T = [1;0];

P1 = P0 + alpha*T;
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

function y = equation_liu(x,a,f,g,lambda,mu,sigma)
%
% y = equationLiu(x,a,f,g,lambda,mu,sigma)
%
% Input arguments:
%   a,f,g,lambda,mu,sigma: model coefficients
%   x: evaluating point
%
% Output arguments:
%   y: evaluated value

G = exp(-(x-mu).^2/sigma^2);
centre = (x>=(mu-lambda/2) & x<=(mu+lambda/2));

G1 = G;
G1(centre) = max(G);
G2 = g - (a*G1+f*x); 

y = G2;
end
