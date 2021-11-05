function wedge_plot_pip(rho,theta,TRT,param)


% Get wedge delimitators
rho_wedge = linspace(param.RhoCenter,param.RhoMax,param.nRing+1);
theta_wedge = linspace(0,2*pi,param.nWedge+1); % define cut angles

% Set points per wedge side
N = 50; % al aumentarla sube

% Plot center circle
mask = rho <= rho_wedge(2); % mask of values
meanThick = mean(TRT(mask)); % compute value
[X,Y] = pol2cart(linspace(0,2*pi,N),rho_wedge(1)*ones(1,N));
patch(X,Y,meanThick); % Note: patch function takes text or matrix color def
        
% Loop throug wedges
for n=1:param.nRing
    for w=1:param.nWedge
        
        % Get mask for that wedge
        mask_rho = rho > rho_wedge(n) & rho <=rho_wedge(n+1); 
        mask_theta = theta >= theta_wedge(w) & theta < theta_wedge(w+1); 
        mask = mask_rho & mask_theta;
        
        % Apply mask to get mean thickness 
        meanThick = mean(TRT(mask));
        
        % Draw wedge        
        th0 = theta_wedge(w);
        th1 = theta_wedge(w+1);        
        rh0 = rho_wedge(n);
        rh1 = rho_wedge(n+1);
        
        a1 = linspace(th0,th0,N);
        r1 = linspace(rh0,rh1,N);
        a2 = linspace(th0,th1,N);
        r2 = linspace(rh1,rh1,N);
        a3 = linspace(th1,th1,N);
        r3 = linspace(rh1,rh0,N);
        a4 = linspace(th1,th0,N);
        r4 = linspace(rh0,rh0,N);
        
        [X,Y]=pol2cart([a1,a2,a3,a4],[r1,r2,r3,r4]);
        
        patch(X,Y,meanThick); % Note: patch function takes text or matrix color def
        axis equal
    end
end

% colormap(hot);

