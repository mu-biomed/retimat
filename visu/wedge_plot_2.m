function wedge_plot_2(Grid,Data)

% center circle
polsect(0,2*pi,0,Grid.rho(1),Data(1));

% Plot wedges looping through rings and angles
i = 2;
for n=1:length(Grid.rho)-1
    for t=1:Grid.n_angle
        polsect(Grid.theta(t)+Grid.theta0,Grid.theta(t+1)+Grid.theta0,...
            Grid.rho(n),Grid.rho(n+1),Data(i));
        i = i+1;
    end
end

xlim([-1 1]*max(abs(Grid.rho(:))));
ylim([-1 1]*max(abs(Grid.rho(:))));


function polsect(th0,th1,rh0,rh1,cl)
% This function creates a patch from polar coordinates
% Cite As
% Eric (2020). donut (https://www.github.com/ericspivey/donut), GitHub. Retrieved January 29, 2020.

N = 50; % al aumentarla sube resoución y tiempo de procesado
a1 = linspace(th0,th0,N);
r1 = linspace(rh0,rh1,N);
a2 = linspace(th0,th1,N);
r2 = linspace(rh1,rh1,N);
a3 = linspace(th1,th1,N);
r3 = linspace(rh1,rh0,N);
a4 = linspace(th1,th0,N);
r4 = linspace(rh0,rh0,N);
[X,Y]=pol2cart([a1,a2,a3,a4],[r1,r2,r3,r4]);
p=patch(X,Y,cl); % Note: patch function takes text or matrix color def
axis equal
axis off
