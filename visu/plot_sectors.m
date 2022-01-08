function plot_sectors(Z, Sectors, varargin)
%PLOT_SECTORS Visualize sectorization results
%
%   plot_sectors(Z, Sectors, n_point)
%   Visualize sectorization results
%
%   Input arguments (mandatory):
%  
%   'Z'              Values of each sector.
%                    
%   'Sectors'        Struct defining the sectorization. Either manually defined
%                    or obtained after using sectorize_map().
%
%
%   Input arguments (name/value pair):
%
%   'n_point'        Number of points in each non-straight line segment. Higher
%                    number increases resolution and computation time.
%                    Default: 50
%
%   'alpha'          Transparency of patches [0,1]. 
%                    Default: 1
%
%   'edge_color'     Color of the edges of each patch.
%                    Default: 'k' (black)
%
%   
%   Notes
%   -----
%   
%
%
%   References
%   ----------
%   [1] Eric (2020). donut (https://www.github.com/ericspivey/donut), GitHub.
%   Retrieved January 29, 2020. (Used to create the plot_patch function).
%
%
%   Example 1
%   ---------      
%   % Example description
%
%   [header, seg, ~, ~] = read_vol(file,'verbose', 'coordinates');
%   Thickness = compute_thickness(seg, 'TRT', header.scale_z);
%   [X, Y, TRT] = resample_map(header.X_oct, header.Y_oct, Thickness.TRT, ...
%   'regular', 'n_point', 100, 'max_d', 2.5);
%   [TRT, Sectors] = sectorize_map(X, Y, TRT, 'mean', 'etdrs');
%
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2021


N = 50; 

switch Sectors.type    
    case 'regular'
        % Compute coordinates of each square center
        X_center = (Sectors.X_edge(1:end-1) + Sectors.X_edge(2:end))/2;
        Y_center = (Sectors.Y_edge(1:end-1) + Sectors.Y_edge(2:end))/2;        
        
        I = reshape(Z, [Sectors.n_y Sectors.n_x]);
        imagesc(X_center, Y_center, I);
        axis equal;
        xlim([X_center(1) X_center(end)]);
        ylim([Y_center(1) Y_center(end)]);        
    case 'ring'
    
    case 'pie'
        
    case 'wedge'
        radius = Sectors.radius;
        n_angle = Sectors.n_angle;
        theta = linspace(0, 2*pi, n_angle+1) + Sectors.theta_0;
        
        % center circle
        plot_circle(radius(1), Z(1), N);
        
        % plot wedges looping through rings and angles
        i = 2;
        for n=1:length(radius)-1
            for t=1:n_angle
                plot_wedge(theta(t), theta(t+1), radius(n), radius(n+1), Z(i), N);
                i = i+1;
            end
        end

        xlim(max(radius) * [-1 1]);
        ylim(max(radius) * [-1 1]);
        axis equal
        axis off
        
    otherwise
        error("Unknown sectorization type. Accepted values: ['regular','ring','pie','wedge','etdrs']");
end

end

function plot_ring(rad_int, rad_ext, z, N)

theta = linspace(0, 2*pi, N);
rho0 = rad_int * ones(1,N);
rho1 = rad_ext * ones(1,N);

[x0, y0] = pol2cart(theta, rho0);
[x1, y1] = pol2cart(theta, rho1);

patch([x0 x1], [y0 y1], z, 'linestyle', 'none');

line(x0, y0, 'color','k');
line(x1, y1, 'color','k');

end

function plot_circle(radius, z, N)

theta = linspace(0, 2*pi, N);
rho = radius * ones(1,N);
[x,y] = pol2cart(theta, rho);
patch(x, y, z, 'linestyle', 'none');
line(x, y, 'Color','black');
end

function plot_pie(th0, th1, radius, N)

% A pie is made of 3 segments

% Straight line (inner to outer)
a1 = [th0 th0];
r1 = [0 radius];

% Outer arc
a2 = linspace(th0, th1, N);
r2 = radius * ones(1,N);

% Straight line (outer to inner)
a3 = [th1 th1];
r3 = [radius 0];

[x,y]= pol2cart([a1,a2,a3],[r1,r2,r3]);
patch(x, y, z); 

end

function plot_wedge(th0, th1, rh0, rh1, z, N)

% A wedge is made of four segments

% Straight line (inner to outer)
a1 = [th0 th0];
r1 = [rh0 rh1];

% Outer arc
a2 = linspace(th0, th1, N);
r2 = rh1 * ones(1,N);

% Straight line (outer to inner)
a3 = [th1 th1];
r3 = [rh1 rh0];

% Inner arc
a4 = linspace(th1, th0, N);
r4 = rh0 * ones(1,N);

[x,y]= pol2cart([a1,a2,a3,a4],[r1,r2,r3,r4]);
patch(x, y, z); 
end
