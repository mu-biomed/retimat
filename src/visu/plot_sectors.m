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
%   'alpha'          Transparency of patches. Range: [0,1]. 
%                    Default: 1 
%
%   'edge_color'     Color of the edges of each patch.
%                    Default: 'k' (black)
%
%   'axis_off'       Set it to true to hide axis.
%                    Default: 'true'
%
%   'axis_equal'     Set it to true to set a aspect ratio of 1.
%                    Default: 'true'
%
%
%   References
%   ----------
%   [1] Eric (2020). donut (https://www.github.com/ericspivey/donut), GitHub.
%   Retrieved January 29, 2020. (Used to create the plot_patch function).
%
%
%   Example
%   ---------      
%   % Plot ETDRS sectorization
%
%   [header, seg, ~, ~] = read_vol(file,'verbose', 'coordinates');
%   Thickness = compute_thickness(seg, 'TRT', header.scale_z);
%   [TRT, Sectors] = sectorize_map(X, Y, TRT, 'mean', 'etdrs');
%   plot_sectors(TRT, Sectors);
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2021

parsed = parse_inputs(varargin);

n_point    = parsed.n_point;
alpha      = parsed.alpha;
edge_color = parsed.edge_color;

switch Sectors.type    
    case 'regular'
        % Compute coordinates of each square center
        X_center = (Sectors.X_edge(1:end-1) + Sectors.X_edge(2:end))/2;
        Y_center = (Sectors.Y_edge(1:end-1) + Sectors.Y_edge(2:end))/2;        
        
        I = reshape(Z, [Sectors.n_y Sectors.n_x]);
        im = imagesc(X_center, Y_center, I);
        im.AlphaData = alpha;
        xlim([Sectors.X_edge(1) Sectors.X_edge(end)]);
        ylim([Sectors.Y_edge(1) Sectors.Y_edge(end)]);        
    
    case 'ring'
        n_sect = Sectors.n_sect;
        radius = Sectors.radius;
        for n=1:n_sect
            plot_ring(radius(n), radius(n+1), Z(n), n_point, alpha, edge_color);
        end
        
    case 'pie'
        n_sect  = Sectors.n_sect;
        radius  = Sectors.radius;
        theta_0 = Sectors.theta_0;
        
        theta_edge = linspace(theta_0, theta_0+2*pi, n_sect+1);
        for i_sect=1:n_sect
            plot_pie(theta_edge(i_sect), theta_edge(i_sect+1), radius,...
                Z(i_sect), n_point, alpha, edge_color);                                     
        end
        
        xlim(radius(end) * [-1 1]);
        ylim(radius(end) * [-1 1]);
        
    case 'wedge'
        radius  = Sectors.radius;
        n_angle = Sectors.n_angle;
        
        theta = linspace(0, 2*pi, n_angle+1) + Sectors.theta_0;
        
        % center circle
        plot_circle(radius(1), Z(1), n_point, alpha, edge_color);
        
        % plot wedges looping through rings and angles
        i = 2;
        for n=1:length(radius)-1
            for t=1:n_angle
                plot_wedge(theta(t), theta(t+1), radius(n), radius(n+1), ...
                    Z(i), n_point, alpha, edge_color);
                i = i+1;
            end
        end

        xlim(max(radius) * [-1 1]);
        ylim(max(radius) * [-1 1]);
        
    otherwise
        error("Unknown sectorization type. Accepted values: ['regular','ring','pie','wedge','etdrs']");
end

if parsed.axis_equal
    daspect([1 1 1]);
    % axis equal % does not work well
end

if parsed.axis_off
    axis off;
end

function parsed = parse_inputs(args)

parsed.n_point    = 50;
parsed.axis_off   = true;
parsed.axis_equal = true;
parsed.alpha      = 1;
parsed.edge_color = 'k';

n_arg = floor(length(args)/2);

for i_arg=1:n_arg
    switch args{2*i_arg - 1}
        case 'n_point'
            parsed.n_point = args{2*i_arg};
        case 'alpha'
            parsed.alpha = args{2*i_arg};            
        case 'axis_off'
            parsed.axis_off = args{2*i_arg};                        
        case 'axis_equal'
            parsed.axis_equal = args{2*i_arg};
        case 'edge_color'
            parsed.edge_color = args{2*i_arg};            
        otherwise
            warning('Unknown parameter');
    end
end

function plot_ring(rad_int, rad_ext, z, n_point, alpha, edge_color)

theta = linspace(0, 2*pi, n_point);
rho0 = rad_int * ones(1,n_point);
rho1 = rad_ext * ones(1,n_point);

[x0, y0] = pol2cart(theta, rho0);
[x1, y1] = pol2cart(theta, rho1);

patch([x0 x1], [y0 y1], z, 'linestyle', 'none', 'FaceAlpha', alpha);

line(x0, y0, 'color', edge_color);
line(x1, y1, 'color', edge_color);

function plot_circle(radius, z, n_point, alpha, edge_color)

theta = linspace(0, 2*pi, n_point);
rho = radius * ones(1, n_point);
[x,y] = pol2cart(theta, rho);
patch(x, y, z, 'linestyle', 'none', 'FaceAlpha', alpha);
line(x, y, 'Color', edge_color);

function plot_pie(th0, th1, radius, z, n_point, alpha, edge_color)

% A pie is made of 3 segments

% Straight line (inner to outer)
a1 = [th0 th0];
r1 = [0 radius];

% Outer arc
a2 = linspace(th0, th1, n_point);
r2 = radius * ones(1, n_point);

% Straight line (outer to inner)
a3 = [th1 th1];
r3 = [radius 0];

[x,y]= pol2cart([a1,a2,a3],[r1,r2,r3]);
patch(x, y, z, 'EdgeColor', edge_color, 'FaceAlpha', alpha); 

function plot_wedge(th0, th1, rh0, rh1, z, n_point, alpha, edge_color)

% A wedge is made of four segments

% Straight line (inner to outer)
a1 = [th0 th0];
r1 = [rh0 rh1];

% Outer arc
a2 = linspace(th0, th1, n_point);
r2 = rh1 * ones(1, n_point);

% Straight line (outer to inner)
a3 = [th1 th1];
r3 = [rh1 rh0];

% Inner arc
a4 = linspace(th1, th0, n_point);
r4 = rh0 * ones(1, n_point);

[x,y]= pol2cart([a1,a2,a3,a4], [r1,r2,r3,r4]);
patch(x, y, z, 'EdgeColor', edge_color, 'FaceAlpha', alpha); 
