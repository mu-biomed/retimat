function plot_grid(conf)


% Set default params
circle_radius = conf.circle_radius;
grid_width = 0.1;
n_point = 50;
n_ring = 3;
n_angle = 24;
color_background = 'none'; 
color_grid_back = [220 220 220]./255;
color_grid = [151 172 191]./255;
color_text = [90 100 91]./255;%ColorText = [120 150 160]./255;
color_fill= [239 92 84]./255;

    % Adjust for user params
if isfield(conf,'circle_radius')
    circle_radius = conf.circle_radius;
end
if isfield(conf,'grid_width')
    grid_width = conf.grid_width;
end
if isfield(conf,'n_point')
    n_point = conf.n_point;
end
if isfield(conf,'n_angle')
    n_angle = conf.n_angle;
end

% Prepare grid
set(gca,'Color',color_background);
set(gca,'XColor','none');
set(gca,'YColor','none');
axis equal;

% Draw outer circle
theta = linspace(0,2*pi - 2*pi/n_angle,n_angle);
[Xr,Yr] = pol2cart(linspace(0,2*pi - 2*pi/n_point,n_point),repmat(circle_radius,1,n_point));
patch(Xr,Yr,color_grid_back,'EdgeColor',color_grid,'LineWidth',grid_width);hold on;
Ytext = circle_radius - circle_radius/n_ring*0.1;
Xtext = (circle_radius/n_ring*(n_ring+1) - circle_radius/n_ring*n_ring)*0.05;
text(Xtext,Ytext,num2str(round(circle_radius/n_ring*n_ring,3,'significant')),'Color',color_text,'FontWeight','Bold');

% Draw radial lines in each direction
[Xg,Yg] = pol2cart(theta,circle_radius*ones(1,n_angle));
for n=1:n_angle
    plot([0 Xg(n)],[0 Yg(n)],'Color',color_grid,'LineWidth',grid_width); 
end

% Draw Rings with numeric text
for n=1:n_ring-1
    [Xr,Yr] = pol2cart(linspace(0,2*pi - 2*pi/n_point,n_point),repmat(circle_radius/n_ring*n,1,n_point));
    plot([Xr Xr(1)],[Yr Yr(1)],'Color',color_grid,'LineWidth',grid_width);

    Ytext = circle_radius/n_ring*n - circle_radius/n_ring*0.1;
    Xtext = (circle_radius/n_ring*(n+1) - circle_radius/n_ring*n)*0.05;
    text(Xtext,Ytext,num2str(round(circle_radius/n_ring*n,3,'significant')),'Color',color_text,'FontWeight','Bold');
end

