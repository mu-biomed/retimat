function correct_display_distortion(I, header, bm, ilm, visu, correct_thickness)

% Get scan characteristics
[n_axial, n_ascan] = size(I);

sy = header.scale_z;
sx = header.scale_x;
fov = header.fov;

% Define nodal distance
d_nodal = 22;

% Field of view and fanning angle of each A-Scan
theta = linspace(-fov/2,fov/2,n_ascan); 

% Define real-world coordinates of B-Scan and layers
% B-Scan coordinates
max_x = sx*(n_ascan-1)/2;
min_x = -max_x;

min_y = -sy*(n_axial-1)/2;
max_y = sy*(n_axial-1)/2;

x_range = min_x:sx:max_x;
y_range = min_y:sy:max_y;

[X, Y] = meshgrid(x_range,-y_range);

% Compute real world segmentation coordinates
y0 = Y(1,1);
y1 = Y(end,1);

bm_c = y0 + (y1-y0)/(n_axial-1)*bm;
ilm_c = y0 + (y1-y0)/(n_axial-1)*ilm;

% Get the point at central bm point
ind = [floor(n_ascan/2) ceil(n_ascan/2)];
y_fovea = mean(bm_c(ind));

% 1D translate all y points (set origin at the fixation point)
Y = Y - y_fovea;
ilm_c = ilm_c - y_fovea;
bm_c = bm_c - y_fovea;

% Corrected A-scan coordinates
X_a = nan(n_axial,n_ascan);
Y_a = nan(n_axial,n_ascan);

for a=1:n_ascan    
    % Distance from the fovea to the first pixel in depth of central A-Scan
    d_start = max(Y(:));
    
    % Distance from A-Scan start to nodal point
    d_ilm = d_nodal - d_start;
    
    % Compute rho coordinates (centred at nodal point)
    rho = d_ilm + (0:sy:sy*(n_axial-1));
    
    % Angle of A-Scan in polar coordinates
    phi_ascan = (theta(a) - 90)*pi/180;

    % Convert coordinates to cartesian
    [x_ascan,y_ascan] = pol2cart(repmat(phi_ascan,1,n_axial),rho);
    
    % Translate them to the nodal point
    y_ascan = y_ascan + d_nodal;
    
    % Store a-scan coordinates
    X_a(:,a) = x_ascan;
    Y_a(:,a) = y_ascan;    
end

% Correct ILM and BM segmentation
x_ilm = nan(1,n_ascan);
y_ilm = nan(1,n_ascan);
x_bm = nan(1,n_ascan);
y_bm = nan(1,n_ascan);

for a=1:n_ascan    
    % Distance from A-Scan start to nodal point
    d_ilm = d_nodal - ilm_c(a);
    d_bm = d_nodal - bm_c(a);   
    
    % Angle of A-Scan in polar coordinates
    phi_ascan = (theta(a) - 90)*pi/180;
        
    % Convert coordinates to cartesian
    [x_ilm(a),y_ilm(a)] = pol2cart(phi_ascan,d_ilm);    
    [x_bm(a),y_bm(a)] = pol2cart(phi_ascan,d_bm);        
end

% Translate them to the nodal point
y_ilm = y_ilm + d_nodal;
y_bm = y_bm + d_nodal;

% Fit eye over corrected RPE
[Rc, x_cen, y_cen, tilt] = get_retinal_curvature(x_bm, y_bm);

if visu
    % Compute wavefront
    N = 50;
    [x_wave,y_wave] = pol2cart(linspace(-3*pi/8,-5*pi/8,N),d_nodal*ones(1,N));
    y_wave = y_wave + d_nodal;

    y_eye = y_cen - sqrt(Rc^2 - (x_wave - x_cen).^2);

    % Visualize B-Scan with nodal point and A-Scans
    subplot(121);hold on;
    surf(X,Y,bscan,'EdgeColor','none');
    view(0,90);colormap(gray);
    xlabel('x');ylabel('z');
    plot3(x_range,bm_c,ones(1,n_ascan),'r');
    plot3(x_range,ilm_c,ones(1,n_ascan),'g')
    scatter3(mean(x_range(ind)),0,1,'yellow','filled')
    % scatter3(mean(x_range(ind)),d_nodal,1,'blue','filled');
    
    % Visualize B-Scan with nodal point and A-Scans
    subplot(122);hold on;
    % B-Scan
    surf(X_a,Y_a,bscan,'EdgeColor',"none",'HandleVisibility','off');
    % Wavefront
    plot3(x_wave,y_wave,ones(1,N),'--r');
    % Segmentation
    plot3(x_ilm,y_ilm,ones(1,n_ascan),'g');
    plot3(x_bm,y_bm,ones(1,n_ascan),'b');
    % Fitted Eye
    plot3(x_wave,y_eye,ones(1,N),'m');
    % Plot options
    xlabel('x (mm)');ylabel('y(mm)');colormap(gray);view(0,90);
    legend('wavefront','ilm','bm','fitted eye');
    set(gca,'FontSize',14);

end

if correct_thickness
%Correct thickness perpendicularity
f = figure(2);clf;hold on;
set(f,'Position', [0 0 1400 500]);
i_ascan = 288;

x0 = x_bm(i_ascan);
y0 = y_bm(i_ascan);
x1 = x_bm(i_ascan + 1);
y1 = y_bm(i_ascan + 1);

% Tangent line
m = (y1 - y0)/(x1 - x0);  %  slope
x = -4:4;
y = m*(x-x0) + y0;

% Normal line 
mn = -1/m;  %  slope
yn = y0:0.01:max(y_ilm);  %  search range
xn = (yn - y0)/mn + x0;

% Intersect
D = pdist2([xn' yn'],[x_ilm' y_ilm']);
[~,ind] = min(D(:));
[a,b] = ind2sub(size(D),ind);

x2 = x_ilm(b);
y2 = y_ilm(b);
% Visualization
surf(X_a,Y_a,bscan,'EdgeColor',"none",'HandleVisibility','off');  %  B-Scan
plot3(x_wave,y_wave,ones(1,N),'--r');  %  wavefront
plot3(x_ilm,y_ilm,ones(1,n_ascan),'g');  %  ilm
plot3(x_bm,y_bm,ones(1,n_ascan),'b');  %  bm
% plot3(x_wave,curve(x_wave),ones(1,N),'m');  % fitted eye

plot3(x,y,ones(1,length(x)),'--y');
plot3(xn,yn,ones(1,length(xn)),'--c');
scatter3(x0,y1,1,'r','filled','HandleVisibility','off');
scatter3(x2,y2,1,'r','filled','HandleVisibility','off');

xlabel('x (mm)');ylabel('y(mm)');colormap(gray);view(0,90);
legend('wavefront','ilm','bm','tangent','normal','Location','southeast');
set(gca,'FontSize',14);xlim([-3.5 3.5]);ylim([-1 1]);


% Original
TRT_acq = 1e3 * sy * (bm - ilm);
[~,pit] = min(TRT_acq);
x_acq = x_range - x_range(pit);

% Corrected
% TRT_cor = 1e3 * sqrt((x_ilm - x_bm).^2+(y_ilm-y_bm).^2);

% Corrected and perpendicular
TRT_per = nan(n_ascan, 1);
x_per = nan(n_ascan, 1);
for i_ascan=1:(n_ascan-1)
    
   
    x0 = x_bm(i_ascan);
    y0 = y_bm(i_ascan);
    x1 = x_bm(i_ascan + 1);
    y1 = y_bm(i_ascan + 1);
    
    % Tangent line
    m = (y1 - y0)/(x1 - x0);  %  slope
    x = -4:4;
    y = m*(x-x0) + y0;
    
    % Normal line 
    mn = -1/m;  %  slope
    yn = y0:0.001:max(y_ilm);  %  search range
    xn = (yn - y0)/mn + x0;
    
    % Intersect
    D = pdist2([xn' yn'],[x_ilm' y_ilm']);
    [~,ind] = min(D(:));
    [a,b] = ind2sub(size(D),ind);
    
    x2 = x_ilm(b);
    y2 = y_ilm(b);
    
    TRT_per(i_ascan) = 1e3*sqrt((x2-x0).^2 + (y2-y0).^2);
    if i_ascan ==1
        x_per(i_ascan) = x0;
    else
        step = sqrt((x0 - x1).^2 + (y0-y1).^2);
        x_per(i_ascan) = x_per(i_ascan-1) + step;
    end
end
[~,pit] = min(TRT_per);
x_per = x_per - x_per(pit);
mask = ~isnan(x_per);
TRT_per = interp1(x_per(mask),TRT_per(mask),x_acq);

% Visualization
f = figure(3);clf;
subplot(211);hold on;
plot(x_acq,TRT_acq,'r','LineWidth',1.5);
plot(x_acq,TRT_per,'b','LineWidth',1.5);
grid on;xlabel('x');ylabel('\mum'); xlim([-3 3]);%ylim([200 400]);
legend({'regular','correct'},'Location','southeast')
subplot(212);hold on;
plot(x_acq, TRT_acq - TRT_per,'k','LineWidth',1.5);grid on;xlim([-3 3]);
legend('regular-correct')
end
