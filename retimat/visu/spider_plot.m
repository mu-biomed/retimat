% function spiderPlot(Data1,Data2,varargin)
function spider_plot(Data,varargin)
% -------------------------
% Data: matrix with N (observations)x M (angles)
% param: (optional) structure with parameters for plotting configuration
% -- center: mean or median
% -- superior: interval superior value % percentile (e.g. 75)
% -- inferior: interval inferior value % percentile (e.g. 25)
% -- circleRadius: radius of the radar outer circle. Change to hide
% outliers
% -- plotDots: flag to overlay raw data points or not 
%--------------------------
% David Romero 24/04/2020

% Set default params
center = 'mean';
superior = 75;
inferior = 25;
circleRadius = prctile(Data(:),75);
% circleRadius = max([max(prctile(Data1,90)) max(prctile(Data2,90))]);
plotDots = true;
plotInterval = true;
gridWidth = 0.1;
nPoints = 50;
nRing = 3;
jitter = 0.04;
applyColormap = true;

ColorBackground = 'none'; 
ColorGridBack = [1 1 1];%[220 220 220]./255;
ColorGrid = [151 172 191]./255;
ColorText = [90 100 91]./255;%ColorText = [120 150 160]./255;
ColorFill= [239 92 84]./255;

% Check input parameters
if nargin == 0
    error('SpiderPlot requires at least 1 input');
elseif nargin > 1
    param = varargin{end};
    % Adjust for user params
    if isfield(param,'center')
        center = param.center;
    end
    if isfield(param,'superior')
        superior = param.superior;
    end
    if isfield(param,'inferior')
        inferior = param.inferior;
    end
    if isfield(param,'circleRadius')
        circleRadius = param.circleRadius;
    end
    if isfield(param,'plotDots')
        plotDots = param.plotDots;
    end
    if isfield(param,'plotInterval')
        plotInterval = param.plotInterval;
    end
    if isfield(param,'ColorFill')
        ColorFill = param.ColorFill;
    end
    if isfield(param,'applyColormap')
        applyColormap = param.applyColormap;
    end
end

% Get max circle radius and number of directions
rhoMax = circleRadius;
nAngle = size(Data,2);
theta = linspace(0,2*pi - 2*pi/nAngle,nAngle);


% Prepare grid
set(gca,'Color',ColorBackground);
set(gca,'XColor','none');
set(gca,'YColor','none');
axis equal;

% Draw outer circle
[Xr,Yr] = pol2cart(linspace(0,2*pi - 2*pi/nPoints,nPoints),repmat(rhoMax,1,nPoints));
patch(Xr,Yr,ColorGridBack,'EdgeColor',ColorGrid,'LineWidth',gridWidth);hold on;
Ytext = rhoMax - rhoMax/nRing*0.1;
Xtext = (rhoMax/nRing*(nRing+1) - rhoMax/nRing*nRing)*0.05;
text(Xtext,Ytext,num2str(round(rhoMax/nRing*nRing,3,'significant')),'Color',ColorText,'FontWeight','Bold');

% Draw radial lines in each direction
[Xg,Yg] = pol2cart(theta,rhoMax*ones(1,nAngle));
for n=1:nAngle
    plot([0 Xg(n)],[0 Yg(n)],'Color',ColorGrid,'LineWidth',gridWidth); 
end

% Draw Rings with numeric text
for n=1:nRing-1
    [Xr,Yr] = pol2cart(linspace(0,2*pi - 2*pi/nPoints,nPoints),repmat(rhoMax/nRing*n,1,nPoints));
    plot([Xr Xr(1)],[Yr Yr(1)],'Color',ColorGrid,'LineWidth',gridWidth);

    Ytext = rhoMax/nRing*n - rhoMax/nRing*0.1;
    Xtext = (rhoMax/nRing*(n+1) - rhoMax/nRing*n)*0.05;
%     text(Xtext,Ytext,num2str(round(rhoMax/nRing*n,3,'significant')),'Color',ColorText,'FontWeight','Bold');
end

% aux = {Data1,Data2};
aux = {Data};
for d=1:length(aux)
    Data = aux{d};
    % hide points out of the main circle

%     Data(Data>rhoMax) = nan;

    % Compute center line and inferior and superior intervals
    switch center
        case 'mean'
            Valmean = mean(Data);
        case 'median'
            Valmean = median(Data);        
        case 'otherwise'
            error('Wrong value for center parameter');
    end

    ValInf = prctile(Data,inferior);
    ValSup = prctile(Data,superior);

    % Overlay dots
    if plotDots
        for n=1:nAngle
            [Xp,Yp] = pol2cart(theta(n)+jitter*randn(1,size(Data,1)),Data(:,n)');
            Rp = sqrt(Xp.^2+Yp.^2);
            Xp(Rp>rhoMax) = [];
            Yp(Rp>rhoMax) = [];
            scatter(Xp,Yp,5,'k','filled','MarkerFaceAlpha',0.1);
        end
    end
    
    % Plot data
    [X,Y] = pol2cart(theta,Valmean);
    [XInf,YInf] = pol2cart(theta,ValInf);
    [XSup,YSup] = pol2cart(theta,ValSup);
    if plotInterval
        patch([XSup XSup(1) X X(1)],[YSup YSup(1) Y Y(1)],ColorFill(d,:),'FaceAlpha',0.35,'EdgeColor','none');hold on;
        patch([X X(1) XInf XInf(1)],[Y Y(1) YInf YInf(1)],ColorFill(d,:),'FaceAlpha',0.35,'EdgeColor','none');
    end

    % Plotrange
    [XInf,YInf] = pol2cart(theta,prctile(Data,0.5));
    [XSup,YSup] = pol2cart(theta,prctile(Data,99.5));
    
%     [XInf,YInf] = pol2cart(theta,min(Data));
%     [XSup,YSup] = pol2cart(theta,max(Data));        
    if plotInterval
        patch([XSup XSup(1) X X(1)],[YSup YSup(1) Y Y(1)],ColorFill(d,:),'FaceAlpha',0.10,'EdgeColor','none');
        patch([X X(1) XInf XInf(1)],[Y Y(1) YInf YInf(1)],ColorFill(d,:),'FaceAlpha',0.10,'EdgeColor','none');
    end

%     [XInf,YInf] = pol2cart(theta,min(Data));
%     [XSup,YSup] = pol2cart(theta,max(Data));
%     plot([XInf XInf(1)],[YInf YInf(1)],'--r','LineWidth',1);
%     plot([XSup XSup(1)],[YSup YSup(1)],'--r','LineWidth',1);
    
    % Plot data
    p = plot([X X(1)],[Y Y(1)],'Color',ColorFill(d,:),'LineWidth',2.5);

    % Change colormap
    if applyColormap
        cd = [uint8(colormap(hot)*255) uint8(ones(256,1))].';

        minVal = min(ValInf);
        maxVal = max(ValSup);
        range =round((Valmean-minVal)./(maxVal-minVal).*255+1);
        range = [range range(1)];
        pause(1);
        set(p.Edge,'ColorBinding','interpolated','ColorData',cd(:,range));
    end


% Draw Rings with numeric text
for n=1:nRing-1
    Ytext = rhoMax/nRing*n - rhoMax/nRing*0.1;
    Xtext = (rhoMax/nRing*(n+1) - rhoMax/nRing*n)*0.05;
    text(Xtext,Ytext,num2str(round(rhoMax/nRing*n,3,'significant')),'Color',ColorText,'FontWeight','Bold');
end

xlim([-rhoMax rhoMax]);
ylim([-rhoMax rhoMax]);

end

