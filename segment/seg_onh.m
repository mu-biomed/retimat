function mask = seg_onh(I, method)
%SEG_ONH Segment optic nerve head (onh) from gray fundus images
%
%   mask = seg_onh(I, method)
%   Returns a mask corresponding to the onh 
%
%   Input arguments:
%  
%   'I'              Grayscale fundus (slo) images.
%                    
%   'method'         Method to be used to segment the onh.
%                    Default: brute-force            
%  
%  
%   Output arguments:
%  
%   'mask'           Binary mask depicting onh.          
%  
%
%   Notes
%   -----
%   This function is designed to work with slo images acquired during OCT 
%   acquisition and may not work with color funduscopic images.
%
%   Example 1
%   ---------      
%   % Example description
%
%     [~,~,~,slo] = read_vol('myimage.vol');
%     mask = seg_onh(slo);
%     
%
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2021


% Convert to gray if color image is given
I = im2double(rgb2gray(I));

% Get complementary image
I = imcomplement(I);

% Unnecessary I think
% I = I.^2;I = I./max(I);

% Resize for search accuracy
n_point = 500;
I = imresize(I,[n_point n_point]);

% Get coordinates
[X,Y] = meshgrid(linspace(0,500,n_point),linspace(0,500,n_point));

switch method
    case 'brute_force'    
        n_grid = 30;
        [Xc,Yc] = meshgrid(linspace(250,450,n_grid),linspace(50,250,n_grid));
        rx = 30;
        ry = 60;
        RSS = nan(n_grid,n_grid);
        for n=1:n_grid
            for i=1:n_grid
                RSS(n,i) = get_RSS(X, Y, Xc(n,i), Yc(m,i), rx, ry);
            end
        end

    case 'gradient_descent'
        % Iteration configuration
        step = 2*1e4;
        dstep = 0.1;
        nIter = 80;
        
        [xc,yc,rx,ry] = gradient_descent(X,Y,I,n_iter,step,dstep);
        figure;
        subplot(121);
        plot(1:(nIter+1),rss_all);grid on;
        xlabel('Iteration');ylabel('RSS')
        subplot(122);
        
        imagesc(I);colormap('gray');hold on;
        plot_ellipse(X,Y,xc,yc,rx,ry,'g')

    case 'matlab_opt'
        myFun = @(X,Y,x,I)get_RSS(X,Y,x,I);
        if strcmp(eye,'OD')
            x0 = [300 200 30 30];
        else
            x0 = [200 200 30 30];
        end
        fun = @(x)myFun(X,Y,x,I);
        x = fminsearch(fun,x0);        
        error("Not yet fully implemented");
        x = fminsearch(@get_RSS,[350 200 60 60]);
    otherwise
        error("Not supported method");
end

end

function gradient_descent()

rss_all = nan(1,nIter);
plotImage = false;

% --- Initial guess ---
xc = 300;
yc = 200;
rx = 40;
ry = 40;
mask = ellipse_mask(X,Y,xc,yc,rx,ry);   % Compute mask 
rss = 1/(nPoints^2)*sum((I(:)-mask(:)).^2); % Compute error
rss_all(1) = rss;
   
for iter = 1:n_iter
        
    % Delta X
    xc1 = xc + dstep;
    yc1 = yc;
    rx1 = rx;
    ry1 = ry;    
    mask1 = getMask(X,Y,xc1,yc1,rx1,ry1);   % Compute mask 
    rss1 = 1/(nPoints^2)*sum((I(:)-mask1(:)).^2); % Compute error    
    drss1 = rss1-rss;
    dxc = drss1/(xc1-xc);

    % Delta Y
    xc2 = xc;
    yc2 = yc + dstep;
    rx2 = rx;
    ry2 = ry;    
    mask2 = getMask(X,Y,xc2,yc2,rx2,ry2);   % Compute mask 
    rss2 = 1/(nPoints^2)*sum((I(:)-mask2(:)).^2); % Compute error    
    drss2 = rss2-rss;
    dyc = drss2/(yc2-yc);
    
     % Delta rx
    xc3 = xc;
    yc3 = yc;
    rx3 = rx + dstep;
    ry3 = ry;    
    mask3 = getMask(X,Y,xc3,yc3,rx3,ry3);   % Compute mask 
    rss3 = 1/(nPoints^2)*sum((I(:)-mask3(:)).^2); % Compute error    
    drss3 = rss3-rss;
    drx = drss3/(rx3-rx);
 
%     % Delta ry
    xc4 = xc;
    yc4 = yc;
    rx4 = rx;
    ry4 = ry + dstep;    
    mask4 = getMask(X,Y,xc4,yc4,rx4,ry4);   % Compute mask 
    rss4 = 1/(nPoints^2)*sum((I(:)-mask4(:)).^2); % Compute error    
    drss4 = rss4-rss;
    dry = drss4/(ry4-ry);
    
    % Update coefficients
    xc5 = xc - step*dxc;
    yc5 = yc - step*dyc;
%     rx5 = rx - step*drx;
%     ry5 = ry - step*dry;  
    rx5 = rx;
    ry5 = ry;
    
    mask5 = getMask(X,Y,xc5,yc5,rx5,ry5); % compute mask
    rss5 = 1/(nPoints^2)*sum((I(:)-mask5(:)).^2); % compute error  
    rss_all(iter+1) = rss; % save value

    if plotImage
        f = figure(7);
        ax = subplot(121);
        imagesc(I);colormap(ax,'gray');hold on;
        plotEllipse(X,Y,xc,yc,rx,ry,'r');
        plotEllipse(X,Y,xc5,yc5,rx5,ry5,'g');
                
        subplot(122);
        surf(Xc,Yc,RSS,'EdgeColor','none');hold on;
        scatter3(xc,yc,rss,'r','filled');
%         scatter3(xc1,yc1,rss1,'b','filled'); 
%         scatter3(xc2,yc2,rss2,'b','filled');
        scatter3(xc3,yc3,rss5,'g','filled');
        
        xlabel('xc');ylabel('yc');        
        set(f,'position',[0 0 1000 400]);
   end

   
   % Save values from previous iteration
   xc = xc5;
   yc = yc5;
   rx = rx5;
   ry = ry5;
   rss = rss5;  
end
      

end
function mask = ellipse_mask(X,Y,xc,yc,rx,ry)
    R = sqrt((ry.*(X-xc)).^2+(rx.*(Y-yc)).^2);
    Rel = rx*ry;
    mask = double(R<=Rel);
end

function RSS = get_RSS(X,Y,xc, yc, rx, ry, I)

% Get ellipse mask and compute RSS
mask = ellipse_mask(X,Y,xc,yc,rx,ry);
n_point = size(X,1);
RSS = 1/(n_point^2)*sum((I(:)-mask(:)).^2); 
end

function plot_ellipse(X,Y,xc,yc,rx,ry,colour)
%PLOTELIPSE Summary of thXs functXon goes here
scatter(xc,yc,'filled',colour)
t = -pi:0.01:pi;
x = xc + rx*cos(t);
y = yc + ry*sin(t);
plot(x,y,colour);
mask = getMask(X,Y,xc,yc,rx,ry);
switch colour
    case 'r'
        h = imshow(cat(3,ones(size(X)),zeros(size(X)),zeros(size(X))));
    case 'g'
        h = imshow(cat(3,zeros(size(X)),ones(size(X)),zeros(size(X))));
    case 'b'
        h = imshow(cat(3,zeros(size(X)),zeros(size(X)),ones(size(X))));                
end
set(h,'AlphaData',0.2.*mask)
end

