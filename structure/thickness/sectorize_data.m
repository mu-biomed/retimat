function Z1 = sectorize_data(Grid,X0,Y0,Z0)


switch Grid.type
    case 'regular'
        
        % Get grid parameters (for readability)
        n_sect = Grid.n_sect;
        X = Grid.X;
        Y = Grid.Y;
        
        Z1 = nan(1,n_sect);

        % Compute average per sector
        iSect = 1;

        for i = 1:length(X)-1     
            for j=1:length(Y)-1
                % Compute sector mask
                mask1 = X0>=X(i) & X0<X(i+1);
                mask2 = Y0>=Y(j) & Y0<Y(j+1);
                mask = mask1 & mask2;   

                % Average values
                Z1(iSect) = nanmean(Z0(mask));            
                iSect = iSect+1;
            end
        end   
       
        
    case 'wedge'
        
        % Get grid parameters (for readability)
        n_sect = Grid.n_sect;
        rho = Grid.rho;
        theta0 = Grid.theta0;
        theta = Grid.theta;
               
        Z1 = nan(1,n_sect);

        % Change to polar
        [theta_s,rho_s] = cart2pol(X0,Y0);
        
        % Rotate all angles to set to position initial theta as theta = 0
        theta_s = theta_s - theta0;
        
        % Reconfigure to have only positive values
        theta_s(theta_s<0) = theta_s(theta_s<0)+2*pi;
    
        % Central region
        mask = abs(rho_s)<= rho(1);
        Z1(1) = nanmean(Z0(mask));

        % Loop through rings and angles
        i_sect = 2;
        for i = 1:length(rho)-1     
            for t=1:Grid.n_angle
                % Compute roi mask
                mask1 = abs(rho_s)>rho(i) & abs(rho_s)<=rho(i+1);
                mask2 = theta_s >= theta(t) & theta_s<theta(t+1);
                mask = mask1 & mask2;   
                
                Z1(i_sect) = nanmean(Z0(mask));

                i_sect = i_sect+1;
            end
        end   

end



