function [scale_x_fun, scale_y_fun] = compute_biometrics(header, sex, method)

% Get scan focus ~= SE
scan_focus = header.ScanFocus;

% Get fundus parameters
fov = double(header.FieldSizeSlo);
nx = double(header.SizeXSlo);
ny = double(header.SizeYSlo);


switch method
    case 'Spectralis'
        % As per Ctori, 2015
        % Gullstrand eye model 
        % AL = 24.3854 mm
        % Ave_K = 7.7 mm (predefined unless entered at acquisition)
        % Refractive error measured when focusing (ScanFocus)

        scale_x_fun = header.ScaleXSlo;
        scale_y_fun= header.ScaleYSlo;
        
    case 'Littman'
        if strcmp(sex,'M')
            Ave_K = 7.72;  % Based on [HyeongSu, PlOS, 2019]
        elseif strcmp(sex,'F')  
            Ave_K = 7.63;  % Based on [HyeongSu, PlOS, 2019]
        end
        
        AL = 24*Ave_K/7.8 - 0.4*scan_focus;  % Littman formula
        f0 = 7.2; % Gullstrand eye model (24.3854 mm - 17.185 mm)  [Atchinson pp. 250-251]
        f_nodal = AL - f0;
        scale_x_deg = fov/(nx-1); % angular resolution (degrees between two adjacent fundus pixels)
        scale_y_deg = fov/(ny-1); % angular resolution (degrees between two adjacent fundus pixels)
        scale_x_fun = f_nodal * (scale_x_deg)*pi/180; % s = r*phi
        scale_y_fun = f_nodal * (scale_y_deg)*pi/180;
end