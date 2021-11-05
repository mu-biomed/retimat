function Grid = configure_grid(Grid)

switch Grid.type
    case 'regular'
        % Compute coordinates (as 1D arrays)
        Grid.X = linspace(-Grid.maxD,Grid.maxD,Grid.n_square+1);
        Grid.X(end)= Grid.X(end)+0.00000001; % just to get all points
        Grid.Y = Grid.X;        
        % Get sector number
        Grid.n_sect = Grid.n_square^2;
        % Get name
        Grid.name = ['regular_3_' num2str(1e3*2*Grid.maxD/Grid.n_square)];
 
    case 'wedge'
        % Compute polar coordinates
        Grid.theta = linspace(0,2*pi,Grid.n_angle+1);
        Grid.theta(end) = Grid.theta(end)+1e-7;
        if ~isfield(Grid,'theta0')
            Grid.theta0 = -pi/Grid.nAngle;
        end
        % Get sector number
        Grid.n_sect = 1+(length(Grid.theta)-1) *(length(Grid.rho)-1);
        % Get name
        Grid.name = ['wedge_' num2str(Grid.n_sect) '_sect_' num2str(Grid.n_angle) '_ang']; 
        
    otherwise
        error('Not supported grid type');
end