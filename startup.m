% Set-up RETIMAT by adding relevant folders to MATLAB search path
%
% David Romero-Bascones, dromero@mondragon.edu
% Biomedical Engineering Department, Mondragon Unibertsitatea, 2024

try
    wdir = replace(which(mfilename), 'startup.m', '');    
    
    addfcn = @(x) addpath(fullfile(wdir, 'src', x));

    addfcn('io');
    addfcn('reflectance');
    % addfcn('segment');
    addfcn('spatial');
    addfcn('structure');
    addfcn('texture');
    addfcn('utils');
    addfcn('visu');
    addfcn('external');
    
    disp('RETIMAT is ready to be used.');
catch emsg
    disp('Error setting up RETIMAT.');
    disp(emsg.message);
end

