% Set-up RETIMAT by adding relevant folders to MATLAB search path
%
% David Romero-Bascones, dromero@mondragon.edu
% Biomedical Engineering Department, Mondragon Unibertsitatea, 2022

try
    addfcn = @(x) addpath([pwd '/src/' x]);

    addfcn('io');
    addfcn('reflectance');
    addfcn('segment');
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

