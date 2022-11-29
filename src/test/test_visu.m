close all;clc;clearvars;
addpath(genpath('..'));

% Read a raster file
% [header, seg, bscan, fundus] = read_vol('../data/raster.vol');
% T = compute_thickness(seg, 'TRT');
n_col = 4;

%% Sectors: pie
subplot(1,n_col,1);
Sect.n_sect  = 4;
Sect.type    = 'pie';
Sect.theta_0 = -pi/4;
Sect.radius  = 3;
Z = [1 2 3 4];
plot_sectors(Z, Sect, 'axis_off', false);

%% Sectors: ETDRS
subplot(1,n_col,2);
Z = [1 2 3 4 5 6 7 8 9];
Sect.n_sect  = 9;
Sect.type    = 'wedge';
Sect.n_angle = 4;
Sect.theta_0 = -pi/4;
Sect.radius  = [0.5 1.5 3];
plot_sectors(Z, Sect,'axis_off', false);

%% Sectors: rings
subplot(1,n_col,3);
Z = [1 2 3 4 5];
Sect.n_sect  = 5;
Sect.type    = 'ring';
Sect.radius  = linspace(0.5, 3, Sect.n_sect + 1);
plot_sectors(Z, Sect, 'edge_color','none');

%% Sectors: regular
subplot(1,n_col,4);
Z = randn(30,30);
Sect.type    = 'regular';
Sect.n_x     = 30;
Sect.n_y     = 30;
Sect.X_edge  = linspace(-3, 3, 30+1);
Sect.Y_edge  = linspace(-3, 3, 30+1);
plot_sectors(Z, Sect, 'edge_color','none','alpha',0.1);

%% Generate report
% generate_report(bscan, seg, fundus, {'TRT','GCIPL','INL','RNFL'}, 'n_plot_bscan',10, 'n_col_max', 7);