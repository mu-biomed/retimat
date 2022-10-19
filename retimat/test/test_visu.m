close all;clc;clearvars;

addpath(genpath('..'));

%% Plot sectors (pie)
Z = [1 2 3 4];
Sect.n_sect = 4;
Sect.type = 'pie';
Sect.theta_0 = -pi/4;
Sect.radius = 3;

plot_sectors(Z,Sect);