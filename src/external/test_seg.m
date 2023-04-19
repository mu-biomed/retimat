close all;clc;clearvars;

% file = '/home/drombas/github/retimat/tutorials/data/example_raster.vol';

file = 'C:\Users\dromero\Desktop\GITHUB\retimat\tutorials\data\example_raster.vol';
addpath(genpath('C:\Users\dromero\Desktop\GITHUB\retimat'));

% addpath(genpath('/home/drombas/github/retimat/src'));
% [header, seg, bscan] = read_vol(file,'raw_pixel');
% 
% header.size_x = header.scale_x * (header.n_ascan - 1);
% header.size_y = header.scale_y * (header.n_bscan - 1);
% 
% segment_layers_aura(bscan, header);

%% AURA
[h, ~, bscan] = read_vol(file, 'raw_pixel');

h.SizeX = h.n_ascan;
h.SizeZ = h.n_axial;
h.NumBScans = h.n_bscan;
h.ScaleX = h.scale_x;
h.ScaleZ = h.scale_z;
h.Distance = h.scale_y;
h.ScanPosition = h.eye;
header = h;

params.resizedata = true;
params.minseg = false;
params.smooth = true;
params.segmethod = 1;
params.printtoscreen = true;
params.skip_completed = true;
params.displayresult = true;

segment_layers(bscan, h, params)
