close all;clc;clearvars;

file = '/home/drombas/github/retimat/tutorials/data/example_raster.vol';

% addpath(genpath('/home/drombas/github/retimat/src'));
% [header, seg, bscan] = read_vol(file,'raw_pixel');
% 
% header.size_x = header.scale_x * (header.n_ascan - 1);
% header.size_y = header.scale_y * (header.n_bscan - 1);
% 
% segment_layers_aura(bscan, header);

%% AURA
filenames = {file};

params.resizedata = true;
params.minseg = false;
params.smooth = true;
params.segmethod = 1;
params.gridradii = [500 1500 2500];
params.logfile = false;
params.printtoscreen = true;
params.resultfolder = './Results';
params.overwrite_results = true;
params.saveXMLfiles = true;
params.displaygrid = true;
params.skip_completed = true;
params.displayresult = true;

seg_aura(file, params)
