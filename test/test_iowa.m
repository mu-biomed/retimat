close all;clc;clearvars;

addpath(genpath('/home/david/GITHUB/retimat'));

file = '/home/david/Dropbox (MGEP)/projects/retimat/iowa_refsurfaces.xml';

[h,seg] = read_xml_iowa(file, true);

surf(h.X,h.Y,h.scale_z*(seg.BM-seg.ILM),'EdgeColor','none')

axis([-3 3 -3 3]);