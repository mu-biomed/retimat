close all;clc;clearvars;

% addpath(genpath('/home/david/GITHUB/retimat'));
addpath(genpath('C:\Users\dromero\Desktop\GITHUB\retimat'));

sub_eye = 'E07_HD_E_7574';
% Spectralis loading

layers = {'TRT','RNFL','GCIP','INL'};
n_layer = length(layers);

file = ['C:\Users\dromero\Desktop\test_iowa\' sub_eye '.vol'];
[h1,seg1,bscan,slo] = read_vol(file,'coordinates');

T1 = compute_thickness(seg1,layers,h1.scale_z);

% IOWA loading
file = ['C:\Users\dromero\Desktop\test_iowa\' sub_eye '\' sub_eye '_Surfaces_Iowa.xml'];
[h2,seg2] = read_xml_iowa(file, true);

T2 = compute_thickness(seg2, layers,h2.scale_z);

for i_layer=1:n_layer
    layer = layers{i_layer};
    
    subplot(2,n_layer,i_layer);
    surf(h1.X_oct,h1.Y_oct,T1.(layer),'EdgeColor','none')
    view(0,90);
    axis([-3 3 -3 3]);
    title([layer ' - Spectralis']);
    
    subplot(2,n_layer,i_layer+n_layer);
    surf(h2.X_oct,h2.Y_oct,T2.(layer),'EdgeColor','none')
    view(0,90);
    axis([-3 3 -3 3]);
    title([layer ' - AURA']);
end