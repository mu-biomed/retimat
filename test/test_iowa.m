close all;clc;clearvars;
addpath(genpath('C:\Users\dromero\Desktop\GITHUB\retimat'));

%% Star
in_dir = '../data_private/iowa/vol_star';
images = {'HC001AM_H_2341', 'HC001AM_H_2342',...
          'HC002AM_H_2347', 'HC002AM_H_2348',...
          'HC002AF_A_2357', 'HC002AF_A_2358'};
n_image = length(images);
      
layers = {'TRT','RNFL','GCIP','INL'};
n_layer = length(layers);

for i_image=1:n_image
    image = images{i_image};
    
    file_vol  = [in_dir '/' image '/' image '_OCT_Iowa.vol'];
    file_iowa = [in_dir '/' image '/' image '_Surfaces_Iowa.xml'];
        
    [h1,seg1,bscan,slo] = read_vol(file_vol, 'coordinates');
    T1 = compute_thickness(seg1,layers,h1.scale_z);

    % IOWA loading
    [h2,seg2] = read_xml_iowa(file_iowa, 'coordinates');
    T2 = compute_thickness(seg2, layers,h2.scale_z);

    for i_layer=1:n_layer
        layer = layers{i_layer};

        subplot(2,n_layer,i_layer);
        surf(h1.X_oct,h1.Y_oct,T1.(layer),'EdgeColor','none')
        view(0,90);
        axis([-3 3 -3 3]);
        title([layer ' - Spectralis']);

        subplot(2,n_layer,i_layer+n_layer);
        surf(h1.X_oct,h1.Y_oct,T2.(layer),'EdgeColor','none')
        view(0,90);
        axis([-3 3 -3 3]);
        title([layer ' - IOWA']);
    end
end
      
%% Raster
in_dir = '../data_private/iowa/vol_raster';

images = {'HC001AF_H_2644', 'HC001AF_H_2645',...
          'HC001AM_H_2339', 'HC001AM_H_2340',...
          'HC002AF_A_2888', 'HC002AF_A_2889'};
n_image = length(images);

layers = {'TRT','RNFL','GCIP','INL'};
n_layer = length(layers);

for i_image=1:n_image
    image = images{i_image};
    
    file_vol  = [in_dir '/' image '/' image '_OCT_Iowa.vol'];
    file_iowa = [in_dir '/' image '/' image '_Surfaces_Iowa.xml'];
        
    [h1,seg1,bscan,slo] = read_vol(file_vol, 'coordinates');
    T1 = compute_thickness(seg1,layers,h1.scale_z);

    % IOWA loading
    [h2,seg2] = read_xml_iowa(file_iowa, 'coordinates');
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
        title([layer ' - IOWA']);
    end
end