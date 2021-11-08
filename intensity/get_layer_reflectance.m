function intensity = get_layer_intensity(I, top, bottom)
% Get layer intensity
% get_layer_intensity(I, top, bottom)
%
% Biomedical Engineering Department, Mondragon Unibertsitatea, 2021
% David Romero Bascones dromero@mondragon.edu

n_col = size(I, 2);

intensity = nan(1,n_col);

for i_col=1:n_col
    roi = round(top(i_col):bottom(i_col));
    intensity(i_col) = mean(I(roi, i_col));
    
end