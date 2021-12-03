function Inorm = normalize_reflectance(bscan, top, bottom)


[~, n_ascan, n_bscan] = size(bscan);

I = nan(n_bscan, n_ascan);

for i_bscan=1:n_bscan
   
    % Get reference layer reflectance
    ref_reflectance = get_layer_intensity(bscan(:,:,i_bscan), top(i_bscan,:), ...
        bottom(i_bscan,:));
    
    ref_reflectance = mean(ref_reflectance);
    
    bscan(:,:,i_bscan) = bscan(:,:,i_bscan)/ref_reflectance;
    
end