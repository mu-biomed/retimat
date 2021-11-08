function R = get_reflectance_map(bscan, top, bottom, metric)


[~, n_ascan, n_bscan] = size(bscan);

MR = nan(n_bscan, n_ascan);

for i_bscan=1:n_bscan
    
    MR(i_bscan,:) = get_layer_intensity(bscan(:,:,i_bscan), top(i_bscan,:), ...
        bottom(i_bscan,:));    
end

switch metric
    case 'mean'
        R = MR;
    
    case 'total'
        Thickness = get_thickness_map(top, bottom, header);
        
        TR = MR.*Thickness./header.scale_z;        
        R = TR;                
    
    case 'layer_index'
        LI = nan(n_bscan, n_ascan);
        
        Thickness = get_thickness_map(top, bottom, header);
        
        for i_bscan=1:n_bscan
            Isa = prctile(reshape(bscan(:,:,i_bscan),1,[]), 99);
            LI(i_bscan,:) = MR(i_bscan,:).*Thickness(i_bscan,:)/Isa;
        end        
        R = LI;
end

