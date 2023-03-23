function [I_flat, shift, mask_retina] = flatten_retina(I, method, scale_z)
% Pilot rpe based flattening
M = size(I, 2);

rpe = nan(1,M);
switch method
    case 'middle'
        If = imgaussfilt(I, 1.5);

        mid = nan(1,M);
        for i=1:M
            pdf = If(:,i)/sum(If(:,i));
            cdf = cumsum(pdf);
            [~, mid(i)] = min(abs(cdf - 0.4));
            
            mid_point = round(mid(i));
            [~, rpe(i)] = max([zeros(mid_point-1,1); If(mid_point:end,i)]);
        end
        mask_retina = [];
    case 'mask'
        mask_retina = seg_retina(I, scale_z, 50, 'mean', false);
        n_pixel_rnfl = 1e-3*80/scale_z;
        for i=1:M            
            first_mask = find(mask_retina(:,i), 1);
            mid_point = round(first_mask + n_pixel_rnfl);
            [~, rpe(i)] = max([zeros(mid_point-1,1); I(mid_point:end,i)]);
        end        
    otherwise
        error("Not supported retina flattening method");        
end

% Second order polynomial
p = polyfit(1:M, rpe, 2);
rpe_2 = p(3) + p(2)*(1:M) + p(1)*(1:M).^2;

% Flatten the retina
shift = round(min(rpe_2) - rpe_2);

I_flat = I;
for i=1:M
    I_flat(:,i) = circshift(I_flat(:,i), shift(i));
end
