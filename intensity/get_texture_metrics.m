function X = get_texture_metrics(P)


% Definitions from Tazarjani, 2021

% P: GLCM (2D probability density function)

energy = sum(P(:).^2);
entropy = -sum(P(:).*log(P(:)));

[i,j] = meshgrid(1:size(P,1), 1:size(P,2));

contrast = sum((i(:) - j(:)).^2 .* P(:));
homogeneity = sum(P(:) ./ (1 + (i(:) - j(:)).^2));

% correlation = sum(i(:)
% sum_of_squares = 

