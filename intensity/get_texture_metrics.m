function X = get_texture_metrics(P)


% Definitions from Tazarjani, 2021

% P: GLCM (2D probability density function)

N = size(P,1);
Pn = P./sum(P(:)); % Normalize GLCM matrix

[i,j] = meshgrid(1:N, 1:N);

mu = sum(i(:) .* Pn(:));
sd = sqrt(sum(Pn(:) .* (i(:) - mu).^2));

mux = sum(i(:) .* Pn(:));
muy = sum(j(:) .* Pn(:));
sdx = sqrt(sum(Pn(:) .* (i(:) - mux).^2));
sdy = sqrt(sum(Pn(:) .* (j(:) - muy).^2));

P_sum = zeros(1, 2*N-1);
k = 2:2*N;
for ik=1:length(k)
    mask = double((i + j) == k(ik));
    P_sum(ik) = sum(Pn(:).*mask(:));
end

P_dif = zeros(1, N);
d = 0:N-1;
for id=1:length(d)
    mask = double(abs(i - j) == d(id));
    P_dif(id) = sum(Pn(:).*mask(:));
end


% ----------------------------- Features ------------------------------------
energy = sum(Pn(:).^2);
% energy = graycoprops(P, 'energy').Energy;

entropy = -sum(Pn(:).*log(Pn(:)));
contrast = sum((i(:) - j(:)).^2 .* Pn(:)); 
% contrast = graycoprops(P, 'contrast').Contrast;

homogeneity = sum(Pn(:) ./ (1 + (i(:) - j(:)).^2));  % Likely to be wrong on Tazarjani,2021 (same dfinition as inverse moment)
homogeneity2 = graycoprops(P, 'homogeneity').Homogeneity; % different definition (abs vs squared)

correlation = sum(Pn(:) .* (i(:) - mux)./sdx .* (j(:) - muy)./sdy); 
% correlation = graycoprops(P, 'correlation').Correlation;

sum_of_squares = sum(Pn .* (i(:) - mu).^2);
cluster_shade = sum(Pn(:) .* (i(:) + j(:) -2*mu).^4);
cluster_prominence = sum(Pn(:) .* (i(:) + j(:) -2*mu).^3);
dissimilarity = sum(Pn(:) .* abs(i(:) - j(:)));
autocorrelation = sum(Pn(:) .* i(:) .* j(:));

sum_average = sum(P_sum .* k);
sum_entropy = -sum(P_sum .* log(P_sum));
sum_variance = sum(P_sum .* (k - sum_average).^2);

inverse_difference = sum(Pn(:) ./ (1 + abs(i(:) - j(:))));
inverse_difference_moment = sum(Pn(:) ./ (1 + (i(:) - j(:)).^2));

dif_average = sum(P_dif .* d);
dif_entropy = -sum(P_dif .* log(P_dif));
dif_variance = sum(P_dif .* (d - dif_average).^2);

max_prob = max(Pn(:));

% Hxy = -


