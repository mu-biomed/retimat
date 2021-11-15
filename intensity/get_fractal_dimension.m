function FD_map = get_fractal_dimension(I)


% Omar S. Al-Kadi (e-mail: o.al-kadi@sussex.ac.uk; o.alkadi@ju.edu.jo)
% University of Sussex, Brighton, UK.
%
% Please cite the following paper if you use this file for your research:
% O. S. Al-Kadi and D. Watson, �Texture Analysis of Aggressive and non-Aggressive Lung Tumor CE CT Images,� IEEE Transactions on Biomedical Engineering, vol. 55, pp. 1822-1830, 2008.
%----------------------------------------------------------------------------

I = imread('bricks.jpg');

[M,N]= size(I);

%------- performing non-linear filtering on a varying size pixel block -------%
h = waitbar(0,'Performing 3-D Box Counting...');
B = cell(1,6);

for r=2:7
    rc = @(x) floor(((max(x)-min(x))/r))+ 1; % non-linear filter
    F= colfilt(I, [r r],'sliding', rc);
    B{r}= log(double(F * (49/(r^2))));
    waitbar(r/6)
end
close(h)

i=log(2:7); % Normalised scale range vector
%------- computing the slope using linear regression -------%
Nxx=dot(i,i)-(sum(i)^2)/6;
h = waitbar(0,'Transforming to FD...');

FD_map = nan(M,N);
for m = 1:M
    for n = 1:N
        fd= [B{7}(m,n), B{6}(m,n), B{5}(m,n), B{4}(m,n), B{3}(m,n), B{2}(m,n)]; % Number of boxes multiscale vector
        Nxy=dot(i,fd)-(sum(i)*sum(fd))/6; 
        FD_map(m, n)= (Nxy/Nxx); % slope of the linear regression line
    end
    waitbar(m/M)
end
close(h)

% Fractal dimension features
% mean, std, lacunarity