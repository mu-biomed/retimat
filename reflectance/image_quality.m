function metric = image_quality(I, scanner)
%IMAGE_QUALITY Compute an image quality metric
%
%   metric = image_quality(I)
%   Compute the maximum tissue contrast index (mTCI) from each bscan [1].
%
%   Input arguments:
%  
%   'I'              2D or 3D matrix with bscan data. If 3D data the 3rd
%                    dimension is assumed to be the bscan index.
%            
%   'scanner'        String definint the OCT scanner.
%                    Options: 'Cirrus','RTVue',Spectralis','3D-OCT-1000'
%                    If not provided an arbitrary value will be used.
%
%
%  
%   Output arguments:
%  
%   'metric'         Image quality metric.        
%  
%
%
%   References
%   ----------
%   [1] Huang, Signal Quality Assessment of Retinal Optical Coherence 
%   Tomography Images, IOVS, 2012. 
%   https://dx.doi.org/10.1167%2Fiovs.11-8755
%
%
%   Example
%   ---------      
%   % Example description
%
%     [~,
%     [GLCMS,SI] = graycomatrix(I,'NumLevels',9,'G',[])
%     
%
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2022

if nargin == 1
    scanner = 'unknown';
end

switch scanner
    case 'Cirrus'
        cN1Bs = 0.4;
    case 'RTVue'
        cN1Bs = 0.55;
    case 'Spectralis'
        cN1Bs = 0.45;
    case '3D-OCT-1000'
        cN1Bs = 0.5;
    otherwise
        warning('Unknown scanner. Using 0.45 as an arbitrary default');
        cN1Bs = 0.45;
end

if ndims(I) == 3
    metric = nan(1, n_bscan);
    for i_bscan=1:n_bscan
        metric(i_bscan) = compute_mTCI(I(:,:,i_bscan), cN1Bs);    
    end
    return;
end

metric = compute_mTCI(I, cN1Bs);

function mTCI = compute_mTCI(I, cN1Bs)
I = I(:);
N = length(I);

h = histogram(I);
y_pdf = h.Values;
x_pdf = h.BinEdges(1:end-1) + h.BinWidth/2;

[y_cdf, x_cdf] = ecdf(I);

% Get N1 (histogram mode)
[yN1, ind_N1] = max(y_pdf);
N1 = x_pdf(ind_N1(1));

% Compute cN1s
ind_N1s = find(y_pdf >= 0.95*yN1);
N1s = x_pdf(ind_N1s(1));
cN1s = sum(I <= N1s)/N;

% Compute cN2
cN2 = cN1s*0.999/cN1Bs;

% Compute N2
ind_N2 = find(y_cdf >= cN2);
N2 = x_cdf(ind_N2(1));

% Compute N3
ind_N3 = find(y_cdf >= 0.999);
N3 = x_cdf(ind_N3(1));

% Compute mTCI
mTCI = (N3 - N1)/(N2 - N1);