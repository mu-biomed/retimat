function iq = image_quality(I, metric, varargin)
% Compute an image quality metric for each bscan in I. 
%
%
% Input arguments (mandatory)
% ---------------------------
% * **I**:           2D or 3D matrix with bscan data. If 3D data the 3rd dimension is assumed to be the bscan index.
%            
% * **metric**:      Metric used to compute image quality. Accepted options
%
%   - 'mTCI': maximum tissue contrast index [3].
%   - 'snr': signal to noise ratio.
%   - 'psnr': peak to noise ratio [1].
%   - 'cnr': contrast to noise ratio [2].
%
%
% Input arguments (optional)
% --------------------------
% * **scanner**: String defining the OCT scanner (when using mTCI). If not provided an arbitrary value will be used.
%   
%   - 'Cirrus'
%   - 'RTVue'
%   - 'Spectralis'
%   - '3D-OCT-1000'
%                    
% * **seg**:           struct with segmentation used in 'snr','psnr','cnr'.
%
%  
% Output arguments
% ---------------- 
% * **iq**: Image quality metric.        
%  
%
% Notes
% -----
% The definition of SNR, CNR varies across the literature. Check the
% references for the precise definition implemented here.
%   
% Here SNR, PSNR and CNR rely on an accurate segmentation of both the ILM
% and the BM. If the image is pure nois, highly ocluded, or the 
% segmentation is wrong these metrics might not work.
%
% A more accurate noise estimation procedure may require acquiring a NOISE
% profiling or reference image [4].
%
%
% References
% ----------
%
% [1] Shirasawa, Objective Determination of Optimal Number of Spectral-
% Domain Optical Coherence Tomographic Images of Retina to Average, PLOS
% ONE, 2014. http://dx.doi.org/10.1371/journal.pone.0110550
%
% [2] Rico-Jimenez J J, Real-time OCT image denoising using a self-fusion
% neural network," Biomed. Opt. Express, 2022.
% https://doi.org/10.1364/BOE.451029
%
% [3] Huang, Signal Quality Assessment of Retinal Optical Coherence 
% Tomography Images, IOVS, 2012. Https://dx.doi.org/10.1167%2Fiovs.11-8755
%
% [4] Sahu, Statistical modeling and Gaussianization procedure based 
% de-speckling algorithm for retinal OCT images, Journal of Ambient
% Intelligence and Humanized Computing, 2018.
%
%
% Example
% ---------      
% Compute mTCI for a whole volume
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% .. code-block:: matlab
%
%   [~, ~, bscan] = read_vol(file);
%   mTCI = image_quality(bscan, 'mTCI', 'Spectralis');

if nargin == 1
    metric = 'mTCI';
end

if nargin > 3
    error('This function accepts a maximum of 3 arguments');
end

I = double(I);
% I(I==0) = nan; % Spectralis data has a peak in zero. The detection level
% is above the real 0 and therefore many values are mapped to 0 (ceil).

% Metrics relying on the segmentation (informed)
if any(strcmp(metric, {'psnr','cnr','snr'}))
    if nargin < 3
     error([metric 'requires segmentation data as an input.']);
    end
    
    seg = varargin{1};
        
    % Computation
    if ismatrix(I)
        iq = compute_informed(I, metric, seg.ILM, seg.BM);
    elseif ndims(I) == 3
        n_bscan = size(I, 3);
        iq = nan(1, n_bscan);
        for i_bscan=1:n_bscan
            ilm = seg.ILM(i_bscan,:);
            bm = seg.BM(i_bscan,:);
            iq(i_bscan) = compute_informed(I(:,:,i_bscan), metric, ilm, bm);    
        end
    end    
    
elseif strcmp(metric, 'mTCI')
    if nargin == 3
        scanner = varargin{1};            
    else
        scanner = 'unknown';
    end

    % Get noise profiling defaults depending on scanner
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
            warning(['Unknown scanner. Valid options: Cirrus, RTVue, Spectralis',...
                ', 3D-OCT-1000. Using 0.45 as an arbitrary default.']);
            cN1Bs = 0.45;
    end

    % Metric computation
    if ismatrix(I)
        iq = compute_mTCI(I, cN1Bs);
    elseif ndims(I) == 3
        n_bscan = size(I, 3);
        iq = nan(1, n_bscan);
        for i_bscan=1:n_bscan
            iq(i_bscan) = compute_mTCI(I(:,:,i_bscan), cN1Bs);    
        end
    end            
else
    error("Unknown metric. Valid options are: 'mTCI', 'psnr','cnr','snr'");
end

function met = compute_informed(I, metric, ilm, bm)
% Efficiently compute the retinal mask (foreground)
idx = repmat(1:size(I,1),size(I,2),1)';
mask_in = (idx >= ilm) & (idx <= bm);

switch metric
    case 'snr'
        signal = mean(I(mask_in), 'omitnan');
        noise = std(I(~mask_in), 'omitnan');
        met = 10*log10(signal^2/noise^2);
    
    case 'psnr'
        signal = max(I(mask_in));
        noise = std(I(~mask_in), 'omitnan');
        met = 10*log10(signal^2/noise^2);
        
    case 'cnr'
        mu_fore = mean(I(mask_in), 'omitnan');
        mu_back = mean(I(~mask_in), 'omitnan');

        var_fore = nanvar(I(mask_in));
        var_back = nanvar(I(~mask_in));
        
        met = abs(mu_fore-mu_back) / sqrt(var_fore+var_back); 
end

function mTCI = compute_mTCI(I, cN1Bs)
I = I(~isnan(I));

N = length(I);

% Get pdf estimation from histogram
[h, edges] = histcounts(I,200);
y_pdf = h/N;
x_pdf = (edges(1:end-1) + edges(2:end))/2;

% Get cdf estimation
[y_cdf, x_cdf] = ecdf(I);

% Get N1 (histogram mode)
[yN1, ind_N1] = max(y_pdf);
N1 = x_pdf(ind_N1(1));
% cN1 = sum(I <= N1)/N;

% Compute cN1s
ind_N1s = find(y_pdf >= 0.95*yN1);
N1s = x_pdf(ind_N1s(1));
cN1s = sum(I <= N1s)/N;

% Compute cN2
cN2 = cN1s*0.999/cN1Bs;

if cN2 > 1
    error("Unable to compute mTCI due to a cN2 threshold > 1.");
end

% Compute N2
ind_N2 = find(y_cdf >= cN2);

N2 = x_cdf(ind_N2(1));

% Compute N3
cN3 = 0.999;
ind_N3 = find(y_cdf >= cN3);
N3 = x_cdf(ind_N3(1));

% Compute mTCI
mTCI = (N3 - N1)/(N2 - N1);

% subplot(121);hold on;
% max_val = max(y_pdf);
% plot(x_pdf,y_pdf,'k','LineWidth',1.5);
% plot([N1 N1],[0 max_val],'--r'); text(N1,max_val/2,'N1');
% plot([N1s N1s],[0 max_val],'--r');text(N1s,max_val/2,'N1s');
% plot([N2 N2],[0 max_val],'--g');text(N2,max_val/2,'N2');
% plot([N3 N3],[0 max_val],'--m');text(N3,max_val/2,'N3');
% 
% subplot(122);hold on;
% plot(x_cdf,y_cdf,'k');
% plot([N1 N1 0],[0 cN1s cN1s],'--r');
% plot([N2 N2 0],[0 cN2 cN2],'--g');
% plot([N3 N3 0],[0 cN3 cN3],'--m');
% 
% disp('');