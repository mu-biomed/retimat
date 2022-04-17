function iq = image_quality(I, metric, varargin)
%IMAGE_QUALITY Compute an image quality metric
%
%   metric = image_quality(I)
%   Compute an image quality metric for each bscan in I. Supported options
%   are signal to noise ratio (snr) and maximum tissue contrast index 
%   (mTCI) [1].
%
%   Input arguments:
%  
%   'I'              2D or 3D matrix with bscan data. If 3D data the 3rd
%                    dimension is assumed to be the bscan index.
%            
%   'metric'         Metric used to compute image quality. Accepted options
%                    Options: 'snr','snr_approx','mTCI'
%                    Default: 'mTCI'
%
%   Optional input arguments (varargin):
%
%   'scanner'        When using mTCI. String defining the OCT scanner.
%                    Options: 'Cirrus','RTVue',Spectralis','3D-OCT-1000'
%                    If not provided an arbitrary value will be used.
%   
%   'seg'            Segmentation structure. When using 'snr'.
%
%   'scale_z'        Axial resolution. Used in 'snr_approx' to approximate
%                    the number of pixels of the retina.
%
%  
%   Output arguments:
%  
%   'iq'         Image quality metric.        
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
%   % Compute mTCI for a whole volume
%
%   [~, ~, bscan] = read_vol(file);
%   mTCI = image_quality(bscan, 'Spectralis');
%
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2022

if nargin == 1
    metric = 'mTCI';
end

if nargin > 3
    error('This function accepts a maximum of 3 arguments');
end

I = double(I);
% I(I==0) = nan; % Spectralis data has a peak in zero. The detection level
% is above the real 0 and therefore many values are mapped to 0 (ceil).

switch metric
    case 'mTCI'
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
        
    case 'snr'
        if nargin < 3
            error('snr requires segmentation data as an input');
        end        
        seg = varargin{1};
        
        % SNR computation
        if ismatrix(I)
            iq = compute_snr(I, seg.ILM, seg.BM);
        elseif ndims(I) == 3
            n_bscan = size(I, 3);
            iq = nan(1, n_bscan);
            for i_bscan=1:n_bscan
                ilm = seg.ILM(i_bscan,:);
                bm = seg.BM(i_bscan,:);
                iq(i_bscan) = compute_snr(I(:,:,i_bscan), ilm, bm);    
            end
        end
        
    otherwise
        error("Unknown metric. Valid options are: 'mTCI', 'snr'");
end

function snr = compute_snr(I, ilm, bm)

% Efficiently compute the retinal mask
idx = repmat(1:size(I,1),size(I,2),1)';
mask_in = (idx >= ilm) & (idx <= bm);

signal = nanmean(I(mask_in));
noise = nanmean(I(~mask_in));
snr = signal/noise;

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