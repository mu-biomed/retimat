function X = get_LBP_features(I, n_neighbor)

LBP = extractLBPFeatures(I, 'NumNeighbors', n_neighbor);
% extra parameters can be entered (radius, rotation-invariance)

X.mean = mean(LBP);
X.median = median(LBP);
X.std = std(LBP);
X.iqr = iqr(LBP);
X.range = max(LBP) - min(LBP);
X.skewness = skewness(LBP);
X.kurtosis = kurtosis(LBP);
