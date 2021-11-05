function Segment = remove_segmentation_outliers(Segment,th)

f = fields(Segment);

for i = 1:length(f)
    Segment.(f{i})(abs(Segment.(f{i})) > th) = nan;
end
