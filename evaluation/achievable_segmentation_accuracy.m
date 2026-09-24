function asa = achievable_segmentation_accuracy(labels, segmentation)
%ACHIEVABLE_SEGMENTATION_ACCURACY Compute the optional ASA metric.
%
% ASA = ACHIEVABLE_SEGMENTATION_ACCURACY(LABELS, SEGMENTATION) computes ASA
% for one superpixel label map and one ground-truth segmentation. BR and IBR
% are the default metrics used by the bundled evaluation demo.

validateattributes(labels, {'numeric', 'logical'}, ...
    {'2d', 'nonempty', 'finite'}, mfilename, 'labels', 1);
validateattributes(segmentation, {'numeric', 'logical'}, ...
    {'2d', 'nonempty', 'finite'}, mfilename, 'segmentation', 2);
if ~isequal(size(labels), size(segmentation))
    error('achievable_segmentation_accuracy:SizeMismatch', ...
        'labels and segmentation must have the same size.');
end

[~, ~, superpixels] = unique(labels(:), 'sorted');
[~, ~, regions] = unique(segmentation(:), 'sorted');
overlap = accumarray([superpixels, regions], 1);
asa = sum(max(overlap, [], 2)) / numel(labels);
end
