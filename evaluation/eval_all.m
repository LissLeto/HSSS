function [BR, IBR] = eval_all(labels, Seg, fig_ori, is_disp)
%EVAL_ALL Compute boundary recall metrics for superpixel segmentation.
%
%   [BR, IBR] = EVAL_ALL(labels, Seg) computes:
%       BR  - Boundary Recall over all ground-truth boundary pixels.
%       IBR - Instance Boundary Recall, i.e. the average Boundary Recall
%             computed separately for each ground-truth region.
%
%   [BR, IBR] = EVAL_ALL(labels, Seg, fig_ori, is_disp) optionally
%   displays a boundary-recall visualization when is_disp is true. In the
%   visualization, missed GT boundaries are red, recalled GT boundaries are
%   black, and superpixel-only boundaries are green.
%
%   Inputs:
%       labels  - H-by-W superpixel label map. Labels may be non-consecutive
%                 and may start from 0 or 1.
%       Seg     - H-by-W ground-truth segmentation label map. Labels may be
%                 non-consecutive and may start from 0 or 1.
%       fig_ori - Optional H-by-W or H-by-W-by-3 image for visualization.
%       is_disp - Optional logical flag. Default: false.
%
%   The boundary tolerance follows the original implementation:
%       r = round(0.0025 * sqrt(H^2 + W^2)).

if nargin < 4 || isempty(is_disp)
    is_disp = false;
end
if nargin < 3
    fig_ori = [];
end

validateattributes(labels, {'numeric', 'logical'}, {'2d', 'nonempty'}, mfilename, 'labels', 1);
validateattributes(Seg, {'numeric', 'logical'}, {'2d', 'nonempty'}, mfilename, 'Seg', 2);

if ~isequal(size(labels), size(Seg))
    error('eval_all:SizeMismatch', 'labels and Seg must have the same height and width.');
end

[M, N] = size(Seg);
r = round(0.0025 * sqrt(M^2 + N^2));

gtBoundary = boundarymask(Seg, 4);
spBoundary = boundarymask(labels, 4);

kernel = true(2 * r + 1);
matchedBoundary = conv2(double(spBoundary), double(kernel), 'same') > 0;

gtBoundaryCount = nnz(gtBoundary);
if gtBoundaryCount == 0
    BR = NaN;
else
    BR = nnz(gtBoundary & matchedBoundary) / gtBoundaryCount;
end

gtLabels = unique(Seg(:));
BR_each = NaN(numel(gtLabels), 1);

for k = 1:numel(gtLabels)
    instanceBoundary = gtBoundary & (Seg == gtLabels(k));
    boundaryCount = nnz(instanceBoundary);

    if boundaryCount > 0
        BR_each(k) = nnz(instanceBoundary & matchedBoundary) / boundaryCount;
    end
end

IBR = mean(BR_each, 'omitnan');

if is_disp
    if isempty(fig_ori)
        error('eval_all:MissingImage', 'fig_ori is required when is_disp is true.');
    end
    show_boundary_recall(fig_ori, gtBoundary, spBoundary, matchedBoundary);
end

end

function show_boundary_recall(fig_ori, gtBoundary, spBoundary, matchedBoundary)
rgb = prepare_display_image(fig_ori);

missedGt = gtBoundary & ~matchedBoundary;
recalledGt = gtBoundary & matchedBoundary;
spOnly = spBoundary & ~gtBoundary;

rgb = set_rgb_pixels(rgb, spOnly, [0, 255, 0]);
rgb = set_rgb_pixels(rgb, recalledGt, [0, 0, 0]);
rgb = set_rgb_pixels(rgb, missedGt, [255, 0, 0]);

figure;
imshow(rgb);
title('Boundary Recall');
end

function rgb = prepare_display_image(img)
if ndims(img) == 2
    img = repmat(img, 1, 1, 3);
end

if size(img, 3) ~= 3
    error('eval_all:InvalidImage', 'fig_ori must be H-by-W or H-by-W-by-3.');
end

if isfloat(img)
    if max(img(:)) <= 1
        rgb = uint8(round(255 * min(max(img, 0), 1)));
    else
        rgb = uint8(min(max(round(img), 0), 255));
    end
else
    rgb = uint8(img);
end
end

function rgb = set_rgb_pixels(rgb, mask, color)
for c = 1:3
    channel = rgb(:, :, c);
    channel(mask) = color(c);
    rgb(:, :, c) = channel;
end
end
