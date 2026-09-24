function [label_maps, info] = hsss_segment(rgb_image, target_counts, options)
%HSSS_SEGMENT Hierarchical Superpixel Segmentation by Searching Seeds.
%
% [LABEL_MAPS, INFO] = HSSS_SEGMENT(IMAGE, TARGET_COUNTS, OPTIONS)
% returns an H-by-W-by-numel(TARGET_COUNTS) array. IMAGE may be an RGB
% array or an image path. TARGET_COUNTS are returned in descending order.
%
% Important OPTIONS fields:
%   EdgeMap          External edge map array/path. Empty uses Sobel.
%   EdgePercentile   Edge threshold percentile (default 80).
%   MaxL             Upper bound used by Getpara2 (default 1).
%   DGap             Hierarchy pruning multiplier (default 1).
%   RandomSeed       Deterministic seed for tie breaking (default 0).
%   SpatialWeight    Pixel-to-seed spatial weight (default 0.5).
%   SearchRadiusMultiplier Candidate-seed radius multiplier (default 5).
%   MaxMergeRounds   Safety limit (default 30).
%   SeriesLevel      Level passed to the original hierarchy (default 1).

if nargin < 3
    options = struct();
end
options = apply_defaults(options);
target_counts = validate_targets(target_counts);

if ischar(rgb_image) || isstring(rgb_image)
    rgb_for_preparation = imread(rgb_image);
else
    rgb_for_preparation = rgb_image;
end

timer_id = tic;
[labels, lab_image, para, seeds, edge_map, foundation_info] = ...
    prepare_hsss_image(rgb_for_preparation, options);
if any(target_counts > foundation_info.InitialSuperpixels)
    error('HSSS:TargetAboveFoundation', ...
        'A requested count exceeds the %d foundation superpixels.', ...
        foundation_info.InitialSuperpixels);
end

requested = target_counts;
label_maps = zeros(size(labels, 1), size(labels, 2), 0);
rounds = 0;
while ~isempty(requested) && rounds < options.MaxMergeRounds
    rounds = rounds + 1;
    labels_before = labels;
    requested_before = requested;
    [labels, label_maps, requested, para] = series_SP_fast( ...
        labels, lab_image, requested, para, options.SeriesLevel, ...
        label_maps, options.DGap);
    if isequal(labels, labels_before) && isequal(requested, requested_before)
        break
    end
end
if ~isempty(requested)
    error('HSSS:TargetNotReached', ...
        'Could not reach target count(s): %s.', mat2str(requested));
end

info = foundation_info;
info.TargetCounts = target_counts;
info.ActualCounts = squeeze(max(max(label_maps, [], 1), [], 2)).';
info.MergeRounds = rounds;
info.RuntimeSeconds = toc(timer_id);
info.SeedMask = seeds;
info.EdgeMap = edge_map;
info.Method = "original";
end

function options = apply_defaults(options)
defaults = struct('EdgeMap', [], 'EdgePercentile', 80, 'MaxL', 1, ...
    'DGap', 1, 'RandomSeed', 0, 'MaxMergeRounds', 30, ...
    'SeriesLevel', 1, 'SpatialWeight', 0.5, ...
    'SearchRadiusMultiplier', 5);
names = fieldnames(defaults);
for idx = 1:numel(names)
    if ~isfield(options, names{idx})
        options.(names{idx}) = defaults.(names{idx});
    end
end
end

function targets = validate_targets(target_counts)
validateattributes(target_counts, {'numeric'}, ...
    {'vector', 'integer', 'positive', 'finite'});
targets = sort(unique(double(target_counts(:).'), 'stable'), 'descend');
end
