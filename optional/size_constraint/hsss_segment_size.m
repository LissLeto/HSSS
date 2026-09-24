function [label_maps, info] = hsss_segment_size(rgb_image, target_counts, options)
%HSSS_SEGMENT_SIZE HSSS with an optional soft superpixel-size constraint.
%
% [LABEL_MAPS, INFO] = HSSS_SEGMENT_SIZE(IMAGE, TARGET_COUNTS, OPTIONS)
% uses the same foundation segmentation as HSSS_SEGMENT. The merge score
% receives a soft penalty when a candidate region would exceed SizeTau
% times the target average area.
%
% Additional OPTIONS fields:
%   LambdaSize       Size-penalty weight (default 2).
%   SizeTau          Free size ratio before penalization (default 1.25).

% This is an experimental add-on. It does not replace the original method.

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
used_original_fallback = false;
while ~isempty(requested) && rounds < options.MaxMergeRounds
    rounds = rounds + 1;
    labels_before = labels;
    requested_before = requested;
    [labels, label_maps, requested] = series_SP_fast_size( ...
        labels, lab_image, requested, para, options.SeriesLevel, ...
        label_maps, options.DGap, options.LambdaSize, options.SizeTau);
    if isequal(labels, labels_before) && isequal(requested, requested_before)
        % Strong penalties can leave no admissible size-constrained merge
        % for some edge maps. Finish the requested hierarchy with the
        % original rule instead of returning an inexact superpixel count.
        used_original_fallback = true;
        [labels, label_maps, requested, para] = series_SP_fast( ...
            labels, lab_image, requested, para, options.SeriesLevel, ...
            label_maps, options.DGap);
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
info.Method = "size_constraint";
info.LambdaSize = options.LambdaSize;
info.SizeTau = options.SizeTau;
info.UsedOriginalFallback = used_original_fallback;
end

function options = apply_defaults(options)
defaults = struct('EdgeMap', [], 'EdgePercentile', 80, 'MaxL', 1, ...
    'DGap', 1, 'RandomSeed', 0, 'MaxMergeRounds', 30, ...
    'SeriesLevel', 1, 'LambdaSize', 2, 'SizeTau', 1.25, ...
    'SpatialWeight', 0.5, 'SearchRadiusMultiplier', 5);
names = fieldnames(defaults);
for idx = 1:numel(names)
    if ~isfield(options, names{idx})
        options.(names{idx}) = defaults.(names{idx});
    end
end
validateattributes(options.LambdaSize, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'nonnegative'});
validateattributes(options.SizeTau, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});
end

function targets = validate_targets(target_counts)
validateattributes(target_counts, {'numeric'}, ...
    {'vector', 'integer', 'positive', 'finite'});
targets = sort(unique(double(target_counts(:).'), 'stable'), 'descend');
end
