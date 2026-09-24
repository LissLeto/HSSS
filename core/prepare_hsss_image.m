function [labels, lab_image, para, seeds, edge_map, info] = ...
    prepare_hsss_image(rgb_image, options)
%PREPARE_HSSS_IMAGE Build connected foundation superpixels for HSSS.

% This is a thin public API around the original foundation implementation.

if nargin < 2
    options = struct();
end
options = apply_defaults(options);

if ischar(rgb_image) || isstring(rgb_image)
    rgb_image = imread(rgb_image);
end
if ndims(rgb_image) == 2
    rgb_image = repmat(rgb_image, 1, 1, 3);
end
if size(rgb_image, 3) ~= 3
    error('HSSS:InvalidImage', 'Input must be an RGB or grayscale image.');
end

lab = rgb2lab(rgb_image);
[L, A, B] = imsplit(lab);
L = normalize_channel(L, false);
A = normalize_channel(A, true);
B = normalize_channel(B, true);
lab_image = cat(3, L, A, B);
para = Getpara2(lab_image, options.MaxL);

if isempty(options.EdgeMap)
    [edge_map, ~, ~] = GetGar(rgb_image);
else
    edge_map = options.EdgeMap;
    if ischar(edge_map) || isstring(edge_map)
        edge_map = imread(edge_map);
    end
    if ndims(edge_map) == 3
        edge_map = rgb2gray(edge_map);
    end
end
if ~isequal(size(edge_map, 1), size(rgb_image, 1)) || ...
        ~isequal(size(edge_map, 2), size(rgb_image, 2))
    error('HSSS:EdgeSizeMismatch', ...
        'The edge map and input image must have the same height and width.');
end

edge_threshold = prctile(double(edge_map(:)), options.EdgePercentile);
rng(options.RandomSeed, 'twister');
[labels, ~, seeds] = fundationSP_accelerated( ...
    lab_image, rgb_image, edge_map, edge_threshold, options.MaxL, ...
    options.SpatialWeight, options.SearchRadiusMultiplier);
labels = connect_SP_singlepass_mex(labels, lab_image, para);
labels = compact_labels(labels);

info = struct();
info.InitialSuperpixels = max(labels(:));
info.SeedCount = nnz(seeds);
info.EdgeThreshold = edge_threshold;
info.EdgePercentile = options.EdgePercentile;
info.SpatialWeight = options.SpatialWeight;
info.SearchRadiusMultiplier = options.SearchRadiusMultiplier;
info.Parameters = para;
end

function options = apply_defaults(options)
defaults = struct('EdgeMap', [], 'EdgePercentile', 80, ...
    'MaxL', 1, 'RandomSeed', 0, 'SpatialWeight', 0.5, ...
    'SearchRadiusMultiplier', 5);
names = fieldnames(defaults);
for idx = 1:numel(names)
    if ~isfield(options, names{idx})
        options.(names{idx}) = defaults.(names{idx});
    end
end
validateattributes(options.SpatialWeight, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'nonnegative'});
validateattributes(options.SearchRadiusMultiplier, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});
end

function out = normalize_channel(channel, subtract_minimum)
channel = double(channel);
if subtract_minimum
    channel = channel - min(channel(:));
end
denominator = max(channel(:));
if denominator == 0
    out = zeros(size(channel));
else
    out = channel / denominator;
end
end

function labels_out = compact_labels(labels_in)
[~, ~, compact] = unique(labels_in(:), 'sorted');
labels_out = reshape(compact, size(labels_in));
labels_out = double(labels_out);
end
