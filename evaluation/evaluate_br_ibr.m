function metrics = evaluate_br_ibr(labels, ground_truth)
%EVALUATE_BR_IBR Evaluate BR and IBR against one or more annotations.
%
% METRICS = EVALUATE_BR_IBR(LABELS, GROUND_TRUTH) accepts one segmentation
% matrix, a BSDS groundTruth cell array, or a structure loaded from a BSDS
% ground-truth MAT file. The returned table contains one row per annotation
% and a final mean row.

validateattributes(labels, {'numeric', 'logical'}, ...
    {'2d', 'nonempty', 'finite'}, mfilename, 'labels', 1);
segmentations = unpack_ground_truth(ground_truth);
num_gt = numel(segmentations);
br = zeros(num_gt, 1);
ibr = zeros(num_gt, 1);

for idx = 1:num_gt
    segmentation = segmentations{idx};
    if ~isequal(size(labels), size(segmentation))
        error('evaluate_br_ibr:SizeMismatch', ...
            'labels and ground truth %d have different sizes.', idx);
    end
    [br(idx), ibr(idx)] = eval_all(labels, segmentation);
end

ground_truth_index = [string((1:num_gt).'); "mean"];
metrics = table(ground_truth_index, [br; mean(br)], [ibr; mean(ibr)], ...
    'VariableNames', {'GroundTruth', 'BR', 'IBR'});
end

function segmentations = unpack_ground_truth(input_gt)
if isnumeric(input_gt) || islogical(input_gt)
    segmentations = {input_gt};
    return
end
if isstruct(input_gt) && isfield(input_gt, 'groundTruth')
    input_gt = input_gt.groundTruth;
end
if ~iscell(input_gt) || isempty(input_gt)
    error('evaluate_br_ibr:InvalidGroundTruth', ...
        'Ground truth must be a matrix, BSDS cell array, or loaded structure.');
end
segmentations = cell(size(input_gt));
for idx = 1:numel(input_gt)
    item = input_gt{idx};
    if isstruct(item) && isfield(item, 'Segmentation')
        segmentations{idx} = item.Segmentation;
    elseif isnumeric(item) || islogical(item)
        segmentations{idx} = item;
    else
        error('evaluate_br_ibr:InvalidGroundTruthItem', ...
            'Ground-truth item %d has no Segmentation field.', idx);
    end
end
end
