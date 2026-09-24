%% Run and evaluate HSSS on BSDS500 image 3063

example_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(example_dir));
addpath(root_dir);
setup_hsss();

image_id = '3063';
rgb_image = imread(fullfile(example_dir, [image_id, '.jpg']));

options = struct();
options.EdgeMap = [];
options.EdgePercentile = 60;
options.RandomSeed = 0;
options.SpatialWeight = 0.5;
options.SearchRadiusMultiplier = 5;
[label_maps, info] = hsss_segment(rgb_image, 400, options);
labels = label_maps(:, :, 1);

loaded_gt = load(fullfile(example_dir, [image_id, '_gt.mat']));
metrics = evaluate_br_ibr(labels, loaded_gt);

fprintf('Superpixel number: %d\n', info.ActualCounts(1));
disp(metrics);

% Metrics above use every annotation. The visualization uses annotation 4.
eval_all(labels, loaded_gt.groundTruth{4}.Segmentation, rgb_image, true);
