%% HSSS single-image demo
% Change image_id to select one of the bundled BSDS500 examples:
% 3063, 22090, 104055, 124084, or 228076.

clear;
clc;
close all;

root_dir = fileparts(mfilename('fullpath'));
setup_hsss();

image_id = '3063';
edge_mode = 'sobel';  % 'rcf' reads a saved map; 'sobel' computes it now.
run_evaluation_demo = true;  % Also show BR/IBR results for image 3063.
sample_dir = fullfile(root_dir, 'examples', 'BSDS500');
image_path = fullfile(sample_dir, [image_id, '.jpg']);
edge_path = fullfile(sample_dir, [image_id, '_rcf.png']);
rgb_image = imread(image_path);
target_count = [1200,1000,800];

options = struct();
switch lower(edge_mode)
    case 'rcf'
        if ~isfile(edge_path)
            error('HSSS:MissingRCFMap', ...
                'RCF edge map not found: %s', edge_path);
        end
        options.EdgeMap = edge_path;
        edge_name = 'RCF';
    case 'sobel'
        % Empty EdgeMap makes prepare_hsss_image compute a Sobel map from
        % the current RGB image by calling GetGar.
        options.EdgeMap = [];
        edge_name = 'Sobel';
    otherwise
        error('HSSS:InvalidEdgeMode', ...
            'edge_mode should be ''rcf'' or ''sobel''.');
end
options.EdgePercentile = 60;
options.RandomSeed = 0;
options.SpatialWeight = 0.5;
options.SearchRadiusMultiplier = 5;

[label_maps, info] = hsss_segment(rgb_image, target_count, options);

figure('Name', 'HSSS input image', 'NumberTitle', 'off');
imshow(rgb_image);
title('Input image');

% hsss_segment sorts requested counts from large to small. Each slice of
% label_maps corresponds to the count at the same position in
% info.TargetCounts. Display every hierarchy level in a separate window.
for idx = 1:size(label_maps, 3)
    labels = label_maps(:, :, idx);
    actual_count = info.ActualCounts(idx);
    overlay = overlay_superpixel_boundaries( ...
        rgb_image, labels, [255, 0, 0]);

    figure('Name', sprintf('HSSS %s K=%d', edge_name, actual_count), ...
        'NumberTitle', 'off');
    imshow(overlay);
    title(sprintf('HSSS (%s): K = %d', edge_name, actual_count));
end
fprintf('Edge mode:              %s\n', edge_name);
fprintf('Requested superpixels: %s\n', mat2str(target_count));
fprintf('Returned order:        %s\n', mat2str(info.ActualCounts));
fprintf('Foundation superpixels:     %d\n', info.InitialSuperpixels);
fprintf('Runtime:                %.3f s\n', info.RuntimeSeconds);

if run_evaluation_demo
    fprintf('\nRunning the BR/IBR evaluation example...\n');
    run(fullfile(sample_dir, 'demo_eval.m'));
end
