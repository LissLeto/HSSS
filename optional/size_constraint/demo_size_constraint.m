%% Experimental size-constrained HSSS demo

clear;
clc;
close all;

root_dir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root_dir);
setup_hsss();

image_id = '22090';
sample_dir = fullfile(root_dir, 'examples', 'BSDS500');
image_path = fullfile(sample_dir, [image_id, '.jpg']);
edge_path = fullfile(sample_dir, [image_id, '_rcf.png']);
rgb_image = imread(image_path);
target_count = 400;

options = struct();
options.EdgeMap = edge_path;
options.EdgePercentile = 80;
options.RandomSeed = 0;
options.SpatialWeight = 0.5;
options.SearchRadiusMultiplier = 5;
options.LambdaSize = 2;
options.SizeTau = 1.25;

[label_maps, info] = hsss_segment_size(rgb_image, target_count, options);
labels = label_maps(:, :, 1);
overlay = overlay_superpixel_boundaries(rgb_image, labels, [255, 0, 0]);

figure('Name', 'HSSS size-constraint demo');
subplot(1, 2, 1);
imshow(rgb_image);
title('Input image');
subplot(1, 2, 2);
imshow(overlay);
title(sprintf('Size-constrained HSSS: K = %d', max(labels(:))));

fprintf('Requested superpixels: %d\n', target_count);
fprintf('Returned superpixels:  %d\n', max(labels(:)));
fprintf('LambdaSize:             %.3g\n', info.LambdaSize);
fprintf('SizeTau:                %.3g\n', info.SizeTau);
fprintf('Used original fallback: %d\n', info.UsedOriginalFallback);
fprintf('Runtime:                %.3f s\n', info.RuntimeSeconds);
