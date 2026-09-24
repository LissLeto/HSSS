function output = overlay_superpixel_boundaries(rgb_image, labels, color)
%OVERLAY_SUPERPIXEL_BOUNDARIES Draw superpixel boundaries on an image.

if nargin < 3
    color = [255, 0, 0];
end
if ischar(rgb_image) || isstring(rgb_image)
    rgb_image = imread(rgb_image);
end
if ndims(rgb_image) == 2
    rgb_image = repmat(rgb_image, 1, 1, 3);
end

output = im2uint8(rgb_image);
boundary = boundarymask(labels, 4);
for channel = 1:3
    layer = output(:, :, channel);
    layer(boundary) = color(channel);
    output(:, :, channel) = layer;
end
end
