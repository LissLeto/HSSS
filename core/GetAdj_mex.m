function adjacency = GetAdj_mex(labels, superpixel_count)
%GETADJ_MEX Portable MATLAB fallback for the optional Windows MEX file.
%
% The returned symmetric sparse matrix stores shared 4-neighbor boundary
% length between adjacent labels. MATLAB gives a platform-compatible MEX
% file precedence over this .m implementation when one is available.

labels = double(labels);
if nargin < 2
    superpixel_count = max(labels(:));
end

left = labels(:, 1:end-1);
right = labels(:, 2:end);
horizontal = left ~= right;

top = labels(1:end-1, :);
bottom = labels(2:end, :);
vertical = top ~= bottom;

node_i = [left(horizontal); top(vertical)];
node_j = [right(horizontal); bottom(vertical)];
if isempty(node_i)
    adjacency = sparse(superpixel_count, superpixel_count);
    return
end

valid = node_i > 0 & node_j > 0;
node_i = node_i(valid);
node_j = node_j(valid);
adjacency = sparse([node_i; node_j], [node_j; node_i], 1, ...
    superpixel_count, superpixel_count);
end
