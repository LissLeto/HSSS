function setup_hsss()
%SETUP_HSSS Add the public HSSS source folders to the MATLAB path.

root_dir = fileparts(mfilename('fullpath'));
addpath(root_dir);
addpath(fullfile(root_dir, 'core'));
addpath(fullfile(root_dir, 'evaluation'));
addpath(fullfile(root_dir, 'optional', 'size_constraint'));
end
