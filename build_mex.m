function build_mex()
%BUILD_MEX Compile the platform-specific HSSS acceleration MEX files.

root_dir = fileparts(mfilename('fullpath'));
core_dir = fullfile(root_dir, 'core');
targets = { ...
    'assign_pixels_grid_mex', ...
    'connect_SP_singlepass_mex'};

fprintf('Building HSSS MEX files for %s...\n', mexext);
for target_index = 1:numel(targets)
    target = targets{target_index};
    source_file = fullfile(core_dir, [target, '.cpp']);
    if ~isfile(source_file)
        error('HSSS:MissingMexSource', ...
            'MEX source not found: %s', source_file);
    end
    mex('-R2018a', '-O', '-outdir', core_dir, ...
        '-output', target, source_file);
    fprintf('Built: %s\n', fullfile(core_dir, [target, '.', mexext]));
end
rehash;
end
