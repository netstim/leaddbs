function ea_lc_store_coreg(options, moving, fixed, output, transforms)
% Publish exactly one current forward/inverse B0 transform set.
% Called only after coregistration and output reslicing have succeeded.
directory = fullfile(options.root, options.patientname);
folder = fullfile(directory, 'coregistration', 'transformations');
protocol = fullfile(directory, 'coregistration', 'log', 'ea_lc_b0_coreg.mat');
if numel(transforms) ~= 2 || ~all(cellfun(@isfile, transforms)) || ~isfile(output)
    error('LeadDBS:IncompleteB0Coreg', 'B0 registration did not return both transforms and its output image.');
end
method = lower(options.coregmr.method);
if contains(method, 'nonlinear')
    record.method = 'ANTs';
    record.nonlinear = true;
    suffix = 'ants.nii.gz';
elseif contains(method, 'hybrid')
    error('LeadDBS:UnsupportedB0Coreg', 'Hybrid B0 registration requires a composed transform.');
elseif contains(method, 'ants')
    record.method = 'ANTs';
    record.nonlinear = false;
    suffix = 'ants.mat';
elseif contains(method, 'spm')
    record.method = 'SPM';
    record.nonlinear = false;
    suffix = 'spm.mat';
elseif contains(method, 'flirt')
    record.method = 'FLIRT';
    record.nonlinear = false;
    suffix = 'flirt.mat';
else
    error('LeadDBS:UnsupportedB0Coreg', 'Unsupported B0 registration method: %s', options.coregmr.method);
end
ea_mkdir(folder);
ea_mkdir(fileparts(protocol));
record.forward = fullfile('coregistration', 'transformations', ...
    [options.patientname, '_from-b0_to-anchorNative_desc-', suffix]);
record.inverse = fullfile('coregistration', 'transformations', ...
    [options.patientname, '_from-anchorNative_to-b0_desc-', suffix]);
record.moving = relative(moving, directory);
record.fixed = relative(fixed, directory);
record.output = relative(output, directory);
[~, anchorRel] = ea_lc_resolve_anat_anchor(options);
if isempty(anchorRel)
    error('LeadDBS:MissingB0Anchor', 'Cannot record B0 registration without its anatomical anchor.');
end
record.anchor = anchorRel;
record.registrationMethod = options.coregmr.method;
previous = [];
if isfile(protocol)
    previous = load(protocol, 'record');
end
copyfile(transforms{1}, fullfile(directory, record.forward));
copyfile(transforms{2}, fullfile(directory, record.inverse));
save(protocol, 'record');
% Delete only explicitly recorded superseded files after publishing success.
if ~isempty(previous)
    for field = {'forward', 'inverse'}
        old = fullfile(directory, previous.record.(field{1}));
        if ~strcmp(old, fullfile(directory, record.forward)) && ...
                ~strcmp(old, fullfile(directory, record.inverse)) && isfile(old)
            ea_delete(old);
        end
    end
end
for i = 1:numel(transforms)
    if ~strcmp(transforms{i}, fullfile(directory, record.forward)) && ...
            ~strcmp(transforms{i}, fullfile(directory, record.inverse))
        ea_delete(transforms{i});
    end
end
% Remove exact legacy filenames for this pair, never use them for selection.
[~, mov] = ea_niifileparts(moving);
[~, fix] = ea_niifileparts(fixed);
prefixes = {[mov, '2', fix], [fix, '2', mov]};
suffixes = {'_spm.mat', '_ants1.mat', '_ants.mat', '_flirt.mat', ...
    'Composite.nii.gz', 'InverseComposite.nii.gz'};
folders = {fileparts(moving), fileparts(fixed), folder, directory};
for p = 1:numel(prefixes)
    for s = 1:numel(suffixes)
        for d = 1:numel(folders)
            old = fullfile(folders{d}, [prefixes{p}, suffixes{s}]);
            if isfile(old)
                ea_delete(old);
            end
        end
    end
end
fprintf('Current B0 registration recorded: %s\n', record.forward);
end

function path = relative(path, directory)
prefix = [directory, filesep];
if ~startsWith(path, prefix)
    error('LeadDBS:InvalidB0CoregPath', 'B0 registration files must belong to the subject directory.');
end
path = path(numel(prefix)+1:end);
end
