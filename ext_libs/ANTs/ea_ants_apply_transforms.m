function ea_ants_apply_transforms(varargin)
% Wrapper for antsApplyTransforms used for reapplying coregistration
% (linear registration) and normalization (non-linear registration)

ea_libs_helper;

options = varargin{1};

useinverse = 0;

if nargin > 1 % manual application
    input = varargin{2};
    output = varargin{3};
    if ischar(input)
        input = {input};
    end
    if ischar(output)
        output = {output};
    end
    useinverse = varargin{4};
end

if nargin >= 5
    ref = varargin{5};
else
    ref = ''; % use defaults
end

if nargin >= 6
    transform = varargin{6};
else
    transform = ''; % use defaults
end

% Linear, NearestNeighbor, MultiLabel, Gaussian, BSpline
% CosineWindowedSinc, WelchWindowedSinc, HammingWindowedSinc, LanczosWindowedSinc
% GenericLabel (Recommanded for label image)
if nargin >= 7
    interp = varargin{7};
    if ~ischar(interp)
        switch interp
            case 0
                interp = 'NearestNeighbor';
            case 1
                interp = 'Linear';
            case -1
                interp = 'GenericLabel';
            otherwise
                interp = 'BSpline';
        end
    end
else
    if useinverse
        interp = 'GenericLabel';
    else
        interp = 'LanczosWindowedSinc';
    end
end

if nargin >= 8
    imageDim = varargin{8};
    if isnumeric(imageDim)
        imageDim = num2str(imageDim);
    end
else
    imageDim = '3';
end

% Image type (-e 0/1/2/3/4): scalar/vector/tensor/time-series/multichannel
if nargin == 9
    imageType = varargin{9};
    if isnumeric(imageType)
        imageType = num2str(imageType);
    end
else
    imageType = '0';
end

if nargin == 1
    input{1} = options.subj.coreg.anat.preop.(options.subj.AnchorModality);
    output{1} = options.subj.norm.anat.preop.(options.subj.AnchorModality);

    json = loadjson(options.subj.norm.log.method);
    if contains(json.method, 'affine')
        % Three-step affine normalization (Schonecker 2009) used
        warpSuffix = 'ants.mat';
    else
        warpSuffix = 'ants.nii.gz';
    end

    switch options.modality
        case 1 % MR
            input = [input; struct2cell(options.subj.coreg.anat.postop)];
            output = [output; struct2cell(options.subj.norm.anat.postop)];
        case 2 % CT
            input = [input; options.subj.coreg.anat.postop.CT];
            output = [output; options.subj.norm.anat.postop.CT];
    end
end

basedir = [fileparts(mfilename('fullpath')), filesep];
if ispc
    applyTransforms = ea_path_helper([basedir, 'antsApplyTransforms.exe']);
else
    applyTransforms = [basedir, 'antsApplyTransforms.', computer('arch')];
end

for i = 1:length(input)
    if ~exist(input{i},'file')   % skip if unnormalized file doesn't exist
        fprintf('%s not found. Skip normalization...\n',input{i});
        continue
    end

	cmd = [applyTransforms, ...
           ' --verbose 1' ...
           ' --dimensionality ', imageDim, ...
           ' --input-image-type ', imageType, ...
           ' --float 1' ...
           ' --input ',ea_path_helper(input{i}), ...
           ' --output ',ea_path_helper(output{i})];

    if useinverse
        if isempty(ref)
            ref = options.subj.coreg.anat.preop.(options.subj.AnchorModality);
        end

        if isempty(transform)
            cmd = [cmd, ...
                   ' --reference-image ',ea_path_helper(ref),...
                   ' --transform [',ea_path_helper([options.subj.norm.transform.inverseBaseName, warpSuffix]),',0]'];
        else
            cmd = [cmd, ...
                   ' --reference-image ',ea_path_helper(ref),...
                   ' --transform [',ea_path_helper(transform),',0]'];
        end
    else
        if isempty(ref)
            spacedef = ea_getspacedef;
            ref = [ea_space, spacedef.templates{1}, '.nii'];
        end

        if isempty(transform)
            cmd = [cmd, ...
                   ' --reference-image ',ea_path_helper(ref),...
                   ' --transform [',ea_path_helper([options.subj.norm.transform.forwardBaseName, warpSuffix]),',0]'];
        else
            cmd = [cmd, ...
                   ' --reference-image ',ea_path_helper(ref),...
                   ' --transform [',ea_path_helper(transform),',0]'];
        end
    end

    if ~isempty(interp)
        cmd = [cmd, ' --interpolation ', interp];
    end

    [~, inputFileName] = ea_niifileparts(input{i});
    fprintf('\nNormalizing %s ...\n', inputFileName);

    if ~ispc
        system(['bash -c "', cmd, '"']);
    else
        system(cmd);
    end
end
