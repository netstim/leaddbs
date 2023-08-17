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

if isempty(transform)
    json = loadjson(options.subj.norm.log.method);
    if isfield(json, 'custom') && json.custom
        % Custom full path of the transformation supplied.
        warpSuffix = '';
    elseif contains(json.method, 'affine', 'IgnoreCase', true)
        % Three-step affine normalization (Schonecker 2009) used
        warpSuffix = 'ants.mat';
    else
        warpSuffix = 'ants.nii.gz';
    end
    if useinverse
        transform = [options.subj.norm.transform.inverseBaseName, warpSuffix];
    else
        transform = [options.subj.norm.transform.forwardBaseName, warpSuffix];
    end
end

if nargin == 1
    input{1} = options.subj.coreg.anat.preop.(options.subj.AnchorModality);
    output{1} = options.subj.norm.anat.preop.(options.subj.AnchorModality);

    if strcmp(options.subj.postopModality, 'MRI')
        if exist(options.subj.brainshift.transform.scrf,'file') % apply brainshift correction to postop files on the fly.
            fn=fieldnames(options.subj.coreg.anat.postop);
            for postopfile=1:length(fn)
                uuid=ea_generate_uuid;
                copyfile(options.subj.coreg.anat.postop.(fn{postopfile}),[ea_getleadtempdir,uuid,'.nii']);
                nii=ea_load_nii([ea_getleadtempdir,uuid,'.nii']);
                scrf=load(options.subj.brainshift.transform.scrf);
                nii.mat=scrf.mat*nii.mat;
                ea_write_nii(nii);
                input = [input; [ea_getleadtempdir,uuid,'.nii']];
            end
        else
            input = [input; struct2cell(options.subj.coreg.anat.postop)];
        end
        output = [output; struct2cell(options.subj.norm.anat.postop)];
    elseif strcmp(options.subj.postopModality, 'CT')
        if exist(options.subj.brainshift.transform.scrf,'file') % apply brainshift correction to postop files on the fly.
            uuid=ea_generate_uuid;
            copyfile(options.subj.coreg.anat.postop.CT,[ea_getleadtempdir,uuid,'.nii']);
            nii=ea_load_nii([ea_getleadtempdir,uuid,'.nii']);
            scrf=load(options.subj.brainshift.transform.scrf);
            nii.mat=scrf.mat*nii.mat;
            ea_write_nii(nii);
            input = [input; [ea_getleadtempdir,uuid,'.nii']];
        else
            input = [input; options.subj.coreg.anat.postop.CT];
        end
        output = [output; options.subj.norm.anat.postop.CT];
    end
end

basedir = [fileparts(mfilename('fullpath')), filesep];
applyTransforms = ea_getExec([basedir, 'antsApplyTransforms'], escapePath = 1);


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

    if isempty(ref)
        if useinverse
            ref = options.subj.coreg.anat.preop.(options.subj.AnchorModality);
        else
            spacedef = ea_getspacedef;
            ref = [ea_space, spacedef.templates{1}, '.nii'];
        end
    end

    cmd = [cmd, ...
        ' --reference-image ',ea_path_helper(ref),...
        ' --transform [',ea_path_helper(transform),',0]'];

    if ~isempty(interp)
        cmd = [cmd, ' --interpolation ', interp];
    end

    [~, inputFileName] = ea_niifileparts(input{i});
    fprintf('\nNormalizing %s ...\n', inputFileName);

    ea_runcmd(cmd);
end
