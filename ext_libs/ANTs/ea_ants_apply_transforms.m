function ea_ants_apply_transforms(varargin)
% Wrapper for antsApplyTransforms used for reapplying coregistration
% (linear registration) and normalization (non-linear registration)

ea_libs_helper;

options = varargin{1};

useinverse = 0;

if nargin > 1 % manual application
    fis = varargin{2};
    ofis = varargin{3};
    if ischar(fis)
        fis = {fis};
    end
    if ischar(ofis)
        ofis = {ofis};
    end
    useinverse = varargin{4};
end

if nargin >= 5
    refim = varargin{5};
else
    refim = ''; % use defaults
end

if nargin >= 6
    transformfile = varargin{6};
else
    transformfile = ''; % use defaults
end

% Linear, NearestNeighbor, MultiLabel, Gaussian, BSpline
% CosineWindowedSinc, WelchWindowedSinc, HammingWindowedSinc, LanczosWindowedSinc
% GenericLabel (Recommanded for label image)
if nargin >= 7
    interp = varargin{7};
    if ~ischar(interp)
        switch interp
            case 0
                interp='NearestNeighbor';
            case 1
                interp='Linear';
            case -1
                interp='GenericLabel';
            otherwise
                interp='BSpline';
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

if ~isempty(options) && ~isempty(fieldnames(options))
    directory = [options.root,options.patientname,filesep];
    warpsuffix = ea_getantstransformext(directory);
    [~,anatpresent]=ea_assignpretra(options);
end
if nargin == 1
    switch options.modality
        case 1 % MR
            fis = {ea_niigz([directory,options.prefs.prenii_unnormalized])};
            ofis = {ea_niigz([directory,options.prefs.gprenii])};


            if isfield(options.prefs,'tranii_unnormalized')
                fis = [fis,{ea_niigz([directory,options.prefs.tranii_unnormalized])}];
                ofis = [ofis,{ea_niigz([directory,options.prefs.gtranii])}];
            end

            if isfield(options.prefs,'cornii_unnormalized')
                fis = [fis,{ea_niigz([directory,options.prefs.cornii_unnormalized])}];
                ofis = [ofis,{ea_niigz([directory,options.prefs.gcornii])}];
            end

            if isfield(options.prefs,'sagnii_unnormalized')
                fis = [fis,{ea_niigz([directory,options.prefs.sagnii_unnormalized])}];
                ofis = [ofis,{ea_niigz([directory,options.prefs.gsagnii])}];
            end
            if isfield(options.prefs,'fa2anat')
                if exist([directory,options.prefs.fa2anat],'file')
                    fis = [fis,{ea_niigz([directory,options.prefs.fa2anat])}];
                    ofis = [ofis,{ea_niigz([directory,'gl',options.prefs.fa2anat])}];
                end
            end
        case 2 % CT
            fis{1} = ea_niigz([directory,options.prefs.prenii_unnormalized]);
            fis{2} = ea_niigz([directory,options.prefs.ctnii_coregistered]);
            ofis{1} = ea_niigz([directory,options.prefs.gprenii]);
            ofis{2} = ea_niigz([directory,options.prefs.gctnii]);
            if exist([directory,options.prefs.fa2anat],'file')
                fis{3} = ea_niigz([directory,options.prefs.fa2anat]);
                ofis{3} = ea_niigz([directory,'gl',options.prefs.fa2anat]);
            end
    end

    [fis,ofis] = ea_appendgrid(options,fis,ofis,1);
end

basedir = [fileparts(mfilename('fullpath')), filesep];
if ispc
    applyTransforms = ea_path_helper([basedir, 'antsApplyTransforms.exe']);
else
    applyTransforms = [basedir, 'antsApplyTransforms.', computer('arch')];
end

for fi = 1:length(fis)
    if ~exist(fis{fi},'file')   % skip if unnormalized file doesn't exist
        fprintf('%s not found. Skip normalization...\n',fis{fi});
        continue
    end

	cmd = [applyTransforms, ...
           ' --verbose 1' ...
           ' --dimensionality ', imageDim, ...
           ' --input-image-type ', imageType, ...
           ' --float 1' ...
           ' --input ',ea_path_helper(fis{fi}), ...
           ' --output ',ea_path_helper(ofis{fi})];

    if useinverse
        if isempty(refim)
            refim = [directory,anatpresent{1}];
        end

        if isempty(transformfile)
            [~,glprebase] = fileparts(options.prefs.gprenii);
            cmd = [cmd, ...
                   ' --reference-image ',ea_path_helper(refim),...
                   ' --transform [',ea_path_helper([directory,glprebase,'InverseComposite',warpsuffix]),',0]'];
        else
            cmd = [cmd, ...
                   ' --reference-image ',ea_path_helper(refim),...
                   ' --transform [',ea_path_helper(transformfile),',0]'];
        end
    else
        if isempty(refim)
            spacedef=ea_getspacedef;
            refim = [ea_space,spacedef.templates{1},'.nii'];
        end

        if isempty(transformfile)
            [~,glprebase] = fileparts(options.prefs.gprenii);
            cmd = [cmd, ...
                   ' --reference-image ',ea_path_helper(refim),...
                   ' --transform [',ea_path_helper([directory,glprebase,'Composite',warpsuffix]),',0]'];
        else
            cmd = [cmd, ...
                   ' --reference-image ',ea_path_helper(refim),...
                   ' --transform [',ea_path_helper(transformfile),',0]'];
        end
    end

    if ~isempty(interp)
        cmd = [cmd, ' --interpolation ', interp];
    end
    if ~ispc
        system(['bash -c "', cmd, '"']);
    else
        system(cmd);
    end
end
