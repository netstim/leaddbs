function ea_ants_applytransforms(varargin)
% Wrapper for antsApplyTransforms in terms of reapplying normalizations to
% pre- and postop imaging.

ea_libs_helper;

options = varargin{1};

useinverse = 0;
if nargin > 1 % manual application
    fis = varargin{2};
    ofis = varargin{3};
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
    % Linear, NearestNeighbor, MultiLabel, Gaussian, BSpline
    % CosineWindowedSinc, WelchWindowedSinc, HammingWindowedSinc, LanczosWindowedSinc
    % GenericLabel (Recommanded for label image)
    if useinverse
        interp = 'GenericLabel';
    else
        interp = 'BSpline';
    end
end

directory = [options.root,options.patientname,filesep];
warpsuffix = ea_conv_antswarps(directory);

if nargin == 1
    switch options.modality
        case 1 % MR
            fis = {ea_niigz([directory,options.prefs.prenii_unnormalized])};
            ofis = {ea_niigz([directory,options.prefs.gprenii])};
            if isfield(options.prefs,'prenii')
                lfis = {ea_niigz([directory,options.prefs.prenii])};
            end

            if isfield(options.prefs,'tranii_unnormalized')
                fis = [fis,{ea_niigz([directory,options.prefs.tranii_unnormalized])}];
                ofis = [ofis,{ea_niigz([directory,options.prefs.gtranii])}];
                lfis = [lfis,{ea_niigz([directory,options.prefs.tranii])}];
            end

            if isfield(options.prefs,'cornii_unnormalized')
                fis = [fis,{ea_niigz([directory,options.prefs.cornii_unnormalized])}];
                ofis = [ofis,{ea_niigz([directory,options.prefs.gcornii])}];
                lfis = [lfis,{ea_niigz([directory,options.prefs.cornii])}];
            end

            if isfield(options.prefs,'sagnii_unnormalized')
                fis = [fis,{ea_niigz([directory,options.prefs.sagnii_unnormalized])}];
                ofis = [ofis,{ea_niigz([directory,options.prefs.gsagnii])}];
                lfis = [lfis,{ea_niigz([directory,options.prefs.sagnii])}];
            end
            if isfield(options.prefs,'fa2anat')
                if exist([directory,options.prefs.fa2anat],'file')
                    fis = [fis,{ea_niigz([directory,options.prefs.fa2anat])}];
                    ofis = [ofis,{ea_niigz([directory,'gl',options.prefs.fa2anat])}];
                    lfis = [lfis,{ea_niigz([directory,'l',options.prefs.fa2anat])}];
                end
            end
        case 2 % CT
            fis{1} = ea_niigz([directory,options.prefs.prenii_unnormalized]);
            fis{2} = ea_niigz([directory,options.prefs.ctnii_coregistered]);
            ofis{1} = ea_niigz([directory,options.prefs.gprenii]);
            ofis{2} = ea_niigz([directory,options.prefs.gctnii]);
            lfis{1} = ea_niigz([directory,options.prefs.prenii]);
            lfis{2} = ea_niigz([directory,options.prefs.ctnii]);
            if exist([directory,options.prefs.fa2anat],'file')
                fis{3} = ea_niigz([directory,options.prefs.fa2anat]);
                ofis{3} = ea_niigz([directory,'gl',options.prefs.fa2anat]);
                lfis{3} = ea_niigz([directory,'l',options.prefs.fa2anat]);
            end
    end

            [fis,ofis] = ea_appendgrid(options,fis,ofis,1);
end

basedir = [fileparts(mfilename('fullpath')), filesep];
if ispc
    applyTransforms = [basedir, 'antsApplyTransforms.exe'];
else
    applyTransforms = [basedir, 'antsApplyTransforms.', computer('arch')];
end

for fi = 1:length(fis)
    if ~exist(fis{fi},'file')   % skip if unnormalized file doesn't exist
        fprintf('%s not found. Skip normalization...\n',fis{fi});
        continue
    end

    % generate gl*.nii files
    [~,glprebase] = fileparts(options.prefs.gprenii);

	cmd = [applyTransforms, ...
           ' --verbose 1' ...
           ' --dimensionality 3 --float 1' ...
           ' --input ',ea_path_helper(fis{fi}), ...
           ' --output ',ea_path_helper(ofis{fi})];

    if useinverse
        if isempty(refim)
            refim = [directory,options.prefs.prenii_unnormalized];
        end

        if isempty(transformfile)
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

% generate l*.nii files
if nargin == 1 % standard case
    for fi = 1:length(ofis)
        try
            if ~exist(ofis{fi}, 'file')   % skip if normalized file doesn't exist
                fprintf('%s not found. Skip generating l*.nii files (small bounding box)...\n',ofis{fi});
                continue
            end

            matlabbatch{1}.spm.util.imcalc.input = {[ea_space(options),'bb.nii,1'];
                                                    [ofis{fi},',1']};
            matlabbatch{1}.spm.util.imcalc.output = lfis{fi};
            matlabbatch{1}.spm.util.imcalc.outdir = {directory};
            matlabbatch{1}.spm.util.imcalc.expression = 'i2';
            matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

            try spm_jobman('run',{matlabbatch}); end
            clear matlabbatch
        end
    end
end
