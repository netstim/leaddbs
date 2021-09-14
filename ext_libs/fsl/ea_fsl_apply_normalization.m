function ea_fsl_apply_normalization(varargin)
% Wrapper for FSL applywarp used for reapplying normalization (non-linear
% registration)

options = varargin{1};

useinverse=0;
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
    interp=varargin{7};
    if ~ischar(interp)
        switch interp
            case 0
                interp='nn';
            otherwise
                interp='trilinear';

        end
    end
else
    % nn, trilinear, sinc, spline
    interp='trilinear';
end

if ~isempty(options) && ~isempty(fieldnames(options))
    directory = [options.root,options.patientname,filesep];
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
end

basedir = [fileparts(mfilename('fullpath')), filesep];
if ispc
    APPLYWARP = ea_path_helper([basedir, 'applywarp.exe']);
else
    APPLYWARP = [basedir, 'applywarp.', computer('arch')];
end

for fi = 1:length(fis)
    if ~exist(fis{fi}, 'file')   % skip if unnormalized file doesn't exist
        fprintf('%s not found. Skip normalization...\n',fis{fi});
        continue
    end

	cmd = [APPLYWARP, ...
           ' --verbose' ...
           ' --in=', ea_path_helper(fis{fi}), ...
           ' --out=', ea_path_helper(ofis{fi})];

    if useinverse
        if isempty(refim)
           refim = [directory,options.prefs.prenii_unnormalized];
        end

        if isempty(transformfile)
            [~,warpprefix] = fileparts(options.prefs.gprenii);
            cmd = [cmd, ...
                   ' --ref=', ea_path_helper(refim),...
                   ' --warp=',ea_path_helper([directory,warpprefix,'InverseWarpField.nii'])];
        else
            cmd = [cmd, ...
                   ' --ref=', ea_path_helper(refim), ...
                   ' --warp ', ea_path_helper(transformfile)];
        end
    else
        if isempty(refim)
            spacedef=ea_getspacedef;
            refim = [ea_space,spacedef.templates{1},'.nii'];
        end

        if isempty(transformfile)
            [~,warpprefix] = fileparts(options.prefs.gprenii);
            cmd = [cmd, ...
                   ' --ref=', ea_path_helper(refim),...
                   ' --warp=', ea_path_helper([directory,warpprefix,'WarpField.nii'])];
        else
        	cmd = [cmd, ...
                   ' --ref=', ea_path_helper(refim),...
                   ' --warp=', ea_path_helper(transformfile)];
        end
    end

    if ~isempty(interp)
        cmd = [cmd, ' --interp=', interp];
    end

    if strcmp(ofis{fi}(end-2:end),'.gz')
        setenv('FSLOUTPUTTYPE','NIFTI_GZ');
    else
        setenv('FSLOUTPUTTYPE','NIFTI');
    end

    if ~ispc
        system(['bash -c "', cmd, '"']);
    else
        system(cmd);
    end
end
