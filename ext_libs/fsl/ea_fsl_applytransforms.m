function ea_fsl_applytransforms(varargin)
% Wrapper for FSL applywarp in terms of reapplying normalizations to
% pre- and postop imaging.

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

directory = [options.root,options.patientname,filesep];
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
    APPLYWARP = [basedir, 'applywarp.exe'];
else
    APPLYWARP = [basedir, 'applywarp.', computer('arch')];
end

[~,warpprefix] = fileparts(options.prefs.gprenii); % Prefix of the FNIRT warp field file

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
           refim = [ea_space,'t2.nii'];
        end

        if isempty(transformfile)
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

% generate l*.nii files
if nargin == 1 % standard case
    for fi = 1:length(ofis)
        try
            if ~exist(ofis{fi}, 'file')   % skip if normalized file doesn't exist
                fprintf('%s not found. Skip generating l*.nii files (small bounding box)...\n',ofis{fi});
                continue
            end

            matlabbatch{1}.spm.util.imcalc.input = {[ea_space,'bb.nii,1'];
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
