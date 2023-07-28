function ea_fsl_apply_normalization(varargin)
% Wrapper for FSL applywarp used for reapplying normalization (non-linear
% registration)

options = varargin{1};

useinverse = 0;

if nargin > 1 % manual application
    input = varargin{2};
    output = varargin{3};
    useinverse = varargin{4};
end

if nargin >= 5
    ref = varargin{5};
else
    ref = '';
end

if nargin >= 6
    transform = varargin{6};
else
    transform = '';
end

if nargin >= 7
    interp = varargin{7};
    if ~ischar(interp)
        switch interp
            case 0
                interp = 'nn';
            otherwise
                interp = 'trilinear';

        end
    end
else
    % nn, trilinear, sinc, spline
    interp = 'trilinear';
end

if nargin == 1
    input{1} = options.subj.coreg.anat.preop.(options.subj.AnchorModality);
    output{1} = options.subj.norm.anat.preop.(options.subj.AnchorModality);

    switch options.subj.postopModality
        case 'MRI'
            input = [input; struct2cell(options.subj.coreg.anat.postop)];
            output = [output; struct2cell(options.subj.norm.anat.postop)];
        case 'CT'
            input = [input; options.subj.coreg.anat.postop.CT];
            output = [output; options.subj.norm.anat.postop.CT];
    end
end

basedir = [fileparts(mfilename('fullpath')), filesep];
APPLYWARP = ea_getExec([basedir, 'applywarp'], escapePath = 1);


for i = 1:length(input)
    if isfile(input{i})
        cmd = [APPLYWARP, ...
               ' --verbose' ...
               ' --in=', ea_path_helper(input{i}), ...
               ' --out=', ea_path_helper(output{i})];

        if useinverse
            if isempty(ref)
               ref = options.subj.coreg.anat.preop.(options.subj.AnchorModality);
            end

            if isempty(transform)
                cmd = [cmd, ...
                       ' --ref=', ea_path_helper(ref),...
                       ' --warp=', ea_path_helper([options.subj.norm.transform.inverseBaseName, 'fnirt'])];
            else
                cmd = [cmd, ...
                       ' --ref=', ea_path_helper(ref), ...
                       ' --warp ', ea_path_helper(transform)];
            end
        else
            if isempty(ref)
                spacedef = ea_getspacedef;
                ref = [ea_space, spacedef.templates{1}, '.nii'];
            end

            if isempty(transform)
                cmd = [cmd, ...
                       ' --ref=', ea_path_helper(ref),...
                       ' --warp=', ea_path_helper([options.subj.norm.transform.forwardBaseName, 'fnirt'])];
            else
                cmd = [cmd, ...
                       ' --ref=', ea_path_helper(ref),...
                       ' --warp=', ea_path_helper(transform)];
            end
        end

        if ~isempty(interp)
            cmd = [cmd, ' --interp=', interp];
        end

        if strcmp(output{i}(end-2:end),'.gz')
            setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
        else
            setenv('FSLOUTPUTTYPE', 'NIFTI');
        end

        [~, inputFileName] = ea_niifileparts(input{i});
        fprintf('\nNormalizing %s ...\n', inputFileName);

        ea_runcmd(cmd);
    end
end
