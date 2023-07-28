function ea_fsl_reslice(input, reference, output, interp, verbose)
% Reslice images which are in the same space but with different resolutions
% using flirt.

% Overwrite the input image if no output name specified.
if ~exist('output', 'var')
    output = input;
end

% Set default interp method to trilinear
% Can choose from trilinear, nearestneighbour, sinc or spline
if ~exist('interp', 'var')
    interp = 'trilinear';
end

% Verbose by default
if ~exist('verbose', 'var')
    verbose = 1;
end

basedir = [fileparts(mfilename('fullpath')), filesep];
FLIRT = ea_getExec([basedir, 'flirt'], escapePath = 1);


cmd = [FLIRT, ...
    ' -in ', ea_path_helper(input), ...
    ' -ref ', ea_path_helper(reference), ...
    ' -applyxfm -usesqform -interp ', interp, ...
    ' -out ', output];

if verbose
    fprintf('Reslicing %s to match %s ...\n', input, reference);
end

setenv('FSLOUTPUTTYPE','NIFTI');
ea_runcmd(cmd);
