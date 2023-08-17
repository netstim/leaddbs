function ea_fsl_apply_coregistration(varargin)
% Wrapper to apply FSL flirt affine transformation

fixedimage = varargin{1};
movingimage = varargin{2};
outputimage = varargin{3};

volumedir = [fileparts(ea_niifileparts(movingimage)), filesep];

if nargin >= 4
    affine = varargin{4};
else
    % determine the affine matrix to be used
    [~, mov] = ea_niifileparts(movingimage);
    [~, fix] = ea_niifileparts(fixedimage);
    xfm = [mov, '2', fix, '_flirt'];
    affine = dir([volumedir, xfm, '*.mat']);

    if numel(affine) == 0
        error('Please run ea_flirt first before apply the transformation!');
    else
        affine = [volumedir, affine(end).name];
    end
end

% Available interpolation method:
% trilinear, nearestneighbour, sinc, spline
if nargin >= 5
    interp = varargin{5};
else
    interp = 'spline'; % spline interpolition by default
end

basedir = [fileparts(mfilename('fullpath')), filesep];
FLIRT = ea_getExec([basedir, 'flirt'], escapePath = 1);


cmd = [FLIRT, ...
       ' -in ', ea_path_helper(movingimage), ...
       ' -ref ', ea_path_helper(fixedimage), ...
       ' -out ', ea_path_helper(outputimage), ...
       ' -init ', ea_path_helper(affine), ...
       ' -applyxfm', ...
       ' -interp ', interp, ...
       ' -v'];

setenv('FSLOUTPUTTYPE','NIFTI');
ea_runcmd(cmd);
