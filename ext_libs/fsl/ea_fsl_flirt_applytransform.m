function ea_fsl_flirt_applytransform(varargin)
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

basedir = [fileparts(mfilename('fullpath')), filesep];
if ispc
    APPLYWARP = [basedir, 'applywarp.exe'];
else
    APPLYWARP = [basedir, 'applywarp.', computer('arch')];
end

cmd = [APPLYWARP, ...
       ' -i ', ea_path_helper(movingimage), ...
       ' -r ', ea_path_helper(fixedimage), ...
       ' -o ', ea_path_helper(outputimage), ...
       ' --premat=', ea_path_helper(affine), ...
       ' --interp=spline', ...
       ' -v'];

setenv('FSLOUTPUTTYPE','NIFTI');
if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end
