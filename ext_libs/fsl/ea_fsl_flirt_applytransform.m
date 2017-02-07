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
    FLIRT = [basedir, 'flirt.exe'];
else
    FLIRT = [basedir, 'flirt.', computer('arch')];
end

cmd = [FLIRT, ...
       ' -in ', ea_path_helper(movingimage), ...
       ' -ref ', ea_path_helper(fixedimage), ...
       ' -out ', ea_path_helper(outputimage), ...
       ' -init ', ea_path_helper(affine), ...
       ' -applyxfm', ...
       ' -interp spline', ...
       ' -v'];

setenv('FSLOUTPUTTYPE','NIFTI');
if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end
