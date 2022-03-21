function ea_fnirt(varargin)
% Wrapper for FSL FNIRT nonlinear registration

fixedimage = varargin{1};
movingimage = varargin{2};
outputimage = varargin{3};

% Set output dir
volumedir = [fileparts(ea_niifileparts(outputimage)), filesep];

% Do linear registration first, generate the affine matrix 'fslaffine*.mat'
movPath = ea_niifileparts(movingimage);
movingimage_flirt = [movPath, '_flirt'];

if isempty(dir([movingimage_flirt,'.nii*']))
    affineInit = ea_flirt(fixedimage, movingimage, movingimage_flirt);

    % Move initial affine transformation to 'normalization/transformations'
    % folder in case BIDS dataset.
    if isBIDSFileName(movingimage)
        parsedStruct = parseBIDSFilePath(movingimage);
        transformDir = [volumedir, '..', filesep, 'transformations', filesep];
        ea_mkdir(transformDir);
        movefile(affineInit{1}, [transformDir, 'sub-', parsedStruct.sub, '_from-anchorNative', '_to-', ea_getspace, '_desc-flirt.mat']);
        movefile(affineInit{2}, [transformDir, 'sub-', parsedStruct.sub, '_from-', ea_getspace, '_to-anchorNative', '_desc-flirt.mat']);

        % Set init affine transformation
        affineInit = [transformDir, 'sub-', parsedStruct.sub, '_from-anchorNative', '_to-', ea_getspace, '_desc-flirt.mat'];
    end
end

% Clean up coregistered image (only transformation needed)
ea_delete([movingimage_flirt, '.nii*']);

umachine = load([ea_gethome, '.ea_prefs.mat']);
normsettings = umachine.machine.normsettings;
if normsettings.fsl_skullstrip % skullstripping is on
    if isBIDSFileName(movingimage)
        parsedStruct = parseBIDSFilePath(movingimage);
        movingimage_brainmask = setBIDSEntity(movingimage, 'label', 'Brain', 'mod', parsedStruct.suffix, 'suffix', 'mask');
    else
        movingimage_brainmask = [movPath, '_brainmask'];
    end

    fixedPath = ea_niifileparts(fixedimage);
    if isBIDSFileName(fixedimage)
        parsedStruct = parseBIDSFilePath(fixedimage);
        fixedimage_brainmask = setBIDSEntity(fixedimage, 'label', 'Brain', 'mod', parsedStruct.suffix, 'suffix', 'mask');
    else
        fixedimage_brainmask = [fixedPath, '_brainmask'];
    end

    mask_params = [' --refmask=', ea_path_helper(fixedimage_brainmask), ...
                  ' --inmask=', ea_path_helper(movingimage_brainmask)];
else % skullstripping is off
    mask_params = '';
end

% Additional inital affine transformation provided
if nargin >= 4
    affineInit = varargin{4};
end

fprintf('\n\nRunning FSL FNIRT: %s\n\n', movingimage);

[~, warpprefix] = ea_niifileparts(outputimage);

% template config
% fnirtstage = [' --subsamp=4,4,2,2,1,1' ...
%               ' --warpres=8,8,8' ...
%               ' --miter=5,5,5,5,5,10' ...
%               ' --infwhm=8,6,5,4.5,3,2' ...
%               ' --reffwhm=8,6,5,4,2,0' ...
%               ' --lambda=300,150,100,50,40,30' ...
%               ' --estint=1,1,1,1,1,0' ...
%               ' --ssqlambda=1' ...
%               ' --regmod=bending_energy' ...
%               ' --intmod=global_non_linear_with_bias' ...
%               ' --intorder=5' ...
%               ' --biasres=50,50,50' ...
%               ' --biaslambda=10000' ...
%               ' --refderiv=0' ...
%               ' --applyrefmask=1,1,1,1,1,1' ...
%               ' --applyinmask=1'];

fnirtstage = [' --ref=', ea_path_helper(fixedimage), ...
              ' --in=', ea_path_helper(movingimage), ...
              mask_params, ...
              ' --aff=', ea_path_helper(affineInit), ...
              ' --iout=', ea_path_helper(outputimage), ...
              ' --cout=', ea_path_helper([volumedir, warpprefix, 'WarpCoef.nii']), ...
              ' --fout=', ea_path_helper([volumedir, warpprefix, 'WarpField.nii']), ...
              ' --subsamp=4,2' ...
              ' --warpres=8,8,8' ...
              ' --miter=5,10' ...
              ' --infwhm=8,2' ...
              ' --reffwhm=4,0' ...
              ' --lambda=300,30' ...
              ' --estint=1,0' ...
              ' --ssqlambda=1' ...
              ' --regmod=bending_energy' ...
              ' --intmod=global_non_linear_with_bias' ...
              ' --intorder=5' ...
              ' --biasres=50,50,50' ...
              ' --biaslambda=10000' ...
              ' --refderiv=0' ...
              ' --applyrefmask=1,1' ...
              ' --applyinmask=1,1' ...
              ' --verbose'];

invwarpstage = [' --warp=', ea_path_helper([volumedir, warpprefix, 'WarpField.nii']), ...
                ' --out=', ea_path_helper([volumedir, warpprefix, 'InverseWarpField.nii']), ...
                ' --ref=', ea_path_helper(movingimage), ...
                ' --verbose'];

basedir = [fileparts(mfilename('fullpath')), filesep];
if ispc
    FNIRT = ea_path_helper([basedir, 'fnirt.exe']);
    INVWARP = ea_path_helper([basedir, 'invwarp.exe']);
else
    FNIRT = [basedir, 'fnirt.', computer('arch')];
    INVWARP = [basedir, 'invwarp.', computer('arch')];
end

fnirtcmd = [FNIRT, fnirtstage];
invwarpcmd = [INVWARP, invwarpstage];

setenv('FSLOUTPUTTYPE','NIFTI');
if ~ispc
    system(['bash -c "', fnirtcmd, '"']);
    system(['bash -c "', invwarpcmd, '"']);
else
    system(fnirtcmd);
    system(invwarpcmd);
end

% Clean up waro coef file
ea_delete([volumedir, warpprefix, 'WarpCoef.nii'])

fprintf('\nFSL FNIRT finished\n');
