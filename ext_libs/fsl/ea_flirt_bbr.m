function affinefile = ea_flirt_bbr(varargin)
% Wrapper for FSL BBR linear registration

fixedimage = varargin{1};
movingimage = varargin{2};
outputimage = varargin{3};

if nargin >= 4
    writeoutmat = varargin{4};
else
    writeoutmat = 1;
end

if nargin >= 5
    if isempty(varargin{5}) % [] or {} or ''
        otherfiles = {};
    elseif ischar(varargin{5}) % single file, make it to cell string
        otherfiles = varargin(5);
    else % cell string
        otherfiles = varargin{5};
    end
else
    otherfiles = {};
end

% Prepare bet image for flirt bbr
[fixpath, fixname] = ea_niifileparts(fixedimage);
fixedimage_bet = [fileparts(fixpath), filesep, fixname, '_brain'];
if isempty(dir([fixedimage_bet,'.nii*']))
    ea_bet(fixedimage, 0, fixedimage_bet, 0.5);
end

volumedir = [fileparts(ea_niifileparts(outputimage)), filesep];

% name of the output transformation
[~, movname] = ea_niifileparts(movingimage);
xfm = [movname, '2', fixname, '_flirtbbr'];

basedir = [fileparts(mfilename('fullpath')), filesep];
if ispc
    FAST = ea_path_helper([basedir, 'fast.exe']);
    FSLMATHS = ea_path_helper([basedir, 'fslmaths.exe']);
    FLIRT = ea_path_helper([basedir, 'flirt.exe']);
    COVERT_XFM = ea_path_helper([basedir, 'convert_xfm.exe']);
else
    FAST = [basedir, 'fast.', computer('arch')];
    FSLMATHS = [basedir, 'fslmaths.', computer('arch')];
    FLIRT = [basedir, 'flirt.', computer('arch')];
    COVERT_XFM = [basedir, 'convert_xfm.', computer('arch')];
end

fixedimage_fast = [fileparts(fixpath), filesep, fixname, '_fslfast'];

% Segment the reference image
fastcmd = [FAST, ' -v ', ...
           '-o ', ea_path_helper(fixedimage_fast), ' ', ...
           '-A ', ea_path_helper([basedir, 'data/standard/tissuepriors/avg152T1_csf']), ' ', ...
                  ea_path_helper([basedir, 'data/standard/tissuepriors/avg152T1_gray']), ' ', ...
                  ea_path_helper([basedir, 'data/standard/tissuepriors/avg152T1_white']), ' ', ...
           ea_path_helper(fixedimage_bet)];

% Creat white matter mask
wmsegcmd = [FSLMATHS, ...
            ' ', ea_path_helper(fixedimage_fast), '_pve_2', ...
            ' -thr 0.5 -bin', ...
            ' ', ea_path_helper(fixedimage_fast), '_wmseg'];


% Creat white matter edge map ???just for visualization???
wmedgecmd = [FSLMATHS, ...
             ' ', ea_path_helper(fixedimage_fast), '_wmseg', ...
             ' -edge -bin', ...
             ' -mas ', ea_path_helper(fixedimage_fast), '_wmseg', ...
             ' ', ea_path_helper(fixedimage_fast), '_wmedge'];

% Standard FLIRT pre-alignment
flirtinitcmd = [FLIRT, ...
               ' -ref ', ea_path_helper(fixedimage_bet), ...
               ' -in ', ea_path_helper(movingimage), ...
               ' -omat ', ea_path_helper(volumedir), xfm, '_init.mat', ...
               ' -dof 6', ...
               ' -v'];

% FLIRT BBR
% The output image has single volume only. Ignore 'applywarp' to the input
% 4D image for now, since it will generate an image of very large size.
flirtbbrcmd = [FLIRT, ...
               ' -ref ', ea_path_helper(fixedimage), ...
               ' -in ', ea_path_helper(movingimage), ...
               ' -wmseg ', ea_path_helper(fixedimage_fast), '_wmseg', ...
               ' -init ', ea_path_helper(volumedir), xfm, '_init.mat', ...
               ' -omat ', ea_path_helper(volumedir), xfm, '.mat', ...
               ' -out ', ea_path_helper(outputimage), ...
               ' -dof 6', ...
               ' -cost bbr', ...
               ' -interp spline', ...
               ' -v'];

% Output inverse xfm for possible further use
% FSL won't handle the inversion internally when applying the transformation
invxfm = [fixname, '2', movname, '_flirtbbr'];
if writeoutmat
    convertxfmcmd = [COVERT_XFM, ...
                  ' -omat ', ea_path_helper(volumedir), invxfm, '.mat', ...
                  ' -inverse ', ea_path_helper(volumedir), xfm, '.mat'];
end

setenv('FSLOUTPUTTYPE','NIFTI');
if ~ispc
    fprintf('\n\nRunning FSL FAST Segmentation: %s\n\n', fixedimage);
    system(['bash -c "', fastcmd, '"']);
    system(['bash -c "', wmsegcmd, '"']);
    system(['bash -c "', wmedgecmd, '"']);

    fprintf('\n\nRunning FSL FLIRT Pre-alignment: %s\n\n', movingimage);
    system(['bash -c "', flirtinitcmd, '"']);

    fprintf('\n\nRunning FSL FLIRT BBR: %s\n\n', movingimage);
    system(['bash -c "', flirtbbrcmd, '"']);
    if writeoutmat
        system(['bash -c "', convertxfmcmd, '"']);
    end
else
    fprintf('\n\nRunning FSL FAST Segmentation: %s\n\n', fixedimage);
    system(fastcmd);
    system(wmsegcmd);
    system(wmedgecmd);

    fprintf('\n\nRunning FSL FLIRT Pre-alignment: %s\n\n', movingimage);
    system(flirtinitcmd);

    fprintf('\n\nRunning FSL FLIRT BBR: %s\n\n', movingimage);
    system(flirtbbrcmd);
    if writeoutmat
        system(convertxfmcmd);
    end
end

% Apply the tranformation to other files
% Be aware of the naming of the output, it will override the original image
if ~isempty(otherfiles)
    for fi = 1:numel(otherfiles)
        ea_fsl_apply_coregistration(fixedimage, otherfiles{fi}, otherfiles{fi}, ...
                                    [volumedir, xfm, '.mat']);
    end
end

if ~writeoutmat
    ea_delete([volumedir, xfm, '.mat']);
    ea_delete([volumedir, invxfm, '.mat']);
    affinefile = {};
else
    affinefile = {[volumedir, xfm, '.mat']
                  [volumedir, invxfm, '.mat']};
end

% Clean up initial transform for FSL FLIRT BBR
ea_delete([volumedir, xfm, '_init.mat']);

% Clean up FSL FAST results
ea_delete([fixedimage_fast, '*']);

fprintf('\nFSL FLIRT BBR done.\n');

%% add methods dump:
cits={
    'M. Jenkinson and S.M. Smith. A global optimisation method for robust affine registration of brain images. Medical Image Analysis, 5(2):143-156, 2001.'
    'M. Jenkinson, P.R. Bannister, J.M. Brady, and S.M. Smith. Improved optimisation for the robust and accurate linear registration and motion correction of brain images. NeuroImage, 17(2):825-841, 2002.'
    'D.N. Greve and B. Fischl, Accurate and robust brain image alignment using boundary-based registration, NeuroImage, 48(1):63-72, 2009.'
};

ea_methods(volumedir,[movname,' was linearly co-registered to ',fixname,' using FLIRT BBR as implemented in FSL (Jenkinson 2001; Jenkinson 2002; https://fsl.fmrib.ox.ac.uk/)'],...
    cits);
