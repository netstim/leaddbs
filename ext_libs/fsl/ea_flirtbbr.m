function affinefile = ea_flirtbbr(varargin)
% Wrapper for FSL FLIRT BBR registration

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
fixedPath = ea_niifileparts(fixedimage);
if isBIDSFileName(fixedimage)
    parsedStruct = parseBIDSFilePath(fixedimage);
    if isfield(parsedStruct, 'acq')
        fixedimage_bet = ea_niifileparts(setBIDSEntity(fixedimage, 'acq', [], 'label', 'Brain', 'acq', parsedStruct.acq));
    else
        fixedimage_bet = ea_niifileparts(setBIDSEntity(fixedimage, 'label', 'Brain'));
    end
else
    fixedimage_bet = [fixedPath, '_brain'];
end

fprintf('\nSkullstripping fixed image...\n');
ea_bet(fixedimage, 1, fixedimage_bet, 0.5);

% Rename mask
mask = dir([fixedimage_bet, '_mask*']);
ext = regexp(mask(end).name, '(?<=_mask)\.nii(\.gz)?$', 'match', 'once');
if isBIDSFileName(fixedimage)
    parsedStruct = parseBIDSFilePath(fixedimage);
    movefile([fixedimage_bet, '_mask', ext], setBIDSEntity(fixedimage, 'label', 'Brain', 'mod', parsedStruct.suffix, 'suffix', 'mask'));
else
    movefile([fixedimage_bet, '_mask', ext], [fixedPath, '_brainmask', ext]);
end

volumedir = [fileparts(ea_niifileparts(outputimage)), filesep];

% name of the output transformation
[~, movname] = ea_niifileparts(movingimage);
xfm = [movname, '2', fixname, '_flirtbbr'];

basedir = [fileparts(mfilename('fullpath')), filesep];
FAST = ea_getExec([basedir, 'fast'], escapePath = 1);
FSLMATHS = ea_getExec([basedir, 'fslmaths'], escapePath = 1);
FLIRT = ea_getExec([basedir, 'flirt'], escapePath = 1);
COVERT_XFM = ea_getExec([basedir, 'convert_xfm'], escapePath = 1);

fixedimage_fast = [fixedPath, '_fslfast'];

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

fprintf('\n\nRunning FSL FAST Segmentation: %s\n\n', fixedimage);
ea_runcmd(fastcmd);
ea_runcmd(wmsegcmd);
ea_runcmd(wmedgecmd);

fprintf('\n\nRunning FSL FLIRT Pre-alignment: %s\n\n', movingimage);
ea_runcmd(flirtinitcmd);

fprintf('\n\nRunning FSL FLIRT BBR: %s\n\n', movingimage);
ea_runcmd(flirtbbrcmd);

if writeoutmat
    ea_runcmd(convertxfmcmd);
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

% Clean up BET image
ea_delete(fixedimage_bet);

fprintf('\nFSL FLIRT BBR done.\n');

%% add methods dump:
cits={
    'M. Jenkinson and S.M. Smith. A global optimisation method for robust affine registration of brain images. Medical Image Analysis, 5(2):143-156, 2001.'
    'M. Jenkinson, P.R. Bannister, J.M. Brady, and S.M. Smith. Improved optimisation for the robust and accurate linear registration and motion correction of brain images. NeuroImage, 17(2):825-841, 2002.'
    'D.N. Greve and B. Fischl, Accurate and robust brain image alignment using boundary-based registration, NeuroImage, 48(1):63-72, 2009.'
};

ea_methods(volumedir,[movname,' was linearly co-registered to ',fixname,' using FLIRT BBR as implemented in FSL (Jenkinson 2001; Jenkinson 2002; https://fsl.fmrib.ox.ac.uk/)'],cits);
