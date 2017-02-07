function affinefile = ea_flirt_bbr(varargin)
% Wrapper for FSL BBR linear registration

fixedimage = varargin{1};
movingimage = varargin{2};
outputimage = varargin{3};

if nargin >= 4
    if isempty(varargin{4}) % [] or {} or ''
        otherfiles = {};
    elseif ischar(varargin{4}) % single file, make it to cell string
        otherfiles = varargin(4);
    else % cell string
        otherfiles = varargin{4};
    end
else
    otherfiles = {};
end

% Prepare bet image for flirt bbr
[fixpath, fixname] = ea_niifileparts(fixedimage);
fixedimage_bet = [fileparts(fixpath), filesep, 'bet_', fixname];
if isempty(dir([fixedimage_bet,'.nii*']))
    ea_bet(fixedimage, 0, fixedimage_bet, 0.5);
end

volumedir = [fileparts(ea_niifileparts(movingimage)), filesep];

% name of the output transformation
[~, movname] = ea_niifileparts(movingimage);
xfm = [movname, '2', fixname, '_flirtbbr'];

basedir = [fileparts(mfilename('fullpath')), filesep];
if ispc
    FAST = [basedir, 'fast.exe'];
    FSLMATHS = [basedir, 'fslmaths.exe'];
    FLIRT = [basedir, 'flirt.exe'];
    COVERT_XFM = [basedir, 'convert_xfm.exe'];
else
    FAST = [basedir, 'fast.', computer('arch')];
    FSLMATHS = [basedir, 'fslmaths.', computer('arch')];
    FLIRT = [basedir, 'flirt.', computer('arch')];
    COVERT_XFM = [basedir, 'convert_xfm.', computer('arch')];
end

fixedimage_fast = [fileparts(fixpath), filesep, 'fast_', fixname];

% Segment the reference image
fastcmd = [FAST, ' -v -o ', ea_path_helper(fixedimage_fast), ' ', ea_path_helper(fixedimage_bet)];

% Creat white matter mask
wmsegcmd = [FSLMATHS, ...
            ' ', ea_path_helper(fixedimage_fast), '_pve_2', ...
            ' -thr 0.5 -bin', ...
            ' ', ea_path_helper(fixedimage_fast), '_wmseg'];


% Creat white matter edge map （just for visualization）
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
convertxfmcmd = [COVERT_XFM, ...
              ' -omat ', ea_path_helper(volumedir), invxfm, '.mat', ...
              ' -inverse ', ea_path_helper(volumedir), xfm, '.mat'];

setenv('FSLOUTPUTTYPE','NIFTI');
if ~ispc
    fprintf(['\n\nRunning FSL FAST Segmentation: ', fixedimage, '\n\n']);
    system(['bash -c "', fastcmd, '"']);
    system(['bash -c "', wmsegcmd, '"']);
    system(['bash -c "', wmedgecmd, '"']);

    fprintf(['\n\nRunning FSL FLIRT Pre-alignment: ', movingimage, '\n\n']);
    system(['bash -c "', flirtinitcmd, '"']);

    fprintf(['\n\nRunning FSL FLIRT BBR: ', movingimage, '\n\n']);
    system(['bash -c "', flirtbbrcmd, '"']);
    system(['bash -c "', convertxfmcmd, '"']);
else
    fprintf(['\n\nRunning FSL FAST Segmentation: ', fixedimage, '\n\n']);
    system(fastcmd);
    system(wmsegcmd);
    system(wmedgecmd);

    fprintf(['\n\nRunning FSL FLIRT Pre-alignment: ', movingimage, '\n\n']);
    system(flirtinitcmd);

    fprintf(['\n\nRunning FSL FLIRT BBR: ', movingimage, '\n\n']);
    system(flirtbbrcmd);
    system(convertxfmcmd);
end

% Apply the tranformation to other files
% Be aware of the naming of the output, it will override the original image
if ~isempty(otherfiles)
    for fi = 1:numel(otherfiles)
        ea_fsl_flirt_applytransform(fixedimage, otherfiles{fi}, otherfiles{fi}, ...
                                    [volumedir, xfm, '.mat']);
    end
end

affinefile = {[volumedir, xfm, '.mat'], ...
              [volumedir, invxfm, '.mat']};

fprintf('\nFSL FLIRT BBR done.\n');
