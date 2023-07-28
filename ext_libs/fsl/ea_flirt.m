function affinefile = ea_flirt(varargin)
% Wrapper for FSL linear registration

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

fprintf('\n\nRunning FSL FLIRT: %s\n\n', movingimage);

umachine = load([ea_gethome, '.ea_prefs.mat']);
normsettings = umachine.machine.normsettings;
% Prepare bet image for flirt
if normsettings.fsl_skullstrip % skullstripping is on
    % Set skullstripped file name
    if isBIDSFileName(movingimage)
        parsedStruct = parseBIDSFilePath(movingimage);
        if isfield(parsedStruct, 'acq')
            inimage = ea_niifileparts(setBIDSEntity(movingimage, 'acq', [], 'label', 'Brain', 'acq', parsedStruct.acq));
        else
            inimage = ea_niifileparts(setBIDSEntity(movingimage, 'label', 'Brain'));
        end
    else
        inimage = [ea_niifileparts(movingimage), '_brain'];
    end

    % Run BET2
    fprintf('\nSkullstripping moving image...\n');
    ea_bet(movingimage, 1, inimage);

    % Rename mask file
    mask = dir([inimage, '_mask*']);
    ext = regexp(mask(end).name, '(?<=_mask)\.nii(\.gz)?$', 'match', 'once');
    if isBIDSFileName(movingimage)
        parsedStruct = parseBIDSFilePath(movingimage);
        movefile([inimage, '_mask', ext], setBIDSEntity(movingimage, 'mod', parsedStruct.suffix, 'label', 'Brain', 'suffix', 'mask'));
    else
        movefile([inimage, '_mask', ext], [ea_niifileparts(movingimage), '_brainmask', ext]);
    end

	% Set skullstripped file name
    if isBIDSFileName(fixedimage)
        parsedStruct = parseBIDSFilePath(fixedimage);
        if isfield(parsedStruct, 'acq')
            refimage = ea_niifileparts(setBIDSEntity(fixedimage, 'acq', [], 'label', 'Brain', 'acq', parsedStruct.acq));
        else
            refimage = ea_niifileparts(setBIDSEntity(fixedimage, 'label', 'Brain'));
        end
    else
        refimage = [ea_niifileparts(fixedimage), '_brain'];
    end

    % Run BET2
    fprintf('\nSkullstripping reference image...\n');
    ea_bet(fixedimage, 1, refimage);

    % Rename mask file
    mask = dir([refimage, '_mask*']);
    ext = regexp(mask(end).name, '(?<=_mask)\.nii(\.gz)?$', 'match', 'once');
    if isBIDSFileName(fixedimage)
        parsedStruct = parseBIDSFilePath(fixedimage);
        movefile([refimage, '_mask', ext], setBIDSEntity(fixedimage, 'mod', parsedStruct.suffix, 'label', 'Brain', 'suffix', 'mask'));
    else
        movefile([refimage, '_mask', ext], [ea_niifileparts(fixedimage), '_brainmask', ext]);
    end
else % skullstripping is off
    fprintf('\nSkip skullstripping...\n\n');
    inimage = ea_niifileparts(movingimage);
    refimage = ea_niifileparts(fixedimage);
end

volumedir = [fileparts(ea_niifileparts(outputimage)), filesep];

% name of the output transformation
[~, mov] = ea_niifileparts(movingimage);
[~, fix] = ea_niifileparts(fixedimage);
xfm = [mov, '2', fix, '_flirt'];
% determine how many runs have been performed before
runs = dir([volumedir, xfm, '*.mat']);
if isempty(runs)
    runs = 0;
else
    runs = str2double(runs(end).name(numel(xfm)+1:end-4)); % suppose runs<10
end

% affine params
affinestage = [' -cost mutualinfo' ...
               ' -searchcost mutualinfo' ...
               ' -interp sinc' ...
               ' -verbose 1'];

if runs == 0 % mattes MI affine + rigid
    affinestage = [affinestage, ...
                   ' -omat ', ea_path_helper([volumedir, xfm ,num2str(runs+1), '.mat'])];
elseif runs > 0
    affinestage = [affinestage, ...
                   ' -init ', ea_path_helper([volumedir, xfm, num2str(runs), '.mat']), ...
                   ' -omat ', ea_path_helper([volumedir, xfm, num2str(runs+1), '.mat'])];
end

basedir = [fileparts(mfilename('fullpath')), filesep];
FLIRT = ea_getExec([basedir, 'flirt'], escapePath = 1);
COVERT_XFM = ea_getExec([basedir, 'convert_xfm'], escapePath = 1);

flirtcmd = [FLIRT, ...
            ' -ref ', ea_path_helper(refimage), ...
            ' -in ', ea_path_helper(inimage), ...
            ' -out ', ea_path_helper(outputimage) ...
            affinestage];

% Output inverse xfm for possible further use
% FSL won't handle the inversion internally when applying the transformation
invxfm = [fix, '2', mov, '_flirt'];
convertxfmcmd = [COVERT_XFM, ...
              ' -omat ', ea_path_helper([volumedir, invxfm, num2str(runs+1), '.mat']), ...
              ' -inverse ', ea_path_helper([volumedir, xfm, num2str(runs+1), '.mat'])];

setenv('FSLOUTPUTTYPE','NIFTI');
ea_runcmd(flirtcmd);
ea_runcmd(convertxfmcmd);

if ~isempty(otherfiles)
    for fi = 1:numel(otherfiles)
        ea_fsl_apply_coregistration(fixedimage, otherfiles{fi}, otherfiles{fi}, ...
                                    [volumedir, xfm, num2str(runs+1), '.mat']);
    end
end

if ~writeoutmat
    ea_delete([volumedir, xfm, num2str(runs+1), '.mat'])
    ea_delete([volumedir, invxfm, num2str(runs+1), '.mat'])
    affinefile = {};
else
    affinefile = {[volumedir, xfm, num2str(runs+1), '.mat']
                  [volumedir, invxfm, num2str(runs+1), '.mat']};
end

% Clean up BET image when skullstripping is on
if normsettings.fsl_skullstrip
    ea_delete({inimage, refimage});
end

fprintf('\nFSL FLIRT done.\n');

%% add methods dump:
cits={
    'M. Jenkinson and S.M. Smith. A global optimisation method for robust affine registration of brain images. Medical Image Analysis, 5(2):143-156, 2001.'
    'M. Jenkinson, P.R. Bannister, J.M. Brady, and S.M. Smith. Improved optimisation for the robust and accurate linear registration and motion correction of brain images. NeuroImage, 17(2):825-841, 2002.'
};

ea_methods(volumedir,[mov,' was linearly co-registered to ',fix,' using FLIRT as implemented in FSL (Jenkinson 2001; Jenkinson 2002; https://fsl.fmrib.ox.ac.uk/)'],cits);

