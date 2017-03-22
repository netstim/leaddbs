function affinefile = ea_flirt(varargin)
% Wrapper for FSL linear registration

fixedimage = varargin{1};
movingimage = varargin{2};
outputimage = varargin{3};

if nargin >= 4
    writeoutmat = varargin{4};
else
    writeoutmat = 0;
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

fprintf(['\n\nRunning FSL FLIRT: ', movingimage, '\n\n']);

% Prepare bet image for flirt, generate the brain masks '*_bet_mask.nii'
[movpath, movname] = ea_niifileparts(movingimage);
[fixpath, fixname] = ea_niifileparts(fixedimage);
movingimage_bet = [fileparts(movpath), filesep, 'bet_', movname];
fixedimage_bet = [fileparts(fixpath), filesep, 'bet_', fixname];
if isempty(dir([movingimage_bet,'.nii*']))
    ea_bet(movingimage, 1, movingimage_bet);
end
if isempty(dir([fixedimage_bet,'.nii*']))
    ea_bet(fixedimage, 1, fixedimage_bet);
end

volumedir = [fileparts(ea_niifileparts(movingimage)), filesep];

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
if ispc
    FLIRT = [basedir, 'flirt.exe'];
    COVERT_XFM = [basedir, 'convert_xfm.exe'];
else
    FLIRT = [basedir, 'flirt.', computer('arch')];
    COVERT_XFM = [basedir, 'convert_xfm.', computer('arch')];
end

flirtcmd = [FLIRT, ...
            ' -ref ', ea_path_helper(fixedimage_bet), ...
            ' -in ', ea_path_helper(movingimage_bet), ...
            ' -out ', ea_path_helper(outputimage) ...
            affinestage];

% Output inverse xfm for possible further use
% FSL won't handle the inversion internally when applying the transformation
invxfm = [fix, '2', mov, '_flirt'];
convertxfmcmd = [COVERT_XFM, ...
              ' -omat ', ea_path_helper([volumedir, invxfm, num2str(runs+1), '.mat']), ...
              ' -inverse ', ea_path_helper([volumedir, xfm, num2str(runs+1), '.mat'])];

setenv('FSLOUTPUTTYPE','NIFTI');
if ~ispc
    system(['bash -c "', flirtcmd, '"']);
    system(['bash -c "', convertxfmcmd, '"']);
else
    system(flirtcmd);
    system(convertxfmcmd);
end

if ~isempty(otherfiles)
    for fi = 1:numel(otherfiles)
        ea_fsl_flirt_applytransform(fixedimage, otherfiles{fi}, otherfiles{fi}, ...
                                    [volumedir, xfm, num2str(runs+1), '.mat']);
    end
end

if ~writeoutmat
    ea_delete([volumedir, xfm, num2str(runs+1), '.mat'])
    ea_delete([volumedir, invxfm, num2str(runs+1), '.mat'])
    affinefile = {};
else
    affinefile = {[volumedir, xfm, num2str(runs+1), '.mat'], ...
                  [volumedir, invxfm, num2str(runs+1), '.mat']};
end

fprintf('\nFSL FLIRT done.\n');

%% add methods dump:
cits={
    'M. Jenkinson and S.M. Smith. A global optimisation method for robust affine registration of brain images. Medical Image Analysis, 5(2):143-156, 2001.'
    'M. Jenkinson, P.R. Bannister, J.M. Brady, and S.M. Smith. Improved optimisation for the robust and accurate linear registration and motion correction of brain images. NeuroImage, 17(2):825-841, 2002.'
    };

ea_methods(volumedir,[mov,' was linearly co-registered to ',fix,' using FLIRT as implemented in FSL (Jenkinson 2001; Jenkinson 2002; https://fsl.fmrib.ox.ac.uk/)'],...
    cits);

