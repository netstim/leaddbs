function affinefile = ea_ants(varargin)
% Wrapper for ANTs linear registration


fixedimage=varargin{1};
movingimage=varargin{2};
outputimage=varargin{3};

if nargin >= 4
    writematout=varargin{4};
else
    writematout=1;
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

try
    msks=varargin{6};
    if isempty(msks)
        usemasks=0;
    else
        if ~iscell(msks)
            msks={};
            usemasks=0;
        else
            usemasks=1;
        end
    end
catch
    usemasks=0;
end

try
    options=varargin{7};
catch
    options=struct;
end

try
    interp=varargin{8};
catch
    interp='Linear';
end

outputbase = ea_niifileparts(outputimage);
volumedir = [fileparts(outputbase), filesep];

basedir = [fileparts(mfilename('fullpath')), filesep];
if ispc
    HEADER = [basedir, 'PrintHeader.exe'];
    ANTS = [basedir, 'antsRegistration.exe'];
    antsApplyTransforms = [basedir, 'antsApplyTransforms.exe'];
else
    HEADER = [basedir, 'PrintHeader.', computer('arch')];
    ANTS = [basedir, 'antsRegistration.', computer('arch')];
    antsApplyTransforms = [basedir, 'antsApplyTransforms.', computer('arch')];
end

if ~ispc
    [~, imgsize] = system(['bash -c "', HEADER, ' ', ea_path_helper(fixedimage), ' 2"']);
else
    [~, imgsize] = system([HEADER, ' ', ea_path_helper(fixedimage), ' 2']);
end
imgsize = cellfun(@(x) str2double(x),ea_strsplit(imgsize,'x'));

if any(imgsize>256)
    rigidconvergence='[1000x500x250x0,1e-6,10]';
    rigidshrinkfactors='12x8x4x2';
    rigidsoomthingssigmas='4x3x2x1vox';

    affineconvergence='[1000x500x250x0,1e-6,10]';
    affineshrinkfactors='8x4x2x2';
    affinesoomthingssigmas='4x3x2x1vox';
else
    rigidconvergence='[1000x500x250x0,1e-6,10]';
    rigidshrinkfactors='8x4x2x1';
    rigidsoomthingssigmas='3x2x1x0vox';

    affineconvergence='[1000x500x250x0,1e-6,10]';
    affineshrinkfactors='8x4x2x1';
    affinesoomthingssigmas='3x2x1x0vox';
end

% name of the output transformation
[~, mov] = ea_niifileparts(movingimage);
[~, fix] = ea_niifileparts(fixedimage);
xfm = [mov, '2', fix, '_ants'];
% determine how many runs have been performed before
runs = dir([volumedir, xfm, '*.mat']);
if isempty(runs)
    runs = 0;
else
    runs = str2double(runs(end).name(length(xfm)+1:end-4)); % suppose runs<10
end


if runs==0 % mattes MI affine + rigid
    rigidstage = [' --transform Rigid[0.1]' ...
    ' --convergence ', rigidconvergence, ...
    ' --shrink-factors ', rigidshrinkfactors, ...
    ' --smoothing-sigmas ', rigidsoomthingssigmas, ...
    ' --initial-moving-transform [', ea_path_helper(fixedimage), ',', ea_path_helper(movingimage), ',1]', ...
    ' --metric Mattes[', ea_path_helper(fixedimage), ',', ea_path_helper(movingimage), ',1,32,Regular,0.25]'];

    affinestage = [' --transform Affine[0.1]'...
    ' --metric MI[', ea_path_helper(fixedimage), ',', ea_path_helper(movingimage), ',1,32,Regular,0.25]' ...
    ' --convergence ', affineconvergence, ...
    ' --shrink-factors ', affineshrinkfactors ...
    ' --smoothing-sigmas ', affinesoomthingssigmas];


elseif runs==1
    rigidstage = [' --transform Rigid[0.1]' ...
        ' --convergence ', rigidconvergence, ...
        ' --shrink-factors ', rigidshrinkfactors, ...
        ' --smoothing-sigmas ', rigidsoomthingssigmas, ...
        ' --initial-moving-transform ',ea_path_helper([volumedir, xfm, num2str(runs), '.mat']), ...
        ' --metric GC[', ea_path_helper(fixedimage), ',', ea_path_helper(movingimage), ',1,32,Regular,0.25]'];

    affinestage = [' --transform Affine[0.1]'...
        ' --metric MI[', ea_path_helper(fixedimage), ',', ea_path_helper(movingimage), ',1,32,Regular,0.25]' ...
        ' --convergence ', affineconvergence, ...
        ' --shrink-factors ', affineshrinkfactors ...
        ' --smoothing-sigmas ', affinesoomthingssigmas];

elseif runs==2 % go directly to affine stage, try mattes MI
    rigidstage = [' --transform Rigid[0.1]' ...
    ' --convergence ', rigidconvergence, ...
    ' --shrink-factors ', rigidshrinkfactors, ...
    ' --smoothing-sigmas ', rigidsoomthingssigmas, ...
    ' --initial-moving-transform [', ea_path_helper(fixedimage), ',', ea_path_helper(movingimage), ',1]', ...
    ' --metric Mattes[', ea_path_helper(fixedimage), ',', ea_path_helper(movingimage), ',1,32,Regular,0.25]'];
    affinestage = [
        ' --initial-moving-transform ',ea_path_helper([volumedir, xfm, num2str(runs), '.mat']), ...
        ' --transform Affine[0.1]'...
        ' --metric MI[', ea_path_helper(fixedimage), ',', ea_path_helper(movingimage), ',1,32,Regular,0.25]' ...
        ' --convergence ', affineconvergence, ...
        ' --shrink-factors ', affineshrinkfactors ...
        ' --smoothing-sigmas ', affinesoomthingssigmas];

elseif runs>=3 % go directly to affine stage, try GC again
    rigidstage = [' --transform Rigid[0.1]' ...
    ' --convergence ', rigidconvergence, ...
    ' --shrink-factors ', rigidshrinkfactors, ...
    ' --smoothing-sigmas ', rigidsoomthingssigmas, ...
    ' --initial-moving-transform [', ea_path_helper(fixedimage), ',', ea_path_helper(movingimage), ',1]', ...
    ' --metric Mattes[', ea_path_helper(fixedimage), ',', ea_path_helper(movingimage), ',1,32,Regular,0.25]'];
    affinestage = [
        ' --initial-moving-transform ',ea_path_helper([volumedir, xfm, num2str(runs), '.mat']), ...
        ' --transform Affine[0.1]'...
        ' --metric GC[', ea_path_helper(fixedimage), ',', ea_path_helper(movingimage), ',1,32,Regular,0.25]' ...
        ' --convergence ', affineconvergence, ...
        ' --shrink-factors ', affineshrinkfactors ...
        ' --smoothing-sigmas ', affinesoomthingssigmas];
end




if usemasks
    usereaffine=1; % additional affine step based on mask is probably too much.

    rigidstage=[rigidstage, ... % add nonexisting mask for this stage
        ' --masks [nan,nan]'];
    affinestage=[affinestage, ... % add nonexisting mask for this stage
        ' --masks [nan,nan]'];
    mask1stage = [' --transform Affine[0.1]' ...
        ' --convergence ', rigidconvergence, ...
        ' --shrink-factors ', rigidshrinkfactors, ...
        ' --smoothing-sigmas ', rigidsoomthingssigmas, ...
        ' --initial-moving-transform [', ea_path_helper(fixedimage), ',', ea_path_helper(movingimage), ',1]', ...
        ' --metric Mattes[', ea_path_helper(fixedimage), ',', ea_path_helper(movingimage), ',1,32,Regular,0.25]', ...
        ' --masks [', ea_path_helper(msks{1}),',',ea_path_helper(msks{1}),']'];
    if length(msks)>1
        mask2stage = [' --transform Affine[0.1]'...
            ' --metric MI[', ea_path_helper(fixedimage), ',', ea_path_helper(movingimage), ',1,32,Regular,0.25]' ...
            ' --convergence ', affineconvergence, ...
            ' --shrink-factors ', affineshrinkfactors ...
            ' --smoothing-sigmas ', affinesoomthingssigmas ...
            ' --masks [', ea_path_helper(msks{2}),',',ea_path_helper(msks{2}),']'];
    else
        mask2stage='';
    end
else
    mask1stage='';
    mask2stage='';
end

ea_libs_helper;
antscmd = [ANTS, ' --verbose 1' ...
    ' --dimensionality 3 --float 1' ...
    ' --output [',ea_path_helper(outputbase), ',', ea_path_helper(outputimage), ']' ...
    ' --interpolation ',interp ...
    ' --use-histogram-matching 1' ...
    ' --winsorize-image-intensities [0.005,0.995]', ...
    rigidstage, affinestage,mask1stage,mask2stage];

if writematout % inverse only needed if matrix is written out.
    invaffinecmd = [antsApplyTransforms, ' --verbose 1' ...
        ' --dimensionality 3 --float 1' ...
        ' --reference-image ', ea_path_helper(movingimage), ...
        ' --transform [', ea_path_helper([outputbase, '0GenericAffine.mat']),',1]' ...
        ' --output Linear[', ea_path_helper([outputbase, 'Inverse0GenericAffine.mat']),']'];
end

if ~ispc
    system(['bash -c "', antscmd, '"']);
    if writematout
    system(['bash -c "', invaffinecmd, '"']);
    end
else
    system(antscmd);
    if writematout
    system(invaffinecmd);
    end
end

if ~isempty(otherfiles)
    for ofi=1:length(otherfiles)
        [options.root,options.patientname]=fileparts(fileparts(otherfiles{ofi}));
        options.root=[options.root,filesep];
        options.prefs=ea_prefs(options.patientname);
        ea_ants_applytransforms(options,otherfiles(ofi),otherfiles(ofi),0,fixedimage,[outputbase, '0GenericAffine.mat']);
    end
end

if ~writematout
    delete([outputbase, '0GenericAffine.mat']);
%    delete([outputbase, 'Inverse0GenericAffine.mat']); % will not be
%    generated anymore if not writematout
    affinefile = {};
else
    movefile([outputbase, '0GenericAffine.mat'], [volumedir, xfm, num2str(runs+1), '.mat']);
    invxfm = [fix, '2', mov, '_ants'];
    movefile([outputbase, 'Inverse0GenericAffine.mat'], [volumedir, invxfm, num2str(runs+1), '.mat']);
    affinefile = {[volumedir, xfm, num2str(runs+1), '.mat'], ...
                  [volumedir, invxfm, num2str(runs+1), '.mat']};
end

fprintf('\nANTs LINEAR registration done.\n');

%% add methods dump:
cits={
    'Avants, B. B., Epstein, C. L., Grossman, M., & Gee, J. C. (2008). Symmetric diffeomorphic image registration with cross-correlation: evaluating automated labeling of elderly and neurodegenerative brain. Medical Image Analysis, 12(1), 26?41. http://doi.org/10.1016/j.media.2007.06.004'
    };

ea_methods(volumedir,[mov,' was co-registered to ',fix,' using a two-stage linear registration (rigid followed by affine) as implemented in Advanced Normlization Tools (Avants 2008; http://stnava.github.io/ANTs/)'],...
    cits);

