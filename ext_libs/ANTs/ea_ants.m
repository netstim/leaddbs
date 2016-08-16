function ea_ants(varargin)
% Wrapper for ANTs linear registration

fixedimage=varargin{1};
movingimage=varargin{2};
outputimage=varargin{3};

if nargin>3
    writematout=varargin{4};
else
    writematout=1;
end
if nargin>4
    otherfiles=varargin{5};
end

if fileparts(movingimage)
    volumedir = [fileparts(movingimage), filesep];
else
    volumedir =['.', filesep];
end

[outputdir, outputname, ~] = fileparts(outputimage);
if outputdir
    outputbase = [outputdir, filesep, outputname];
else
    outputbase = ['.', filesep, outputname];
end

fixedimage = ea_path_helper(ea_niigz(fixedimage));
movingimage = ea_path_helper(ea_niigz(movingimage));
outputimage = ea_path_helper(ea_niigz(outputimage));

basedir = [fileparts(mfilename('fullpath')), filesep];

if ispc
    HEADER = [basedir, 'PrintHeader.exe'];
    ANTS = [basedir, 'antsRegistration.exe'];
else
    HEADER = [basedir, 'PrintHeader.', computer('arch')];
    ANTS = [basedir, 'antsRegistration.', computer('arch')];
end

if ~ispc
    [~, imgsize] = system(['bash -c "', HEADER, ' ',fixedimage, ' 2"']);
else
    [~, imgsize] = system([HEADER, ' ', fixedimage, ' 2']);
end
imgsize = cellfun(@(x) str2double(x),ea_strsplit(imgsize,'x'));

if any(imgsize>256)
    rigidconvergence='[1000x500x250x0,1e-6,10]';
    rigidshrinkfactors='12x8x4x2';
    rigidsoomthingssigmas='4x3x2x1vox';

    affineconvergence='[1000x500x250x0,1e-6,10]';
    affineshrinkfactors='12x8x4x2';
    affinesoomthingssigmas='4x3x2x1vox';
else
    rigidconvergence='[1000x500x250x0,1e-6,10]';
    rigidshrinkfactors='8x4x2x1';
    rigidsoomthingssigmas='3x2x1x0vox';

    affineconvergence='[1000x500x250x0,1e-6,10]';
    affineshrinkfactors='8x4x2x1';
    affinesoomthingssigmas='3x2x1x0vox';
end

% determine how many runs have been performed before..
runs=0;
for r=1:5;
    if exist([volumedir, 'ct2anat',num2str(r),'.mat'],'file')
        runs=r;
    else
        break
    end
end

if runs==0 % mattes MI affine + rigid
    rigidstage = [' --transform Rigid[0.1]' ...
    ' --convergence ', rigidconvergence, ...
    ' --shrink-factors ', rigidshrinkfactors, ...
    ' --smoothing-sigmas ', rigidsoomthingssigmas, ...
    ' --initial-moving-transform [', fixedimage, ',', movingimage, ',1]', ...
    ' --metric Mattes[', fixedimage, ',', movingimage, ',1,32,Regular,0.25]'];

    affinestage = [' --transform Affine[0.1]'...
    ' --metric MI[', fixedimage, ',', movingimage, ',1,32,Regular,0.25]' ...
    ' --convergence ', affineconvergence, ...
    ' --shrink-factors ', affineshrinkfactors ...
    ' --smoothing-sigmas ', affinesoomthingssigmas];

elseif runs==1
    rigidstage = [' --transform Rigid[0.1]' ...
        ' --convergence ', rigidconvergence, ...
        ' --shrink-factors ', rigidshrinkfactors, ...
        ' --smoothing-sigmas ', rigidsoomthingssigmas, ...
        ' --initial-moving-transform ',ea_path_helper([volumedir, 'ct2anat',num2str(runs),'.mat']), ...
        ' --metric GC[', fixedimage, ',', movingimage, ',1,32,Regular,0.25]'];

    affinestage = [' --transform Affine[0.1]'...
        ' --metric MI[', fixedimage, ',', movingimage, ',1,32,Regular,0.25]' ...
        ' --convergence ', affineconvergence, ...
        ' --shrink-factors ', affineshrinkfactors ...
        ' --smoothing-sigmas ', affinesoomthingssigmas];

elseif runs==2 % go directly to affine stage, try mattes MI
    rigidstage = '';
        affinestage = [
        ' --initial-moving-transform ',ea_path_helper([volumedir, 'ct2anat',num2str(runs),'.mat']), ...
        ' --transform Affine[0.1]'...
        ' --metric MI[', fixedimage, ',', movingimage, ',1,32,Regular,0.25]' ...
        ' --convergence ', affineconvergence, ...
        ' --shrink-factors ', affineshrinkfactors ...
        ' --smoothing-sigmas ', affinesoomthingssigmas];
elseif runs>=3 % go directly to affine stage, try GC again
    rigidstage = '';
        affinestage = [
        ' --initial-moving-transform ',ea_path_helper([volumedir, 'ct2anat',num2str(runs),'.mat']), ...
        ' --transform Affine[0.1]'...
        ' --metric GC[', fixedimage, ',', movingimage, ',1,32,Regular,0.25]' ...
        ' --convergence ', affineconvergence, ...
        ' --shrink-factors ', affineshrinkfactors ...
        ' --smoothing-sigmas ', affinesoomthingssigmas];
end

ea_libs_helper
cmd = [ANTS, ' --verbose 1' ...
    ' --dimensionality 3 --float 1' ...
    ' --output [',ea_path_helper(outputbase), ',', outputimage, ']' ...
    ' --interpolation Linear' ...
    ' --use-histogram-matching 1' ...
    ' --winsorize-image-intensities [0.005,0.995]', ...
    rigidstage, affinestage];
if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end

if exist('otherfiles','var')
    if ~isempty(otherfiles)
        for ofi=1:length(otherfiles)
        [options.root,options.patientname]=fileparts(fileparts(otherfiles{ofi}));
        options.root=[options.root,filesep];
        options.prefs=ea_prefs(options.patientname);
        ea_ants_applytransforms(options,otherfiles(ofi),otherfiles(ofi),0,fixedimage,[outputbase, '0GenericAffine.mat']);
        end
    end
end

if ~writematout
    delete([outputbase, '0GenericAffine.mat']);
else
    movefile([outputbase, '0GenericAffine.mat'], [volumedir, 'ct2anat',num2str(runs+1),'.mat']);
end
