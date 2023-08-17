function cmd =ea_ants_schoenecker(varargin)
% Wrapper for ANTs nonlinear registration

fixedimage=varargin{1};
movingimage=varargin{2};
outputimage=varargin{3};

if ischar(fixedimage)
    fixedimage={fixedimage};
elseif ~iscell(fixedimage)
    ea_error('Please supply variable fixedimage as either char or cellstring');
end

if ischar(movingimage)
    movingimage={movingimage};
elseif ~iscell(movingimage)
    ea_error('Please supply variable fixedimage as either char or cellstring');
end

if nargin >= 4
    weights = varargin{4};
else
    weights = ones(length(fixedimage),1);
end

if nargin >= 5
    metrics = varargin{5};
else
    metrics = repmat({'MI'},length(fixedimage),1);
end

slabsupport=1; % check for slabs in anat files and treat slabs differently (add additional SyN stage only in which slabs are being used).

[outputdir, outputname, ~] = fileparts(outputimage);
ea_mkdir(outputdir);
if outputdir
    outputbase = [outputdir, filesep, outputname];
else
    outputbase = ['.', filesep, outputname];
end

if slabsupport
    if isBIDSFileName(movingimage{end})
        parsedStruct = parseBIDSFilePath(movingimage{end});
        anchorName = parsedStruct.suffix;
    else
        [~, anchorName] = ea_niifileparts(movingimage{end});
    end
    disp(['Checking for slabs among structural images (assuming anchor image ',anchorName,' is a whole-brain acquisition)...']);

    for mov=1:length(movingimage)
        mnii=ea_load_nii(movingimage{mov});
        mnii.img=~(mnii.img==0);
        if ~exist('AllMX','var')
            AllMX=mnii.img;
        else
            AllMX=AllMX.*mnii.img;
        end
        sums(mov)=sum(mnii.img(:));
    end
    slabspresent=0; % default no slabs present.

    if length(sums)>1 % multispectral warp
        slabs=sums(1:end-1)<(sums(end)*0.7);
        if any(slabs) % one image is smaller than 70% of last (dominant) image, a slab is prevalent.
            slabmovingimage=ea_path_helper(movingimage(slabs)); % move slabs to new cell slabimage
            slabfixedimage=ea_path_helper(fixedimage(slabs));
            movingimage(slabs)=[]; % remove slabs from movingimage
            fixedimage(slabs)=[]; % remove slabs from fixedimage

            % write out slab mask
            slabspresent=1;
            mnii.dt(1) = 4;
            mnii.img=AllMX;

            tmaskdir = fullfile(outputdir, 'tmp');
            if ~exist(tmaskdir, 'dir')
                mkdir(tmaskdir);
            end

            mnii.fname=[tmaskdir,filesep,'slabmask.nii'];
            ea_write_nii(mnii);
            disp('Slabs found. Separating slabs to form an additional SyN stage.');
            % slabmovingimage = mnii.fname;
        else
            disp('No slabs found.');
        end
    end

else
    slabspresent=0;
    impmasks=repmat({'nan'},length(movingimage),1);
end

for fi=1:length(fixedimage)
    fixedimage{fi} = ea_path_helper(ea_niigz(fixedimage{fi}));
end
for fi=1:length(movingimage)
    movingimage{fi} = ea_path_helper(ea_niigz(movingimage{fi}));
end

if length(fixedimage)~=length(movingimage)
    ea_error('Please supply pairs of moving and fixed images (can be repetitive).');
end

outputimage = ea_niigz(outputimage);

basedir = [fileparts(mfilename('fullpath')), filesep];
HEADER = ea_getExec([basedir, 'PrintHeader'], escapePath = 1);
ANTS = ea_getExec([basedir, 'antsRegistration'], escapePath = 1);
applyTransforms = ea_getExec([basedir, 'antsApplyTransforms'], escapePath = 1);

headercmd = [HEADER, ' ', fixedimage{1}, ' 2'];
[~, imgsize] = ea_runcmd(headercmd);

imgsize = cellfun(@(x) str2double(x),ea_strsplit(imgsize,'x'));

% Rigid stage
rigidconvergence='[100x50x25x10,1e-6,10]';
rigidshrinkfactors='8x4x2x1';
rigidsmoothingssigmas='3x2x1x0';

% Affine stage
affineconvergence='[1000x500x250x100,1e-6,10]';
affineshrinkfactors='8x4x2x1';
affinesmoothingssigmas='3x2x1x0';

% 1. Mask stage
affinemask1convergence='[500x250x100,1e-6,10]';
affinemask1shrinkfactors='4x2x1';
affinemask1smoothingssigmas='2x1x0';

% 2. Mask stage
affinemask2convergence='[100,1e-6,10]';
affinemask2shrinkfactors='1';
affinemask2smoothingssigmas='0';


% Rigid stage
rigidstage = [' --initial-moving-transform [', fixedimage{1}, ',', movingimage{1}, ',1]' ...
    ' --transform Rigid[0.1]' ...
    ' --convergence ', rigidconvergence, ...
    ' --shrink-factors ', rigidshrinkfactors, ...
    ' --smoothing-sigmas ', rigidsmoothingssigmas, ...
    ' --masks [NULL,NULL]'];

for fi=1:length(fixedimage)
    switch metrics{fi}
        case 'MI'
            suffx=',32,Regular,0.25';
        case 'CC'
            suffx=',4';
        case 'GC'
            suffx=',15,Random,0.05';
    end

    try
        rigidstage=[rigidstage,...
            ' --metric ','MI','[', fixedimage{fi}, ',', movingimage{fi}, ',',num2str(weights(fi)),suffx,']'];
    catch
        keyboard
    end
end

% Affine stage
affinestage = [' --transform Affine[0.1]'...
    ' --convergence ', affineconvergence, ...
    ' --shrink-factors ', affineshrinkfactors ...
    ' --smoothing-sigmas ', affinesmoothingssigmas, ...
    ' --masks [NULL,NULL]'];

for fi=1:length(fixedimage)
    switch metrics{fi}
        case 'MI'
            suffx=',32,Regular,0.25';
        case 'CC'
            suffx=',4';
        case 'GC'
            suffx=',15,Random,0.05';
    end
    affinestage=[affinestage,...
        ' --metric ','MI','[', fixedimage{fi}, ',', movingimage{fi}, ',',num2str(weights(fi)),suffx,']'];
end

% 1. Mask stage
if slabsupport && slabspresent
    % re-add slabs to the masked stages:

    fixedimage=[fixedimage,slabfixedimage];
    movingimage=[movingimage,slabmovingimage];
end

affinestage_mask1 = [' --transform Affine[0.1]'...
    ' --convergence ', affinemask1convergence, ...
    ' --shrink-factors ', affinemask1shrinkfactors ...
    ' --smoothing-sigmas ', affinemask1smoothingssigmas, ...
    ' --masks [',ea_space([],'subcortical'),'secondstepmask.nii',',NULL]'];

for fi=1:length(fixedimage)
    switch metrics{fi}
        case 'MI'
            suffx=',32,Regular,0.5';
        case 'CC'
            suffx=',4';
        case 'GC'
            suffx=',15,Random,0.05';
    end
    affinestage_mask1=[affinestage_mask1,...
        ' --metric ',metrics{fi},'[', fixedimage{fi}, ',', movingimage{fi}, ',',num2str(weights(fi)),suffx,']'];
end

% 2. Mask stage
affinestage_mask2 = [' --transform Affine[0.1]'...
    ' --convergence ', affinemask2convergence, ...
    ' --shrink-factors ', affinemask2shrinkfactors ...
    ' --smoothing-sigmas ', affinemask2smoothingssigmas, ...
    ' --masks [',ea_space([],'subcortical'),'thirdstepmask.nii',',NULL]'];

for fi=1:length(fixedimage)
    switch metrics{fi}
        case 'MI'
            suffx=',32,Regular,1';
        case 'CC'
            suffx=',4';
        case 'GC'
            suffx=',15,Random,0.05';
    end
    affinestage_mask2=[affinestage_mask2,...
        ' --metric ',metrics{fi},'[', fixedimage{fi}, ',', movingimage{fi}, ',',num2str(weights(fi)),suffx,']'];
end

ea_libs_helper;

cmd = [ANTS, ' --verbose 1', ...
    ' --dimensionality 3', ...
    ' --output [',ea_path_helper(outputbase), ',', ea_path_helper(outputimage), ']', ...
    ' --interpolation Linear', ...
    ' --use-histogram-matching 1', ...
    ' --float 1',...
    ' --write-composite-transform 0', ...
    rigidstage, affinestage, affinestage_mask1, affinestage_mask2];

invcmd = [applyTransforms, ' --verbose 1' ...
    ' --dimensionality 3 --float 1' ...
    ' --reference-image ', movingimage{end}, ...
    ' --transform [', ea_path_helper([outputbase, '0GenericAffine.mat']),',1]' ...
    ' --output Linear[', ea_path_helper([outputbase, 'Inverse0GenericAffine.mat']),']'];

% Log ANTs command
if isBIDSFileName(outputimage)
    parsedStruct = parseBIDSFilePath(outputimage);
    logDir = strrep(parsedStruct.dir, ['normalization', filesep, 'anat'], ['normalization', filesep, 'log']);
    ea_mkdir(logDir);
    antsCMDFile = [logDir, filesep, 'sub-', parsedStruct.sub, '_desc-antscmd.txt'];
else
    antsCMDFile = [outputdir, filesep, 'ea_ants_command.txt'];
end

fid = fopen(antsCMDFile, 'a');
fprintf(fid, '%s:\n%s\n\n', char(datetime('now')), cmd);
fclose(fid);

ea_runcmd(cmd);
ea_runcmd(invcmd);

if exist('tmaskdir', 'var')
    ea_delete(tmaskdir);
end
