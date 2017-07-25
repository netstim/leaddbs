function ea_ants_nonlinear_monostep(varargin)
% Wrapper for ANTs nonlinear registration

fixedimage=varargin{1};
movingimage=varargin{2};
outputimage=varargin{3};

if ischar(movingimage)
    movingimage={movingimage};
elseif ~iscell(movingimage)
    ea_error('Please supply variable fixedimage as either char or cellstring');
end
try
    subcorticalrefine=varargin{6};
catch
    subcorticalrefine=0;
end
slabsupport=1; % check for slabs in anat files and treat slabs differently (add additional SyN stage only in which slabs are being used).

[outputdir, outputname, ~] = fileparts(outputimage);
if outputdir
    outputbase = [outputdir, filesep, outputname];
else
    outputbase = ['.', filesep, outputname];
end

if ischar(fixedimage)
    fixedimage={fixedimage};
elseif ~iscell(fixedimage)
    ea_error('Please supply variable fixedimage as either char or cellstring');
end

if slabsupport
    disp(['Checking for slabs among structural images (assuming dominant structural file ',movingimage{end},' is a whole-brain acquisition)...']);
    tmaskdir=fullfile(tempdir,'lead');
    if ~exist(tmaskdir,'dir')
        mkdir(tmaskdir);
    end
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
        if any(slabs) % one image is smaller than 0.7% of last (dominant) image, a slab is prevalent.
            slabmovingimage=movingimage(slabs); % move slabs to new cell slabimage
            slabfixedimage=fixedimage(slabs);
            movingimage(slabs)=[]; % remove slabs from movingimage
            fixedimage(slabs)=[]; % remove slabs from movingimage
            
            % write out slab mask
            slabspresent=1;
            slabID=ea_generate_guid;
            mnii.dt=[4,0];
            mnii.img=AllMX;
            mnii.fname=[tmaskdir,filesep,'slabmask_',slabID,'.nii'];
            ea_write_nii(mnii);
            disp('Slabs found. Separating slabs to form an additional SyN stage.');
        else
            disp('No slabs found.');
        end
    end
else
    slabspresent=0;
    impmasks=repmat({'nan'},length(movingimage),1);
end

if nargin>3
    weights=varargin{4};
    options=varargin{5};
else
    weights=ones(length(fixedimage),1);
end

% Load preset parameter set
[~,funname]=fileparts(options.prefs.machine.normsettings.ants_preset);
apref=eval(funname);

directory=fileparts(movingimage{1});
directory=[directory,filesep];

for fi=1:length(fixedimage)
    fixedimage{fi} = ea_path_helper(ea_niigz(fixedimage{fi}));
end
for fi=1:length(movingimage)
    movingimage{fi} = ea_path_helper(ea_niigz(movingimage{fi}));
end

if length(fixedimage)~=length(movingimage)
    ea_error('Please supply pairs of moving and fixed images (can be repetitive).');
end

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
    [~, imgsize] = system(['bash -c "', HEADER, ' ',fixedimage{1}, ' 2"']);
else
    [~, imgsize] = system([HEADER, ' ', fixedimage{1}, ' 2']);
end

imgsize = cellfun(@(x) str2double(x),ea_strsplit(imgsize,'x'));

    rigidconvergence=apref.convergence.rigid;
    rigidshrinkfactors=apref.shrinkfactors.rigid;
    rigidsmoothingssigmas=apref.smoothingsigmas.rigid;
    
    affineconvergence=apref.convergence.affine;
    affineshrinkfactors=apref.shrinkfactors.affine;
    affinesmoothingssigmas=apref.smoothingsigmas.affine;
    
    synconvergence=apref.convergence.syn;
    synshrinkfactors=apref.shrinkfactors.syn;
    synsmoothingssigmas=apref.smoothingsigmas.syn;

rigidstage = [' --initial-moving-transform [', fixedimage{1}, ',', movingimage{1}, ',1]' ...
    ' --transform Rigid[0.1]' ...
    ' --convergence ', rigidconvergence, ...
    ' --shrink-factors ', rigidshrinkfactors, ...
    ' --smoothing-sigmas ', rigidsmoothingssigmas, ...
    ' --masks [NULL,NULL]'];

for fi=1:length(fixedimage)
    try
        rigidstage=[rigidstage,...
            ' --metric ',apref.metric,'[', fixedimage{fi}, ',', movingimage{fi}, ',',num2str(weights(fi)),apref.metricsuffix,']'];
    catch
        keyboard
    end
end

affinestage = [' --transform Affine[0.1]'...
    ' --convergence ', affineconvergence, ...
    ' --shrink-factors ', affineshrinkfactors ...
    ' --smoothing-sigmas ', affinesmoothingssigmas, ...
    ' --masks [NULL,NULL]'];

for fi=1:length(fixedimage)
    affinestage=[affinestage,...
        ' --metric ',apref.metric,'[', fixedimage{fi}, ',', movingimage{fi}, ',',num2str(weights(fi)),apref.metricsuffix,']'];
end



synstage = [' --transform ',apref.antsmode,apref.antsmode_suffix...
    ' --convergence ', synconvergence, ...
    ' --shrink-factors ', synshrinkfactors ...
    ' --smoothing-sigmas ', synsmoothingssigmas, ...
    ' --masks [NULL,NULL]'];


for fi=1:length(fixedimage)
    synstage=[synstage,...
        ' --metric ',apref.metric,'[', fixedimage{fi}, ',', movingimage{fi}, ',',num2str(weights(fi)),apref.metricsuffix,']'];
end

% add slab stage
if slabspresent
    slabstage=[' --transform ',apref.antsmode,apref.antsmode_suffix...
        ' --convergence ', synconvergence, ...
        ' --shrink-factors ', synshrinkfactors ...
        ' --smoothing-sigmas ', synsmoothingssigmas, ...
        ' --use-estimate-learning-rate-once ', ...
        ' --masks [NULL,',[tmaskdir,filesep,'slabmask_',slabID,'.nii'],']'];
    fixedimage=[fixedimage,slabfixedimage];
    movingimage=[movingimage,slabmovingimage];
    
    for fi=1:length(fixedimage)
        slabstage=[slabstage,...
            ' --metric ',apref.metric,'[', fixedimage{fi}, ',', movingimage{fi}, ',',num2str(weights(fi)),apref.metricsuffix,']'];
    end
else
    slabstage='';
end

% add subcortical refine stage:
if subcorticalrefine
    synmaskconvergence=apref.convergence.scrf;
    synmaskshrinkfactors=apref.shrinkfactors.scrf;
    synmasksmoothingssigmas=apref.smoothingsigmas.scrf;
    
    if slabspresent
        movingmask=[tmaskdir,filesep,'slabmask_',slabID,'.nii'];
    else
        movingmask='NULL';
    end
    synmaskstage = [' --transform ',apref.antsmode,apref.antsmode_suffix, ...
        ' --convergence ', synmaskconvergence, ...
        ' --shrink-factors ', synmaskshrinkfactors,  ...
        ' --smoothing-sigmas ', synmasksmoothingssigmas, ...
        ' --use-estimate-learning-rate-once ', ...
        ' --masks [',ea_space([],'subcortical'),'secondstepmask','.nii',',',movingmask,']'];
    for fi=1:length(fixedimage)
        synmaskstage=[synmaskstage,...
            ' --metric ',apref.metric,'[', fixedimage{fi}, ',', movingimage{fi}, ',',num2str(weights(fi)),apref.metricsuffix,']'];
    end
else
    synmaskstage='';
end

ea_libs_helper
if options.prefs.machine.normsettings.ants_numcores
    setenv('ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS',prefs.machine.normsettings.ants_numcores) % no num2str needed since stored as string.
end



