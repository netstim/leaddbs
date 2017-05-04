function ea_ants_nonlinear(varargin)
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
subcorticalrefine=varargin{7};
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
    disp(['Checking for slabs among structural images (assuming dominant structural file ',movingimage{end},'is a whole-brain acquisition)...']);
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
            mnii.dt=[4,0];
            mnii.img=AllMX;
            mnii.fname=[tmaskdir,filesep,'slabmask.nii'];
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
    metrics=varargin{5};
    options=varargin{6};
else
    weights=ones(length(fixedimage),1);
    metrics=repmat({'MI'},length(fixedimage),1);
end
useSyN=options.prefs.machine.normsettings.ants_synmode;


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
        applyTransforms = [basedir, 'antsApplyTransforms.exe'];

else
    HEADER = [basedir, 'PrintHeader.', computer('arch')];
    ANTS = [basedir, 'antsRegistration.', computer('arch')];
        applyTransforms = [basedir, 'antsApplyTransforms.', computer('arch')];

end



if ~ispc
    [~, imgsize] = system(['bash -c "', HEADER, ' ',fixedimage{1}, ' 2"']);
else
    [~, imgsize] = system([HEADER, ' ', fixedimage{1}, ' 2']);
end

imgsize = cellfun(@(x) str2double(x),ea_strsplit(imgsize,'x'));

if any(imgsize>256)
    rigidshrinkfactors='12x8x4x2';
    affineshrinkfactors='12x8x4x2';
    synshrinkfactors='12x8x4x2';
else
    rigidshrinkfactors='8x4x2x1';
    affineshrinkfactors='8x4x2x1';
    synshrinkfactors='6x4x2x1';
end
rigidconvergence='[1000x500x250x100,1e-6,10]';
rigidsmoothingssigmas='3x2x1x0';

affineconvergence='[1000x500x250x100,1e-6,10]';
affinesmoothingssigmas='3x2x1x0';

synconvergence='[100x70x50x20,1e-6,10]';
synsmoothingssigmas='3x2x1x0';


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

switch useSyN
    case 'SyN'
        tsuffix='[0.1,3.0,0.0]';
    case 'BSplineSyN'
        tsuffix='[0.1,26,0,3]'; % as in example script in Tustison 2013
    otherwise
        tsuffix='[0.1,26,0,3]'; % need to find proper values.       
end


synstage = [' --transform ',useSyN,tsuffix...
    ' --convergence ', synconvergence, ...
    ' --shrink-factors ', synshrinkfactors ...
    ' --smoothing-sigmas ', synsmoothingssigmas, ...
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
    synstage=[synstage,...
        ' --metric ',metrics{fi},'[', fixedimage{fi}, ',', movingimage{fi}, ',',num2str(weights(fi)),suffx,']'];
end




% add slab stage


if slabspresent
    slabstage=[' --transform ',useSyN,'[0.3,3.0,0.0]'...
        ' --convergence ', synconvergence, ...
        ' --shrink-factors ', synshrinkfactors ...
        ' --smoothing-sigmas ', synsmoothingssigmas, ...
        ' --masks [NULL,',[tmaskdir,filesep,'slabmask.nii'],']'];
    fixedimage=[fixedimage,slabfixedimage];
    movingimage=[movingimage,slabmovingimage];
    
    for fi=1:length(fixedimage)
        switch metrics{fi}
            case 'MI'
                suffx=',32,Regular,0.25';
            case 'CC'
                suffx=',4';
            case 'GC'
                suffx=',15,Random,0.05';
        end
        slabstage=[slabstage,...
            ' --metric ',metrics{fi},'[', fixedimage{fi}, ',', movingimage{fi}, ',',num2str(weights(fi)),suffx,']'];
    end
    
else
    slabstage='';
end




% add subcortical refine stage:
if subcorticalrefine
    synmaskconvergence='[20x5,1e-6,10]';
    synmaskshrinkfactors='2x1';
    synmasksmoothingssigmas='1x0vox';
    
    if slabspresent
        movingmask=[tmaskdir,filesep,'slabmask.nii'];
    else
        movingmask='NULL';
    end
    synmaskstage = [' --transform ',useSyN,'[0.2,3.0,0.0]', ...
        ' --convergence ', synmaskconvergence, ...
        ' --shrink-factors ', synmaskshrinkfactors,  ...
        ' --smoothing-sigmas ', synmasksmoothingssigmas, ...
        ' --masks [',ea_space([],'subcortical'),'secondstepmask','.nii',',',movingmask,']'];
    for fi=1:length(fixedimage)
        switch metrics{fi}
            case 'MI'
                suffx=',32,Regular,0.25';
            case 'CC'
                suffx=',4';
            case 'GC'
                suffx=',15,Random,0.05';
        end
        synmaskstage=[synmaskstage,...
            ' --metric ',metrics{fi},'[', fixedimage{fi}, ',', movingimage{fi}, ',',num2str(weights(fi)),suffx,']'];
    end
else
    synmaskstage=''; 
end

ea_libs_helper
%setenv('ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS','8')




cmd = [ANTS, ' --verbose 1', ...
    ' --dimensionality 3', ...
    ' --output [',ea_path_helper(outputbase), ',', outputimage, ']', ...
    ' --interpolation Linear', ...
    ' --use-histogram-matching 1', ...
    ' --float 1',...
    ' --write-composite-transform 1', ...
    rigidstage, affinestage, synstage, slabstage, synmaskstage];

display(cmd)
fid=fopen([directory,'ea_ants_command.txt'],'a');
fprintf(fid,[datestr(datetime('now')),':\n',cmd,'\n\n']);
fclose(fid);

if ~ispc
    system(['bash -c "', cmd, '"']);

else
    system(cmd);
end


ea_conv_antswarps(directory);
