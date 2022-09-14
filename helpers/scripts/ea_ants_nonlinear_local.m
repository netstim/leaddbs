function ea_ants_nonlinear_local(varargin)
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

if nargin >= 6
    interp = varargin{6};
else
    interp = 'Linear';
end

[outputdir, outputname, ~] = fileparts(outputimage);
if outputdir
    outputbase = [outputdir, filesep, outputname];
else
    outputbase = ['.', filesep, outputname];
end

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

basedir = [ea_getearoot,'ext_libs',filesep,'ANTs',filesep];

if ispc
    HEADER = ea_path_helper([basedir, 'PrintHeader.exe']);
    ANTS = ea_path_helper([basedir, 'antsRegistration.exe']);
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

convergence='[50x25x10x5,1e-7,10]';
shrinkfactors='12x8x4x1';
smoothingssigmas='4x4x2x1vox';

mask=ea_load_nii(movingimage{1});
mask.img=mask.img~=0;
mask.fname=([tempdir,'synmask.nii']);
ea_write_nii(mask);

synstage = [' --transform BSplineSyN[0.1,26,0,3]'...
            ' --convergence ', convergence, ...
            ' --shrink-factors ', shrinkfactors ...
            ' --smoothing-sigmas ', smoothingssigmas ...
            ' --masks [NULL,',mask.fname,']'];

        for fi=1:length(fixedimage)

            suffx=',4,Random,0.5';
            synstage=[synstage,...
                ' --metric ',metrics{fi},'[', fixedimage{fi}, ',', movingimage{fi}, ',',num2str(weights(fi)),suffx,']'];
        end

ea_libs_helper;

cmd = [ANTS, ' --verbose 1' ...
             ' --dimensionality 3 --float 1' ...
             ' --output [',ea_path_helper(outputbase), ',', outputimage, ']' ...
             ' --interpolation ',interp ...
             ' --write-composite-transform 1', ...
             synstage];

if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end

ea_conv_antswarps(directory);
