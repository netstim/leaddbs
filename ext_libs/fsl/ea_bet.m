function ea_bet(inputimage, outputmask, outputimage, fraintthreshold)
% Wrapper for bet2 brain extraction

% Do not output brain mask by default
if nargin < 2
    outputmask = 0;
end

% overwrite the input image
if ~exist('outputimage','var')
    outputimage = inputimage;
else
    if isempty(outputimage)
        outputimage=inputimage;
    end
end

% If no threshold is set, the default fractional intensity threshold (0->1) is set
% to default=0.5 (as set in FSL binary), smaller values give larger brain outline estimates
if nargin < 4
    fraintthreshold = 0.5;
end

fprintf('\n\nRunning FSL BET2: %s\n\n', inputimage);

inputimage = ea_niifileparts(inputimage);
[outputimage, ~, ext] = ea_niifileparts(outputimage);

basedir = [fileparts(mfilename('fullpath')), filesep];
if ispc
    BET = ea_path_helper([basedir, 'bet2.exe']);
else
    BET = ea_path_helper([basedir, 'bet2.', computer('arch')]);
end

cmd = [BET, ...
       ' ', ea_path_helper(inputimage), ...
       ' ', ea_path_helper(outputimage), ...
       ' --verbose'];

if outputmask
    cmd = [cmd, ' -m'];
end

cmd = [cmd, ' -f ' ,num2str(fraintthreshold)];
switch ext
    case '.nii'
        setenv('FSLOUTPUTTYPE','NIFTI');
    case '.nii.gz'
        setenv('FSLOUTPUTTYPE','NIFTI_GZ');
end
if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end

fprintf('\nFSL BET2 finished\n');
