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

inputimage = ea_path_helper(ea_niigz(inputimage));
outputimage = ea_path_helper(ea_niigz(outputimage));

% Remove the '.nii' or '.nii.gz' ext
outputimage = ea_niifileparts(outputimage);

basedir = [fileparts(mfilename('fullpath')), filesep];
if ispc
    BET = [basedir, 'bet2.exe'];
else 
    BET = [basedir, 'bet2.', computer('arch')];    
end

cmd = [BET, ...
       ' ', inputimage, ...
       ' ', outputimage, ...
       ' --verbose'];

if outputmask
    cmd = [cmd, ' -m'];
end

cmd = [cmd, ' -f ' ,num2str(fraintthreshold)];

setenv('FSLOUTPUTTYPE','NIFTI');
if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);    
end
