function ea_bet(inputimage, outputmask, outputimage)
% Wrapper for bet2 brain extraction

% Do not output brain mask by default
if nargin < 2
    outputmask = 0;
end

% overwrite the input image
if nargin < 3
    outputimage = inputimage;
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

setenv('FSLOUTPUTTYPE','NIFTI');
if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);    
end
