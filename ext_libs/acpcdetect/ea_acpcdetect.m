function ea_acpcdetect(inputimage, outputimage)
% Wrapper for acpcdetect
%
% Detect AC and PC locations, set the origin of the input NIfTI image to
% the AC-PC midpoint. Output image will be tilt-corrected and in RAS
% orientation.

if nargin < 2
    outputimage = inputimage;
end

ea_libs_helper;

basedir = fileparts(mfilename('fullpath'));

if ispc
    ACPCDETECT = ea_path_helper([basedir, filesep, 'acpcdetect.exe']);
else
    ACPCDETECT = [basedir, filesep, 'acpcdetect.', computer('arch')];
end

cmd = [ACPCDETECT, ' -v -no-tilt-correction -noppm -nopng -notxt -i ', ea_path_helper(inputimage)];

setenv('ARTHOME', basedir);
fprintf('\nacpcdetect ...\n\n');
if ~ispc
    [~, cmdout] = system(['bash -c "', cmd, '"']);
else
    [~, cmdout] = system(cmd);
end

disp(cmdout);

niipath = ea_niifileparts(inputimage);
ea_delete([niipath, '.mrx']);
ea_delete([niipath, '_FSL.mat']);
if isfile([niipath, '_RAS.nii'])
    movefile([niipath, '_RAS.nii'], outputimage);
else
    ea_cprintf('CmdWinWarnings', 'Failed to do acpcdetect!\n');
end
