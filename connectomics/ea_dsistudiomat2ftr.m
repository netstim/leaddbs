function ea_dsistudiomat2ftr(matFile, refb0, outputName, outputTRK, LPS)
% Convert the mat file exported from DSI-Studio to the FTR format in LEAD
%
% If the trk is going to be visualized in Surf-Ice, LPS should be set to 1
% to fix the orientation.

[directory, matName] = fileparts(matFile);
if isempty(directory)
    directory = '.';
end

if ~exist('outputname', 'var') || isempty(outputName)
    outputName = matName;
else
    outputName = strrep(outputName, '.mat', '');
end

if ~exist('outputtrk', 'var')
    outputTRK = 1;
end

if ~exist('LPS', 'var')
    LPS = 0;
end

fibinfo = load(matFile);
fibers = fibinfo.tracts';
idx = double(fibinfo.length)';
fibers = [fibers, repelem(1:numel(idx), idx)'];
clear fibinfo
b0 = spm_vol(refb0);

% Default orientation in DSI-Studio and TrackVis is LPS. Flip the
% coordinates to make the orientation in the MAT file inline with b0 image.
if b0.mat(1)>0  % 'R' is positive x-axis
    % flip x
    disp('Flip positive X-axis to R...');
    fibers(:,1) = b0.dim(1)-1-fibers(:,1);
end
if b0.mat(6)>0  % 'A' is positive y-axis
    %flip y
    disp('Flip positive Y-axis to A...');
    fibers(:,2) = b0.dim(2)-1-fibers(:,2);
end
if b0.mat(11)<0  % 'I' is positive z-axis
    %flip z
    disp('Flip positive Z-axis to I...');
    fibers(:,3) = b0.dim(3)-1-fibers(:,3);
end

% Change ZERO-BASED indexing to ONE-BASED indexing.
fibers(:,1:3) = fibers(:,1:3) + 1;

ftr.ea_fibformat = '1.0';
ftr.fourindex = 1;
ftr.fibers = fibers;
ftr.idx = idx;
ftr.voxmm = 'vox';
ftr.mat = b0.mat;

fprintf('\nSaving fibers...\n');
save([directory, filesep, outputName, '.mat'], '-struct', 'ftr', '-v7.3');
disp('Done.');

if outputTRK
    fprintf('\nGenerating trk in b0 space...\n');
    ea_ftr2trk([directory, filesep, outputName, '.mat'], refb0, LPS)
    disp('Done.');
end
