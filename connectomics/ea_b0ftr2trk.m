function ea_b0ftr2trk(ftrfile, b0image)
% Convert FTR in b0 space to TRK
%
% Handle the orientation problem, thus the orientation of the TRK will be
% the same as the b0 image.

[directory, ftrname, ext] = fileparts(ftrfile);
if isempty(directory)
    directory = '.';
end
if isempty(ext)
    ftrfile = [ftrfile, '.mat'];
end

b0 = ea_load_nii(b0image);
niisize = size(b0.img);
specs.origin = [0,0,0];
specs.dim = niisize;
specs.affine = b0.mat;

disp('Checking orientation...');
orient = ea_aff2axcodes(specs.affine);
disp(['b0 orientation: ', orient]);

ftr = load(ftrfile);
if ~isfield(ftr, 'voxmm')
    disp('Fixing voxmm field...');
    voxmm = 'vox';
    save(ftrfile, 'voxmm', '-append');
    ftr.voxmm = 'vox';
end

flipped = [0 0 0];
% FTR is in RAS orientation, adapt it to b0 orientation.
if ~strcmp(orient(1), 'R')
    % Flip back X to L if needed.
    disp('Flip back X to L...');
    flipped(1) = 1;
    ftr.fibers(:,1) = b0.dim(1)-ftr.fibers(:,1);
end
if ~strcmp(orient(2), 'A')
    % Flip back Y to P if needed.
    disp('Flip back Y to P...');
    flipped(2) = 1;
    ftr.fibers(:,2) = b0.dim(2)-ftr.fibers(:,2);
end
if ~strcmp(orient(3), 'S')
    % Flip back X to I if needed.
    disp('Flip back Z to I...');
    flipped(3) = 1;
    ftr.fibers(:,3) = b0.dim(3)-ftr.fibers(:,3);
end

if any(flipped)
    ftr.fibers = single(ftr.fibers);
    disp('Save flipped FTR...');
    save(fullfile(directory,['b0', ftrname, '.mat']), '-struct', 'ftr', '-v7.3');

    disp('Start FTR to TRK conversion...');
    ea_ftr2trk(fullfile(directory,['b0',ftrname,'.mat']), specs);
    ea_delete(fullfile(directory,['b0',ftrname,'.mat']));
    movefile(fullfile(directory,['b0',ftrname,'.trk']), fullfile(directory,[ftrname,'.trk']))
else
    disp('Start FTR to TRK conversion...');
    ea_ftr2trk(fullfile(directory,ftrfile), specs);
end
