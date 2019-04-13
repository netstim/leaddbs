function ea_b0ftr2trk(FTRfile, b0image)
% Convert FTR in b0 space to TRK
%
% Handle the orientation problem, thus the orientation of the TRK will be
% the same as the b0 image.

b0 = ea_load_nii(b0image);
niisize = size(b0.img);
specs.origin = [0,0,0];
specs.dim = niisize;
specs.affine = b0.mat;

disp('Checking orientation...');
orient = ea_aff2axcodes(specs.affine);
disp(['b0 orientation: ', orient]);

ftr = load(FTRfile);
if ~isfield(ftr, 'voxmm')
    disp('Fixing voxmm field...');
    voxmm = 'vox';
    save(FTRfile, 'voxmm', '-append');
    ftr.voxmm = 'vox';
end

% FTR is in RAS orientation, adapt it to b0 orientation.
if ~strcmp(orient(1), 'R')
    % Flip back X to L if needed.
    disp('Flip back X to L...');
    ftr.fibers(:,1) = b0.dim(1)-ftr.fibers(:,1);
end
if ~strcmp(orient(2), 'A')
    % Flip back Y to P if needed.
    disp('Flip back Y to P...');
    ftr.fibers(:,2) = b0.dim(2)-ftr.fibers(:,2);
end
if ~strcmp(orient(3), 'S')
    % Flip back X to I if needed.
    disp('Flip back Z to I...');
    ftr.fibers(:,3) = b0.dim(3)-ftr.fibers(:,3);
end

[path, name] = fileparts(FTRfile);
if isempty(path)
    path = '.';
end

save(fullfile(path, ['b0', name, '.mat']), '-struct', 'ftr', '-v7.3');
ea_ftr2trk(['b0', name], path, specs);

ea_delete(fullfile(path, ['b0', name, '.mat']));
movefile(fullfile(path, ['b0', name, '.trk']), fullfile(path, [name, '.trk']))
