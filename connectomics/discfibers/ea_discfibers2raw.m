function ea_discfibers2raw(discfiber, patientdir)
% Warp discriminative fibers to patient T1 (anchor in LeadDBS) and raw T1 space
%
%     discfiber: discriminative fibers including fiber cell and T-scores
%     patientdir: patient directory

% Load discfibers in MNI space
load(discfiber, 'fibcell')
load(discfiber, 'vals')

if size(fibcell,2) == 2
    fibcell = vertcat(fibcell{:});
end

if iscell(vals) && size(vals,2) == 2
    vals = vertcat(vals{:});
end

% Set reference images
refMNI = [ea_space, 't1.nii'];
refAnchor = fullfile(patientdir, 'anat_t1.nii');
refRawAnchor = fullfile(patientdir, 'raw_anat_t1.nii');

% Calculate fibcell indices
idx = cellfun(@(x) size(x,1), fibcell);
endIndex = cumsum(idx);
startIndex = [1;endIndex+1];
startIndex = startIndex(1:end-1);
fibIndex = arrayfun(@(x,y) [x:y]', startIndex, endIndex, 'Uni', 0);

fibersMNIVox = ea_mm2vox(cell2mat(fibcell), refMNI)';
% Convert discfibers from MNI space to anchor space
[fibersAnchor, fibersAnchorVox]  = ea_map_coords(fibersMNIVox, ...
    refMNI, ...
    fullfile(patientdir,'forwardTransform'), ...
    refAnchor);

% Convert fiber matrix to fibcell
fibersAnchor = single(fibersAnchor)';
fibcellAnchor = cellfun(@(x) fibersAnchor(x,:), fibIndex, 'Uni', 0);

% Save discfibers in anchor space
fibcell = fibcellAnchor;
outputMAT = regexprep(discfiber, '\.mat$', 'Anchor.mat');
save(outputMAT, 'fibcell', 'opts', 'vals');

% Register anchor to raw anchor
options = ea_getptopts(patientdir);
options.coregmr.method = 'SPM (Friston 2007)';
affinefile = ea_coregimages(options, refAnchor, refRawAnchor, [patientdir,'tmp.nii'],{},1);
ea_delete([patientdir,'tmp.nii']);

% Convert discfibers from achor space to raw anchor space
fibersRawAnchor = ea_map_coords(fibersAnchorVox, ...
    refAnchor, ...
    affinefile{1}, ...
    refRawAnchor, ...
    options.coregmr.method);

% Convert fiber matrix to fibcell
fibersRawAnchor = single(fibersRawAnchor)';
fibcellRawAnchor = cellfun(@(x) fibersRawAnchor(x,:), fibIndex, 'Uni', 0);

% Save discfibers in anchor space
fibcell = fibcellRawAnchor;
outputMAT = regexprep(discfiber, '\.mat$', 'RawAnchor.mat');
save(outputMAT, 'fibcell', 'vals');
