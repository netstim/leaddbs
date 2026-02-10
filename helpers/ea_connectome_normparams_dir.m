function outDir = ea_connectome_normparams_dir(directory)
% Return directory where Lead-Connectome SPM norm params (y_ea_*.nii) are stored.
% BIDS: normalization/transformations; legacy: subject root.
if nargin < 1 || isempty(directory)
    outDir = '';
    return;
end
if contains(directory, 'derivatives') || contains(directory, 'leaddbs')
    outDir = fullfile(directory, 'normalization', 'transformations');
else
    outDir = directory;
end
