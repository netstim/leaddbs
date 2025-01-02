function modality = ea_getmodality(BIDSFilePath, opts)
% Extract image modality from BIDS file path
arguments
    BIDSFilePath {mustBeText}
    opts.acq {mustBeNumericOrLogical} = true % Keep acq label by default
end

if ~iscell(BIDSFilePath)
    wasChar = 1;
    BIDSFilePath = {BIDSFilePath};
else
    wasChar = 0;
end

modality = cell(size(BIDSFilePath));

for i=1:length(BIDSFilePath)
    if ~opts.acq || ~isempty(regexp(BIDSFilePath{i}, '_CT\.nii(.gz)?$', 'once')) % Skip plane label
        modality{i} = regexp(BIDSFilePath{i}, '(?<=_)([^\W_]+)(?=\.nii(\.gz)?$)', 'match', 'once');
    else % Keep plane label
        modality{i} = regexp(BIDSFilePath{i}, '(?<=_acq-)((ax|sag|cor|iso)\d*_[^\W_]+)(?=\.nii(\.gz)?$)', 'match', 'once');
    end
end

if wasChar
    modality = modality{1};
end
