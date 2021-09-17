function modality = ea_getmodality(BIDSFilePath)
% Extract image modality from BIDS file path and keep [ax|cor|sag] label
% for post-op MRI

if isempty(regexp(BIDSFilePath, '_acq-(ax|cor|sag)_', 'once'))
    modality = regexp(BIDSFilePath, '(?<=_)([a-zA-Z0-9]+)(?=\.nii(\.gz)?$)', 'match', 'once');
else % Keep plane label for post-op MRI
    modality = regexp(BIDSFilePath, '(?<=_acq-)((ax|sag|cor)_[a-zA-Z0-9]+)(?=\.nii(\.gz)?$)', 'match', 'once');
end
