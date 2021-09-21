function modality = ea_getmodality(BIDSFilePath)
% Extract image modality from BIDS file path and keep [ax|cor|sag] label
% for post-op MRI

if isempty(regexp(BIDSFilePath, '_acq-(ax|cor|sag)_', 'once'))
    modality = regexp(BIDSFilePath, '(?<=_)([^\W_]+)(?=\.nii(\.gz)?$)', 'match', 'once');
else % Keep plane label for post-op MRI
    modality = regexp(BIDSFilePath, '(?<=_acq-)((ax|sag|cor)_[^\W_]+)(?=\.nii(\.gz)?$)', 'match', 'once');
end
