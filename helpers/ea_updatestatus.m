function ea_updatestatus(handles, subj)
% subj is the struct returned by BIDSFetcher.getSubj(subjId)

if isfile(fullfile(fileparts(fileparts(fileparts(subj.subjDir))), 'miniset.json'))
    isMiniset = 1;
else
    isMiniset = 0;
end

if strcmp(subj.postopModality, 'CT')
    if isMiniset
        statusone = 'Miniset is used.';
    elseif isfile(subj.postopAnat.CT.coreg)
        statusone = 'Coregistered post-op CT found. Please run normalization.';
    elseif isfile(subj.postopAnat.CT.raw)
        statusone = 'Raw post-op CT found. Please coregister to pre-op MRI first.';
    else
        statusone = 'Seems no post-op images found.';
    end
elseif strcmp(subj.postopModality, 'MRI')
    fields = fieldnames(subj.postopAnat);
    if isMiniset
        statusone = 'Miniset is used.';
    elseif isfile(subj.postopAnat.(fields{1}).coreg)
        statusone = 'Coregistered post-op MRI found. Please run normalization.';
    elseif isfile(subj.postopAnat.(fields{1}).raw)
        statusone = 'Raw post-op MRI found. Please coregister to pre-op MRI first.';
    else
        statusone = 'Seems no post-op images found.';
    end
elseif strcmp(subj.postopModality, 'None')
    statusone = 'No post-op images found.';
end

if isfield(subj, 'preopAnat')
    fields = fieldnames(subj.preopAnat);
    if isfile(subj.preopAnat.(fields{1}).norm) ...
            && ~isempty(dir([subj.norm.transform.forwardBaseName, '*'])) ...
            && ~isempty(dir([subj.norm.transform.inverseBaseName, '*']))
        statusone = 'Normalization has been done.';
    end
end

set(handles.statusone, 'String', statusone);
set(handles.statusone, 'TooltipString', statusone);

if strcmp(subj.postopModality, 'None')
    statustwo = '';
elseif isfile(subj.recon.recon)
    statustwo = 'Reconstruction found.';
else
    statustwo = 'Reconstruction not found. Please run reconstruction.';
end

set(handles.statustwo, 'String', statustwo);
set(handles.statustwo, 'TooltipString', statustwo);
