function ea_dbsegment_menu(handles)
% Run DBSegment from Tools menu.

bids = getappdata(handles.leadfigure, 'bids');
subjId = getappdata(handles.leadfigure, 'subjId');
if isempty(bids) || isempty(subjId)
    ea_cprintf('CmdWinWarnings', 'No patient selected!\n');
    return;
else
    subj = bids.getSubj(subjId{1});

    % Choose modality to be segmented
    preopModalities = fieldnames(subj.coreg.anat.preop);
    preopModalities = preopModalities(endsWith(preopModalities, {'T1w'}));
    if isempty(preopModalities)
        ea_cprintf('CmdWinWarnings', 'Pre-op T1w image not found!\n');
        return;
    elseif isscalar(preopModalities)
        modality = preopModalities{1};
    else
        modality = questdlg('Choose pre-op image to be segmented', '', preopModalities{:}, preopModalities{1});
    end
    ea_cprintf('*Comments', 'Segmenting pre-op %s image...\n', regexprep(modality, '[^\W_]+_', ''));

    % Run DBSegment segmentation
    ea_dbsegment(subj.coreg.anat.preop.(modality));
end
