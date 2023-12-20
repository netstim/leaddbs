function ea_thomas_menu(handles)
% Run THOMAS from Tools menu.

bids = getappdata(handles.leadfigure, 'bids');
subjId = getappdata(handles.leadfigure, 'subjId');
if isempty(bids) || isempty(subjId)
    ea_cprintf('CmdWinWarnings', 'No patient selected!\n');
    return;
else
    subj = bids.getSubj(subjId{1});

    % Choose modality to be segmented
    preopModalities = fieldnames(subj.coreg.anat.preop);
    preopModalities = preopModalities(endsWith(preopModalities, {'T1w', 'FGATIR', 'WMn'}));
    if isempty(preopModalities)
        ea_cprintf('CmdWinWarnings', 'Pre-op T1w/FGATIR/WMn image not found!\n');
        return;
    elseif numel(preopModalities) == 1
        modality = preopModalities{1};
    else
        modality = questdlg('Choose pre-op image to be segmented', '', preopModalities{:}, preopModalities{1});
    end
    ea_cprintf('*Comments', 'Segment pre-op %s image...\n', regexprep(modality, '[^\W_]+_', ''));

    % Prepare file
    ea_mkdir(fullfile(subj.subjDir, 'thomas'));
    copyfile(subj.coreg.anat.preop.(modality), fullfile(subj.subjDir, 'thomas'));
    T1Image = replace(subj.coreg.anat.preop.(modality), fullfile(subj.coregDir, 'anat'), fullfile(subj.subjDir, 'thomas'));

    % Run THOMAS segmentation
    if contains(modality, 'T1w')
        ea_thomas(T1Image, 't1');
    else
        ea_thomas(T1Image);
    end

    % Copy generated atlas into subj atlas folder
    copyfile(fullfile(subj.subjDir, 'thomas', 'atlases', '*'), subj.atlasDir);
end
