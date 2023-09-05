function ea_runthomas_gui(handles)
% Run THOMAS from Tools menu.

bids = getappdata(handles.leadfigure, 'bids');
subjId = getappdata(handles.leadfigure, 'subjId');
if isempty(bids) || isempty(subjId)
    ea_cprintf('CmdWinWarnings', 'No patient selected!\n');
    return;
else
    subj = bids.getSubj(subjId{1});

    % Check T1w modality
    if contains(subj.AnchorModality, 'T1')
        modality = subj.AnchorModality;
    else
        preopModalities = fieldnames(subj.coreg.anat.preop);
        T1Ind = find(contains(preopModalities, 'T1'), 1);
        if isempty(T1Ind)
            ea_cprintf('CmdWinWarnings', 'Pre-op T1w image not found!\n');
            return;
        else
            modality = preopModalities{T1Ind};
        end
    end

    % Prepare file
    ea_mkdir(fullfile(subj.subjDir, 'thomas'));
    copyfile(subj.coreg.anat.preop.(modality), fullfile(subj.subjDir, 'thomas'));
    T1Image = replace(subj.coreg.anat.preop.(modality), fullfile(subj.coregDir, 'anat'), fullfile(subj.subjDir, 'thomas'));

    % Run THOMAS segmentation
    ea_runthomas(T1Image);

    % Copy generated atlas into subj atlas folder
    copyfile(fullfile(subj.subjDir, 'thomas', 'atlases', '*'), subj.atlasDir);
end
