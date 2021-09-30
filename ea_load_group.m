function ea_load_group(handles, groupdir)

ea_busyaction('on', handles.leadfigure, 'group');

set(handles.groupdir_choosebox, 'String', groupdir);
set(handles.groupdir_choosebox, 'TooltipString', groupdir);

analysisFile = ea_getGroupAnalysisFile(folder, groupdir);
load(analysisFile, 'M')
if ~isfield(M.ui, 'mirrorsides')
    % Fix missing 'mirrorsides' field for old analysis
    try
        currentUISetting = getappdata(handles.leadfigure, 'M');
        M.ui.mirrorsides = currentUISetting.ui.mirrorsides;
    catch
        M.ui.mirrorsides = 0;
    end
end

setappdata(handles.leadfigure, 'M', M);
try
    setappdata(handles.leadfigure, 'S', M.S);
    setappdata(handles.leadfigure, 'vatmodel', M.S(1).model);
end

ea_refresh_lg(handles);

ea_addrecentpatient(handles, {groupdir}, 'groups', 'groups')

ea_busyaction('off', handles.leadfigure, 'group');
