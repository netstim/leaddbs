% ===============================
% File: ea_get_selected_electrode.m
% Purpose: Utility for downstream modules to fetch the currently selected
% electrode regardless of how many exist.
% ===============================
function e = ea_get_selected_electrode(handles)
    e = struct();
    if ~isfield(handles,'currentPatient') || ~isfield(handles.currentPatient,'electrodes')
        return
    end
    E = handles.currentPatient.electrodes;
    sel = 1;
    if isfield(handles,'electrode_list') && isgraphics(handles.electrode_list)
        try
            sel = handles.electrode_list.Value;
        catch
            v = get(handles.electrode_list,'Value');
            if ~isempty(v), sel = v; end
        end
    end
    sel = max(1, min(sel, numel(E)));
    e = E(sel);
end