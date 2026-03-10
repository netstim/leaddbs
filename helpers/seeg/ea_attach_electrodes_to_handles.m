
function handles = ea_attach_electrodes_to_handles(handles, subjDir)
    electrodes = ea_list_electrodes(subjDir);
    handles.currentPatient = struct();
    handles.currentPatient.subjDir   = subjDir;
    handles.currentPatient.electrodes = electrodes;

    % Populate a UI control if available (table or list)
    if isfield(handles,'electrode_list') && isgraphics(handles.electrode_list)
        labels = arrayfun(@(e) e.label, electrodes, 'uni', 0);
        if isempty(labels), labels = "none"; end
        try
            handles.electrode_list.Items = labels;
        catch
            set(handles.electrode_list,'String',labels,'Value',min(1,numel(labels)));
        end
    end

    % Optional: set model popup based on first electrode; keep editable
    if isfield(handles,'electrode_model_popup') && isgraphics(handles.electrode_model_popup)
        if ~isempty(electrodes) && ~isempty(electrodes(1).model)
            try
                val = find(strcmp(handles.electrode_model_popup.String, electrodes(1).model), 1);
                if ~isempty(val), handles.electrode_model_popup.Value = val; end
            end
        end
    end
end
