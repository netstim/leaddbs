function ea_options2handles(options,handles)

if  isempty(handles.patientlist.Data) || isempty(handles.patientlist.Selection) || ~isfield(options, 'modality')
    arrayfun(@(x) set(x, 'Enable', 'off'), handles.registrationtab.Children);
    arrayfun(@(x) set(x, 'Enable', 'off'), handles.localizationtab.Children);
    arrayfun(@(x) set(x, 'Enable', 'off'), handles.optionaltab.Children);
    set(handles.overwriteapproved, 'Enable', 'off');
elseif options.modality == 3
    arrayfun(@(x) set(x, 'Enable', 'on'), handles.registrationtab.Children);
    set(handles.coregctmethod, 'Enable', 'off');
    set(handles.scrf, 'Enable', 'off');
    set(handles.scrfmask, 'Enable', 'off');
    arrayfun(@(x) set(x, 'Enable', 'off'), handles.localizationtab.Children);
    arrayfun(@(x) set(x, 'Enable', 'on'), handles.optionaltab.Children);
    set(handles.overwriteapproved, 'Enable', 'on');
else
    arrayfun(@(x) set(x, 'Enable', 'on'), handles.registrationtab.Children);
    arrayfun(@(x) set(x, 'Enable', 'on'), handles.localizationtab.Children);
    arrayfun(@(x) set(x, 'Enable', 'on'), handles.optionaltab.Children);
    set(handles.overwriteapproved, 'Enable', 'on');
end

if isfield(options, 'coregct')
    set(handles.coreg_checkbox,'Value',options.coregct.do);
    
    if ~isfield(options.coregct, 'method')
        index = [];
    else
        index = find(ismember(handles.coregctmethod.String, options.coregct.method), 1);
    end
    
    if ~isempty(index)
        set(handles.coregctmethod, 'Value', index);
    else
        set(handles.coregctmethod, 'Value', 1);
    end
end

if isfield(options, 'normalize')
    set(handles.normalize_checkbox,'Value',options.normalize.do);

    if ~isfield(options.normalize, 'method')
        index = [];
    else
        index = find(ismember(handles.normmethod.String, options.normalize.method), 1);
    end

    if ~isempty(index)
        set(handles.normmethod,'Value',index);
    else
        set(handles.normmethod,'Value',1);
    end
end

if isfield(options, 'checkreg')
    set(handles.checkreg, 'Value', options.checkreg);
else
    set(handles.checkreg, 'Value', 0);
end


if isfield(options, 'normcheck')
    set(handles.normcheck, 'Value', options.normcheck);
end

if isfield(options, 'modality')
    set(handles.MRCT,'Value',options.modality);
end

if isfield(options, 'sides')
    for i=1:15
        if ismember(i,options.sides)
            set(handles.(['side', num2str(i)]), 'Value', 1);
        else
            set(handles.(['side', num2str(i)]), 'Value', 0);
        end
    end
end

if isfield(options, 'doreconstruction')
    set(handles.doreconstruction_checkbox,'Value',options.doreconstruction);
end

if isfield(options, 'automask') && options.automask
    set(handles.maskwindow_txt,'String','auto')
elseif isfield(options, 'maskwindow') && ~isempty(options.maskwindow)
    set(handles.maskwindow_txt,'String',num2str(options.maskwindow));
end

if isfield(options, 'd2')
    set(handles.writeout2d_checkbox,'Value',options.d2.write);
end

if isfield(options, 'd3')
    set(handles.render_checkbox,'Value',options.d3.write);
    set(handles.exportservercheck,'Value',options.d3.autoserver);
end

if isfield(options, 'manualheightcorrection')
    set(handles.manualheight_checkbox,'Value',options.manualheightcorrection);
end

if isfield(options, 'entrypoint')
    set(handles.targetpopup, 'Value', find(ismember(handles.targetpopup.String, options.entrypoint)));
else
    set(handles.targetpopup, 'Value', 1);
end

if isfield(options, 'elmodel')
    set(handles.electrode_model_popup, 'Value', find(ismember(handles.electrode_model_popup.String, options.elmodel)));
end

if isfield(options, 'atlasset')
    set(handles.atlassetpopup, 'Value',  find(ismember(handles.atlassetpopup.String, options.atlasset)));
end
