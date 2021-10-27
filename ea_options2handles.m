function ea_options2handles(options,handles)

% set handles
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

if isfield(options, 'checkreg')
    set(handles.checkreg, 'Value', options.checkreg);
else
    set(handles.checkreg, 'Value', 0);
end

% CT coregistration
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

if isfield(options, 'normcheck')
    set(handles.normcheck, 'Value', options.normcheck);
end

set(handles.MRCT,'Value',options.modality);

for i=1:15
    if ismember(i,options.sides)
        set(handles.(['side', num2str(i)]), 'Value', 1);
    else
        set(handles.(['side', num2str(i)]), 'Value', 0);
    end
end

set(handles.doreconstruction_checkbox,'Value',options.doreconstruction);

if isfield(options, 'automask') && options.automask
    set(handles.maskwindow_txt,'String','auto')
elseif isfield(options, 'maskwindow') && ~isempty(options.maskwindow)
    set(handles.maskwindow_txt,'String',num2str(options.maskwindow));
end

set(handles.writeout2d_checkbox,'Value',options.d2.write);
set(handles.render_checkbox,'Value',options.d3.write);
set(handles.exportservercheck,'Value',options.d3.autoserver);

if isfield(options, 'manualheightcorrection')
    set(handles.manualheight_checkbox,'Value',options.manualheightcorrection);
end

if isfield(options, 'entrypointn')
    set(handles.targetpopup,'Value',options.entrypointn);
end

if isfield(options, 'elmodeln')
    set(handles.electrode_model_popup,'Value',options.elmodeln);
end

if isfield(options, 'atlassetn')
    set(handles.atlassetpopup,'Value',options.atlassetn);
end
