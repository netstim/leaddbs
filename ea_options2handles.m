function ea_options2handles(options,handles)

% set handles
set(handles.dicomcheck,'Value',options.dicomimp.do);
set(handles.normalize_checkbox,'Value',options.normalize.do);
if options.normalize.methodn>length(handles.normmethod,'String')
    set(handles.normmethod,'Value',1);
else
    set(handles.normmethod,'Value',options.normalize.methodn);
end
set(handles.normcheck,'Value',options.normalize.check);

% CT coregistration
set(handles.coregct_checkbox,'Value',options.coregct.do);
if options.coregct.methodn>length(handles.coregctmethod,'String')
    set(handles.coregctmethod,'Value',1);
else
    set(handles.coregctmethod,'Value',options.coregct.methodn);
end
set(handles.coregctcheck,'Value',options.coregctcheck);

set(handles.MRCT,'Value',options.modality);

if ismember(1,options.sides)
    set(handles.right_checkbox,'Value',1);
else
        set(handles.right_checkbox,'Value',0);
end
if ismember(2,options.sides)
    set(handles.left_checkbox,'Value',1);
else
    set(handles.left_checkbox,'Value',0);
end

set(handles.doreconstruction_checkbox,'Value',options.doreconstruction);

if options.automask
    set(handles.maskwindow_txt,'String','auto')
else
    set(handles.maskwindow_txt,'String',num2str(options.maskwindow));
end
set(handles.genptatlascheck,'Value',options.atl.genpt); % generate patient specific atlases
set(handles.writeout2d_checkbox,'Value',options.d2.write);
set(handles.tdcolorscheck,'Value',options.d2.col_overlay);
set(handles.tdcontourcheck,'Value',options.d2.con_overlay);
setappdata(handles.tdcontourcolor,'color',options.d2.con_color);
set(handles.tdlabelcheck,'Value',options.d2.lab_overlay);
set(handles.bbsize,'String',num2str(options.d2.bbsize));
set(handles.manualheight_checkbox,'Value',options.manualheightcorrection);
set(handles.render_checkbox,'Value',options.d3.write);
set(handles.targetpopup,'Value',options.entrypointn);
set(handles.electrode_model_popup,'Value',options.elmodeln);
set(handles.atlassetpopup,'Value',options.atlassetn);
set(handles.exportservercheck,'Value',options.d3.autoserver);
