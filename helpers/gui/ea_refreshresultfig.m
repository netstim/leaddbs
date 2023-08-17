function ea_refreshresultfig(handles,refresh)
% this part makes changes of the figure active:
if isa(handles,'matlab.ui.Figure') % directly supplied figure handle - this way used when calling programmatically
   resultfig=handles; 
else % this way used when called from anatomycontrol GUI
    resultfig=getappdata(handles.acontrolfig,'resultfig');
end
try
    togglestates=getappdata(resultfig,'togglestates');
catch
    resultfig=gcf;
    togglestates=getappdata(handles.acontrolfig,'togglestates');
end
if ~isfield(togglestates,'refreshcuts')
    togglestates.refreshcuts=0;
end

if exist('refresh','var')
    togglestates.refreshview=refresh;
elseif isfield(togglestates,'refreshview')
else
    togglestates.refreshview=0;
end

% reset states based on gui:
if ~isa(handles,'matlab.ui.Figure')
    togglestates.xyzmm=[str2double(get(handles.xval,'String')),str2double(get(handles.yval,'String')),str2double(get(handles.zval,'String'))];
    togglestates.xyztoggles=[get(handles.xtoggle,'Value'),get(handles.ytoggle,'Value'),get(handles.ztoggle,'Value')];
    togglestates.xyztransparencies=[str2double(get(handles.xtrans,'String')),str2double(get(handles.ytrans,'String')),str2double(get(handles.ztrans,'String'))];
    togglestates.template=get(handles.templatepopup,'String');
    togglestates.template=togglestates.template{get(handles.templatepopup,'Value')};
    togglestates.tinvert=0;
    togglestates.customfile=getappdata(resultfig,'customfile');
end
setappdata(resultfig,'togglestates',togglestates); % also store toggle data in resultfig.

try
    options=getappdata(handles.acontrolfig,'options');
catch
    options=struct;
end
ea_anatomyslices(resultfig, togglestates, options, handles);

nativemni = ea_getnativemni;
if togglestates.refreshcuts && nativemni==1
    axis([-100 100 -130 100 -70 100]);
elseif togglestates.refreshcuts && nativemni==2
    axis([-200 200 -50 250 -100 150]);
end

togglestates.refreshcuts=0;
togglestates.refreshview=0;
setappdata(resultfig,'togglestates',togglestates);