function ea_switchnormmethod(handles)

normmethod=getappdata(handles.leadfigure,'normmethod');
ndc=get(handles.normmethod,'String');

currentNormMethod=normmethod{get(handles.normmethod,'Value')};
try
[~,~,hassettings]=feval(currentNormMethod,'prompt');
catch % legacy support
    hassettings=0;
end
if hassettings
    set(handles.normsettings,'enable','on');
    setappdata(handles.normsettings,'currentNormMethod',currentNormMethod);
else
    set(handles.normsettings,'enable','off');
end





