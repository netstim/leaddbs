function ea_switchnormmethod(handles,handlestring)

if ~exist('handlestring','var')
    handlestring='normmethod';
end

normmethod=getappdata(handles.leadfigure,'normmethod');
ndc=get(handles.(handlestring),'String');

currentNormMethod=normmethod{get(handles.(handlestring),'Value')};
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





