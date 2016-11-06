function ea_openin(~,~,appname,handles)

feval(appname,'loadsubs',getappdata(handles.leadfigure,'uipatdir'));

