function ea_gencheckreg(~,~,handles)
% function to simply export checkreg figures for selected patients from
% Menu

disp('Generating check-reg figures for selected patients...');
leadfigure=handles.leadfigure;
ea_busyaction('on',leadfigure,'dbs');


options=ea_handles2options(handles);
options.macaquemodus=getappdata(handles.leadfigure,'macaquemodus');

options.uipatdirs=getappdata(handles.leadfigure,'uipatdir');
options.gencheckreg=1;

options.leadprod = 'dbs';
options.leadfigure=handles.leadfigure;

ea_run('run',options);

ea_busyaction('off',leadfigure,'dbs');
