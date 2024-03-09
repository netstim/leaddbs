function ea_archivefolders(hobj,evt,cmd,handles)



options = ea_handles2options(handles);
options.uipatdirs = getappdata(handles.leadfigure, 'uipatdir');

options.leadprod = 'dbs';

setappdata(handles.leadfigure,'handles', handles);
options.leadfigure = handles.leadfigure;

for pt=1:length(options.uipatdirs)
    disp(['Working on pt ',options.uipatdirs{pt},'...']);
    options=ea_getptopts(options.uipatdirs{pt},options);

    switch cmd
        case 'archive'
            ea_archivefolder(options);
        case 'unarchive'
            ea_unarchivefolder(options);

    end

end







