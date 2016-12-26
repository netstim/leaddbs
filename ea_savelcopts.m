function ea_savelcopts(handles)

isindependent=getappdata(handles.leadfigure,'isindependent');

lc_options=ea_handles2lc(handles);
save([ea_getearoot,'connectomics',filesep,'lc_options.mat'],'-struct','lc_options');

if ~isindependent
    delete(handles.leadfigure);
    
    return
end