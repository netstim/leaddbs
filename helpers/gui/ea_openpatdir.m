function ea_openpatdir(handles)

selectedFolder = getappdata(handles.leadfigure, 'uipatdir');
if isempty(selectedFolder)
    msgbox('Please choose the subject directory first!', '','error');
    return;
elseif length(selectedFolder) > 1
    msgbox('Select only one subject to open the folder.', '', 'warn');
else
    ea_opendir(selectedFolder{1});
end
