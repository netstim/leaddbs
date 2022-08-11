function ea_openpatdir(handles)

outfolder=get(handles.patdir_choosebox,'String');

if ismember(outfolder,{'No Patient Selected','Choose Patient Directory'})
    msgbox('Please set the working directory first!', 'Error','error');
    return;
end

if startsWith(outfolder,'Multiple')
   msgbox('Multiple patients selected');
   return;
end

ea_opendir(outfolder);
