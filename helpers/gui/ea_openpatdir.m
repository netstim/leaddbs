function ea_openpatdir(handles)

BIDSRoot = handles.datasetselect.String;
%get(handles.datasetselect.String,'String');
selectedSubj = strcat('sub-',handles.patientlist.Data.subjId(handles.patientlist.Selection));
if ismember(BIDSRoot,{'No Patient Selected','Choose Dataset Directory'})
    msgbox('Please set the working directory first!', 'Error','error');
    return;
end

if length(selectedSubj) > 1
   msgbox('Multiple patients selected! Please select one patient.');
   return;
end
outfolder = fullfile(BIDSRoot,'derivatives','leaddbs',selectedSubj{1});
ea_opendir(outfolder);
