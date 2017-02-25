function ea_aggregateVTA(~,~,handles,exportwhat)


if ~strcmp(handles.seeddefpopup.String{handles.seeddefpopup.Value}(1:9),'Use VATs:')
    ea_error('Please select a stimulation to export first.');
else
    stimname=handles.seeddefpopup.String{handles.seeddefpopup.Value}(11:end);
end
uipatdir=getappdata(handles.leadfigure,'uipatdir');
fname=uigetdir('','Select where to save files...');
fname=[fname,filesep];

tf=fopen([fname,'export.txt'],'w');
for pt=1:length(uipatdir)
    [pth,ptname]=fileparts(uipatdir{pt});
    copyfile([uipatdir{pt},filesep,'stimulations',filesep,stimname,filesep,exportwhat,'.nii'],[fname,ptname,'.nii']);
fprintf(tf,'%s\n',[fname,ptname,'.nii']);
end
fclose(tf);


msgbox('Files successfully exported.');