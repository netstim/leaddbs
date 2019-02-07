function ea_deletefiducialhelpers(~,~,handles)
answ=questdlg('Are you sure you want to delete fiducial marker(s) for selected patient(s)?','Delete Fiducial Markers','Yes','Abort','Abort');
if strcmp(answ,'Abort')
    return
end


uipatdir=getappdata(handles.leadfigure,'uipatdir');
cnt=1;
for pt=1:length(uipatdir)
    fids=dir([uipatdir{pt},filesep,'fiducials',filesep,'native',filesep,'*.nii.gz']);
    
    for fi=1:length(fids)
        % pt folder
        ea_delete([uipatdir{pt},filesep,'fiducials',filesep,'native',filesep,fids(fi).name]);
        
        % templates
        ea_delete([uipatdir{pt},filesep,'fiducials',filesep,ea_getspace,filesep,fids(fi).name]);
    end
    
end

