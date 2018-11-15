function ea_exportpat(hObj, ~, exptype, handles, target)

uipatdir=getappdata(handles.leadfigure,'uipatdir');
if strcmp(uipatdir{1},'No Patient Selected')
    ea_error('Please select a patient first.');
    return
end
disp('Exporting Results...');
for pt=1:length(uipatdir)
    try
        if ~exist([uipatdir{pt}, filesep, 'export', filesep, lower(exptype)],'dir')
            mkdir([uipatdir{pt}, filesep, 'export', filesep, lower(exptype)]);
        end
        switch exptype
            case 'PDF'
                ea_pat2pdf(uipatdir{pt},handles);
            case 'STL'
                ea_pat2stl(uipatdir{pt},handles);
            case 'PLY'
                ea_pat2ply(uipatdir{pt},handles);
            case 'ZIP'
                ea_pat2ply(uipatdir{pt},handles,target);
                ea_screenshots(uipatdir{pt},handles,target);
                [~,ptname]=fileparts(uipatdir{pt});
                zip([uipatdir{pt},filesep,'export',filesep,'zip',filesep,ptname,'.zip'],...
                    {[uipatdir{pt},filesep,'export',filesep,'ply',filesep,'anatomy.ply'],...
                    [uipatdir{pt},filesep,'export',filesep,'ply',filesep,'combined_electrodes.ply'],...
                    [uipatdir{pt},filesep,'export',filesep,'views']},...
                    [uipatdir{pt},filesep,'export',filesep]);
            case 'LS'
                ea_exportpat([],[],'ZIP',handles,target);
                response=ea_upload_export(uipatdir{pt});
                if strcmp(response.StatusCode,'OK')
                    disp('*** Upload successful.');
                else
                    warndlg('Upload failed.','Upload error.');
                end
        end
    catch
        msgbox(['Export for pt: ',uipatdir{pt},' failed.']);
    end
end
