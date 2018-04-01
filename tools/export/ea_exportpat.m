function ea_exportpat(hObj, ~, exptype, handles)

uipatdir=getappdata(handles.leadfigure,'uipatdir');
if strcmp(uipatdir{1},'No Patient Selected')
    ea_error('Please select a patient first.');
    return
end

for pt=1:length(uipatdir)
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
        case 'LS'
            ea_exportpat([],[],'ZIP',handles);
            %ea_pat2ls(uipatdir{pt},handles);
        case 'ZIP'
            ea_pat2ply(uipatdir{pt},handles);
            ea_screenshots(uipatdir{pt},handles);
            [~,ptname]=fileparts(uipatdir{pt});
            zip([uipatdir{pt},filesep,'export',filesep,'zip',filesep,ptname,'.zip'],...
                {[uipatdir{pt},filesep,'export',filesep,'ply',filesep,'combined_scene.ply'],...
                [uipatdir{pt},filesep,'export',filesep,'views']},...
                [uipatdir{pt},filesep,'export',filesep]);
    end
end
