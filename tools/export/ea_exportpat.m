function ea_exportpat(hObj,~,exptype,handles)

uipatdir=getappdata(handles.leadfigure,'uipatdir');
if strcmp(uipatdir{1},'No Patient Selected')
    ea_error('Please select a patient first.');
    return
end


for pt=1:length(uipatdir)
    if ~exist([uipatdir{pt},filesep,'export'],'dir')
    mkdir([uipatdir{pt},filesep,'export']);
    end
    switch exptype
        case 'PDF'
            mkdir([uipatdir{pt},filesep,'export',filesep,'pdf']);
            ea_pat2pdf(uipatdir{pt},handles);
        case 'STL'
            mkdir([uipatdir{pt},filesep,'export',filesep,'stl']);
            ea_pat2stl(uipatdir{pt},handles);
         case 'PLY'
            mkdir([uipatdir{pt},filesep,'export',filesep,'ply']);
            ea_pat2ply(uipatdir{pt},handles);
    end
    
end