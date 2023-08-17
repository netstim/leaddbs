function ea_showprocessreport(~,~,handles)
exptxt='';
uipatdir=getappdata(handles.leadfigure,'uipatdir');
ea_dispercent(0,'Generating processing report');
for pt=1:length(uipatdir)
    [pth,ptname]=fileparts(uipatdir{pt});
    exptxt=[exptxt,'Subject/Patient: ',ptname,'\n'];
    exptxt=[exptxt,'--------------------------------------','\n','\n'];
    normmethod=ea_whichnormmethod(uipatdir{pt},'coregmrmethod');
    if ~isempty(normmethod)
        exptxt=[exptxt,'Postoperative MRIs were co-registered to preoperative MRIs using: \n',normmethod{end},'\n'];
    end
    
    normmethod=ea_whichnormmethod(uipatdir{pt},'coregctmethod');
    if ~isempty(normmethod)
        exptxt=[exptxt,'A postoperative CT was co-registered to preoperative MRIs using: \n',normmethod{end},'\n'];
    end
    
    
    normmethod=ea_whichnormmethod(uipatdir{pt});
    if ~isempty(normmethod)
        exptxt=[exptxt,'All acquisitions were warped into template space using: \n',normmethod,'\n'];
    end
    
    if exist([uipatdir{pt},filesep,'reconstruction',filesep,ptname,'_desc-reconstruction.mat'],'file')
        options=ea_getptopts(uipatdir{pt});
        options.sides=1:2;
        options.native=1;
        [~,~,~,~,manually_corrected]=ea_load_reconstruction(options);
        exptxt=[exptxt,'DBS electrodes were automatically (pre-)localized.\n'];
        if manually_corrected
            exptxt=[exptxt,'DBS electrodes were manually localized.\n'];
        end
    end
    exptxt=[exptxt,'\n\n\n'];
    ea_dispercent(pt/length(uipatdir));
end
ea_dispercent(1,'end');
ea_textdisp({exptxt});

