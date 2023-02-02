function ea_exportpat(~, ~, exptype, handles, target)

uipatdir=getappdata(handles.leadfigure,'uipatdir');
if strcmp(uipatdir{1},'No Patient Selected')
    ea_error('Please select a patient first.');
    return
end
disp('Exporting Results...');
for pt=1:length(uipatdir)
    %try
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
                ea_screenshots(uipatdir{pt},target);
                [~,ptname]=fileparts(uipatdir{pt});
                zip([uipatdir{pt},filesep,'export',filesep,'zip',filesep,ptname,'.zip'],...
                    {[uipatdir{pt},filesep,'export',filesep,'ply',filesep,'anatomy.ply'],...
                    [uipatdir{pt},filesep,'export',filesep,'ply',filesep,'combined_electrodes.ply'],...
                    [uipatdir{pt},filesep,'export',filesep,'views']},...
                    [uipatdir{pt},filesep,'export',filesep]);
            case 'LS'
                if ~exist([uipatdir{pt},filesep,'ea_pseudonym.mat'],'file')
                    [pth,ptname]=fileparts(uipatdir{pt});
                    patientPseudonym = inputdlg('Please enter patient synonym','Patient Synonym',1,{ptname});
                    save([uipatdir{pt},filesep,'ea_pseudonym.mat'],'patientPseudonym');
                else % if patient has already been exported, we already know the pseudonym.
                    p=load([uipatdir{pt},filesep,'ea_pseudonym.mat']);
                    patientPseudonym=p.patientPseudonym; clear p
                end
                prefs=ea_prefs;
                usercredentials = inputdlg({'Username','Password'},'Enter User Credentials',[1 35],{prefs.tbase.user,prefs.tbase.pw});

                ea_exportpat([],[],'ZIP',handles,target);
                response=ea_upload_export(uipatdir{pt},patientPseudonym,usercredentials);
                if strcmp(response.StatusCode,'OK')
                    disp('*** Upload successful.');
                else
                    warndlg('Upload failed.','Upload error.');
                end
        end
    %catch
    %    keyboard
    %    msgbox(['Export for pt: ',uipatdir{pt},' failed.']);
    %end
end
