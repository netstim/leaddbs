function M=ea_checkremap_lg(M,handles)

if ~exist(M.patient.list{1},'dir')

    [rootpth,patientname]=fileparts(M.patient.list{1});
    defdirloc=fullfile(fileparts(fileparts(fileparts((M.root)))),'leaddbs');
    if exist(fullfile(defdirloc,patientname),'dir')
        answ=questdlg({'Patients folders do not reside within the specified location: ',...
            rootpth,...
            'but seem to reside in the default BIDS location: ',...
            defdirloc,...
            'Do you want to remap lead group pointers to the default location?'},'Remap folders','Yes','No','Yes');
        if strcmp(answ,'Yes')
            % Adapt patient base folder ('dataset/derivatives/leaddbs' folder)

            newBaseFolder = defdirloc;

            for i=1:length(M.patient.list)
                [~, ptsName] = fileparts(M.patient.list{i});
                M.patient.list{i} = fullfile(newBaseFolder, ptsName);
            end
            setappdata(handles.leadfigure, 'M', M);

            analysisFile = ea_regexpdir(M.root, ['^dataset-[^\W_]+_analysis-', M.guid, '\.mat$'], 0);
            save(analysisFile{1}, 'M', '-v7.3');
        end
    end
end

