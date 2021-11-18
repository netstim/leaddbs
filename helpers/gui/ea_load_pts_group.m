function ea_load_pts_group(handles,uipatdir)

ea_busyaction('on',handles.leadfigure,'group');

% check if Lead group is empty at present:
projectexists=0;
if ~strcmp(handles.groupdir_choosebox.String,'Choose Dataset Directory')
    answ=questdlg('A current project is already loaded in Lead-Group. Continuing will ADD the new patients to the current project. Please select Cancel to abort.','Add patients','Continue','Cancel','Cancel');
    switch lower(answ)
        case 'cancel'
            return
        case 'continue'
            projectexists=1;
    end
end

if ~projectexists
    a=questdlg('Please choose a (potentially blank) folder in which you would like to store your new Lead-Group project...','New Lead-Group project','Okay','Cancel','Cancel');
    if strcmpi(a,'cancel')
        return
    end

    % init M and add folder:
    nudir = uigetdir;
    
    if ~nudir % user pressed cancel
        return
    end
    
    nudir=[nudir,filesep];

    if ~isempty(ea_getGroupAnalysisFile(nudir))
        answ=questdlg('In the folder you selected, an old Lead-Group project was already found. Continuing will OVERWRITE the old project. Please press cancel to abort.','Existing Lead-Group project','Overwrite','Cancel','Cancel');
        switch lower(answ)
            case 'cancel'
                return
        end
    end

    set(handles.groupdir_choosebox,'String',nudir);

    % store blank M in folder.
    analysisFile = ea_genGroupAnalysisFile(nudir);

    load(analysisFile, 'M');
    setappdata(handles.leadfigure,'M',M);
    try
        setappdata(handles.leadfigure,'S',M.S);
        setappdata(handles.leadfigure,'vatmodel',M.S(1).model);
    end

    ea_refresh_lg(handles);
end

% add patients:
M = getappdata(handles.leadfigure,'M');
folders = GetFullPath(uipatdir); % folders to add to list.
M.patient.list=[M.patient.list;folders'];
M.patient.group=[M.patient.group;ones(length(folders),1)];

setappdata(handles.leadfigure,'M',M);
ea_refresh_lg(handles);

% save M
M=getappdata(handles.leadfigure,'M');
save(ea_getGroupAnalysisFile(handles.groupdir_choosebox.String),'M','-v7.3');

ea_busyaction('off',handles.leadfigure,'group');
