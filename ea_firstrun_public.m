function ea_firstrun(handles,options)

% check if a newer version is available..
local=ea_getvsn('local',1);
web=ea_getvsn('web',1);
vcheck=0;
if strcmp(local,'Unknown')
    vcheck=1;
elseif ~strcmp(web,'Unknown')
    vcheck=(local<web);
end

try
    if vcheck
        set(handles.updatebutn,'Visible','on');
        set(handles.updatebutn,'BackgroundColor',[0.2,0.8,0.2]);
    else
        set(handles.updatebutn,'Visible','off');
    end
end

if ~strcmp(handles.prod,'dbs_connectome')
    try
        webopts=weboptions('Timeout',5);
        webread('https://www.lead-dbs.org/release/stats.php','id',handles.prod,'ver',['R',version('-release')],webopts);
    catch
        try
            urlread(['https://www.lead-dbs.org/release/stats.php?id=',handles.prod,'&ver=R', version('-release')],'Timeout',5);
        catch
        end
    end
end

if ~isfield(options.prefs,'firstrun') % first run.
    fprintf(['Welcome to LEAD-DBS.\n \n',...
        'This seems to be your first run of the toolbox.\n\n',...
        'Check out the online manual at https://netstim.gitbook.io/leaddbs \n \n',...
        'Join our Slack community at https://leadusers.herokuapp.com \n \n',...
        'We hope that you enjoy using the Lead-DBS toolbox. \n \n',...
        'Any suggestions are more than welcome (andreas.horn@charite.de). \n'
        ]);

    if isdeployed % init .json file. This file has preferences that depend on directories
        fid = fopen([ea_getearoot,'common',filesep,'ea_prefs_default.json'],'wt');
        fwrite(fid, jsonencode(ea_prefs_default('')), 'char'); fclose(fid);
    end

    copyfile([ea_getearoot,'common',filesep,'ea_prefs_default', ea_prefsext],[ea_gethome,'.ea_prefs', ea_prefsext], 'f');

    ea_injectprefstring('firstrun','off');

    if ~isdeployed && ismac
        ea_clear_xattr;
    end

    % check dataset isntallation
    if ~exist([ea_space,'bb.nii'], 'file')
        fprintf(['\nIt seems that you don''t have LEAD dataset installed.\n' ...
                 'You can either install it via ''Install'' --> ''Redownload Data Files'' menu,\n' ...
                 'or download it from https://www.lead-dbs.org/release/download.php?id=data_pcloud or\n' ...
                 'https://www.lead-dbs.org/release/download.php?id=data_onedrive and then extract it into LEAD folder.\n\n']);

        msg = sprintf(['It seems that you don''t have LEAD dataset installed.\nDo you wish to download it now?\n' ...
                       'Alternatively, please check the command window for more information.']);
        choice = questdlg(msg, 'Download Dataset?', 'Yes', 'No', 'No');
        if strcmp(choice, 'Yes')
            ea_update_data('full');
        end
    end
end
