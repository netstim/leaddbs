function ea_editprefs(varargin)

prefsPath = ea_prefspath(ea_prefsext);

if ~isfile(prefsPath)
    copyfile([ea_getearoot,'common',filesep,'ea_prefs_default',ea_prefsext], prefsPath);
end

if ~isdeployed
    edit(prefsPath);
else
    msgbox('Open ~/.leaddbs/ea_prefs_user.json with a text editor to edit preferences.');

    % prefs = ea_prefs;
    % if ~isfield(prefs, 'textedit')
    %     mb = msgbox('Choose app to open .json file');
    %     uiwait(mb);
    %     file = uigetfile('/Applications/*.app');
    %     if file
    %         ea_setprefs('textedit', file(1:end-4), 'user');
    %         prefs = ea_prefs;
    %     else
    %         return
    %     end
    % end
    % system(['open -a ','''',prefs.textedit,'''',' ',prefsPath])
end
