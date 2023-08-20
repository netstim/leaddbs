function ea_editprefs(varargin)

if ~exist([ea_gethome,'.ea_prefs',ea_prefsext],'file')
    copyfile([ea_getearoot,'common',filesep,'ea_prefs_default',ea_prefsext],[ea_gethome,'.ea_prefs',ea_prefsext]);
end

if ~isdeployed
    edit([ea_gethome,'.ea_prefs.m']);
else
    msgbox('Open ~/.ea_prefs.json with a text editor to edit preferences.');

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
    % system(['open -a ','''',prefs.textedit,'''',' ',fullfile(ea_gethome,'.ea_prefs.json')])
end