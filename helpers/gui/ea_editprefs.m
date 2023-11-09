function ea_editprefs(varargin)

if ~isfile(ea_prefspath(ea_prefsext))
    copyfile([ea_getearoot,'common',filesep,'ea_prefs_default',ea_prefsext],ea_prefspath(ea_prefsext));
end
if ~isdeployed
    edit(ea_prefspath);
else

    msgbox('Open ~/.leaddbs/ea_prefs.json with a text editor to edit preferences.');

%     prefs = ea_prefs;
%     if ~isfield(prefs, 'textedit')
%         mb = msgbox('Choose app to open .json file');
%         uiwait(mb);
%         file = uigetfile('/Applications/*.app');
%         if file
%             ea_injectprefstring('textedit',file(1:end-4));
%             prefs = ea_prefs;
%         else
%             return
%         end
%     end
%     system(['open -a ','''',prefs.textedit,'''',' ',ea_prefspath('.json')])
end