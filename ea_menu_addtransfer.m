function ea_menu_addtransfer(handles,prefs)

[leadsuiteapps,leadsuiteids,accs]=ea_getsuiteapps(prefs);
[~,todel]=ismember(handles.callingfunction,leadsuiteids);
leadsuiteids(todel)=[];
leadsuiteapps(todel)=[];
accs(todel)=[];

f = uimenu('Label','Open in...');
for app=1:length(leadsuiteids)
    uimenu(f,'Label',leadsuiteapps{app},'Callback',{@ea_openin,leadsuiteids{app},handles},'Accelerator',accs{app});
end
