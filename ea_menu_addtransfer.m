function ea_menu_addtransfer(handles,callingfunction)

[leadsuiteapps,leadsuiteids]=ea_getsuiteapps;
[~,todel]=ismember(callingfunction,leadsuiteids);
leadsuiteids(todel)=[];
leadsuiteapps(todel)=[];

f = uimenu('Label','Open in...');
for app=1:length(leadsuiteids)
    uimenu(f,'Label',leadsuiteapps{app},'Callback',{@ea_openin,leadsuiteids{app},handles});
end
