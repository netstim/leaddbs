function ea_setprefs(itemname,value)

prefs=ea_prefs;
machine=prefs.machine;
machine.(itemname)=value;
try % may not have write permissions
    save([ea_gethome,'.ea_prefs.mat'],'machine');
catch
    warning('Could not save preferences to user home directory. Please assign permissions correctly.');
end