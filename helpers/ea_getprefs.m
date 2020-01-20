function val=ea_getprefs(itemname)

prefs=ea_prefs;
machine=prefs.machine;
val=eval(['machine.', itemname, ';']);
