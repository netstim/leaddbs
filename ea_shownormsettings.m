function ea_shownormsettings(normfname,handles)

setname=['ea_normsettings_',normfname(14:end)];
feval(setname,handles);
