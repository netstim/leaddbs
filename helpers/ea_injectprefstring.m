function ea_injectprefstring(str)

fid=fopen([ea_gethome, '.ea_prefs.m'],'a');
fprintf(fid,'%s\n',str);
fclose(fid);