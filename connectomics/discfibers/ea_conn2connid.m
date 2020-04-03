function conname=ea_conn2connid(conname)

conname=strrep(conname,' ','');
conname=strrep(conname,'_','');
conname=strrep(conname,'(','');
conname=strrep(conname,')','');
conname=strrep(conname,'-','');
end