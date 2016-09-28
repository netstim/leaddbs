function sfile=ea_getrois(sfile)

fID=fopen(sfile);
sfile=textscan(fID,'%s');
sfile=sfile{1};
fclose(fID);
