function ea_placeIXIages_into_db

load([ea_getearoot,'toolbox',filesep,'IXI',filesep,'IXI_demographics'])
prefs=ea_prefs('');
for i=1:length(IXI.uID)
    thisdir=[prefs.ixi.dir,'IXI',sprintf('%03.f',IXI.uID(i)),filesep];
    if exist(thisdir,'file')
    f=fopen([thisdir,'ea_age.txt'],'w');
    fprintf(f,'%d',IXI.uAGE(i));
    fclose(f);
    end
end