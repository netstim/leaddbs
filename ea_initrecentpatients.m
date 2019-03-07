function ea_initrecentpatients(handles,patsub)

earoot=ea_getearoot;
try
    load([earoot,'common',filesep,'ea_recent',patsub,'.mat']);
catch
    fullrpts={['No recent ',patsub,' found']};
end
save([earoot,'common',filesep,'ea_recent',patsub,'.mat'],'fullrpts');
ea_updaterecentpatients(handles,patsub);